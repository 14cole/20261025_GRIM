"""Reproducible spatial cuts fitted to original complex radar samples.

Arrays use (cross-range, down-range) order, independent of display resolution.
Native fitting avoids trying to undo tapered, interpolated display images.
Include is a spatial forward projection; Exclude retains the fit residual.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from copy import deepcopy
import json
import os
import tempfile

import numpy as np
from matplotlib.path import Path

from GRIM_Backend.datasets.result_storage import allocate_result, operation_result_storage, result_blocks
from GRIM_Backend.plotting.requests import PlotCancelled


RECONSTRUCTION_NOTE = (
    "Native-sample spatial decomposition: a minimum-norm complex least-squares image is fitted "
    "to the original samples before masking. Include is the exact forward projection of the "
    "selected coefficients; Exclude is original minus Include and retains the fit residual. "
    "A full-scene mask is an exact pass-through. Missing samples stay unknown. "
    "The fit is checked on all finite native samples (relative complex error <=0.1%). "
    "A good data fit does not establish unique scattering attribution or calibrated component RCS."
)


@dataclass
class MaskImage:
    name: str
    x: np.ndarray
    y: np.ndarray
    magnitude: np.ndarray
    complex_image: np.ndarray | None = None
    inverse: dict | None = None
    metadata: dict = field(default_factory=dict)

    @property
    def reconstruction_reason(self):
        if self.inverse is not None and self.complex_image is not None:
            return RECONSTRUCTION_NOTE
        if self.metadata.get("reconstruction") == "backprojection":
            return "Replot ISAR to retain the native backprojection source needed for radar reconstruction."
        if self.metadata.get("kind") == "composite":
            return "Image dataset only: Composite combines magnitudes and has no recoverable complex phase."
        return "Image dataset only: radar reconstruction requires a retained coherent ISAR source."


@dataclass
class MaskImageRequest:
    """Captured formation recipe for a lightweight playback frame."""
    name: str
    params: dict
    display: dict = field(default_factory=dict)
    source_digest: str | None = None


def prepare_mask_source(request, cancel=None):
    from GRIM_Backend.plotting.modes.isar_mode import compute_bands
    options = dict(request.params, retain_complex=request.params.get('aperture_mode') != 'composite',
                   retain_mask_source=True, decimate_display=False, cancel_check=cancel)
    result = compute_bands(options)
    _check(cancel)
    if isinstance(result, str):
        raise ValueError(result)
    bands, _elapsed = result
    if request.source_digest is not None and bands[0].get('source_selection_digest') != request.source_digest:
        raise ValueError("The source data changed after this frame was plotted. Replot before making a cut.")
    image = image_from_band(bands[0], options)
    image.metadata['display'] = dict(request.display)
    return image


def reconstruction_extra(dataset):
    """Retain coordinate/field lineage without claiming solver certification."""
    extra = dataset._safe_derived_scalar_extra(include_field_conventions=True)
    for key in ("phase_center_xyz_m", "phase_reference", "time_convention", "polarization_basis",
                "range_phase_convention", "measurement_geometry", "measurement_domain",
                "motion_compensation", "phase_center_motion"):
        if key in (dataset.extra or {}):
            extra[key] = deepcopy(dataset.extra[key])
    dataset._invalidate_assembly_sampling_hash(extra, "ISAR Mask Cut reconstruction")
    return extra


def image_from_band(band, params):
    kind = "composite" if band.get("composite") else "isar"
    az = band["az_values"]
    inverse = band.get("mask_inverse") or {}
    dataset = params.get("dataset")
    polarization = (str(dataset.polarizations[params["pol_idx"]])
                    if dataset is not None and "pol_idx" in params else None)
    name = f"{'Composite' if kind == 'composite' else 'ISAR'} {az[0]:g}°–{az[-1]:g}°"
    return MaskImage(name, np.asarray(band["x_range"]).copy(), np.asarray(band["y_range"]).copy(),
                     band.get("mask_magnitude", band["magnitude"]), band.get("complex_image"),
                     band.get("mask_inverse"),
                     {"kind": kind, "unit": params.get("unit_name", "m"),
                      "contract": deepcopy(band.get("image_contract", {})),
                      "source_name": str(getattr(dataset, "source_path", "") or ""),
                      "source_dataset_id": params.get("source_dataset_id"),
                      "source_dataset_name": params.get("source_dataset_name"),
                      "source_selection_digest": band.get("source_selection_digest"),
                      "frequency_hz": np.asarray(params.get("freq_hz", inverse.get("native_frequency_hz", []))).tolist(),
                      "azimuth_degrees": np.asarray(az).tolist(),
                      "elevation_degrees": params.get("elevation_deg"), "polarization": polarization,
                      "window": params.get("window_name", inverse.get("window")),
                      "reconstruction": band.get("resolved_reconstruction"),
                      "notes": list(band.get("isar_metadata_notes", ())), "cuts": []})


def _check(cancel):
    if cancel is not None and cancel():
        raise PlotCancelled()


def mask_preview(image, cancel=None, max_side=1200):
    """Bounded display with exact block peaks; cuts still use native centers."""
    from .geometry import pool_image_peaks, reduced_axis_edges
    sx, sy = (max(1, int(np.ceil(size / max_side))) for size in image.magnitude.shape)
    if sx == sy == 1:
        values = np.asarray(image.magnitude, dtype=np.float32)
    else:
        blocks = [(start, min(start + sx, len(image.x))) for start in range(0, len(image.x), sx)]
        values = pool_image_peaks(image.magnitude.T, sy, blocks, cancel=lambda: _check(cancel))[0].T
    display = image.metadata.get('display', {})
    level = np.square(values, dtype=np.float32)
    linear = display.get('linear', False)
    if not linear:
        np.maximum(level, np.float32(1e-12), out=level)
        np.log10(level, out=level)
        level *= 10
    return dict(values=level, x=reduced_axis_edges(image.x, max_side), y=reduced_axis_edges(image.y, max_side),
                linear=linear, cmap=display.get('cmap', 'viridis'), limits=display.get('limits'))


def cut_image(source, vertices, *, include, name, cancel=None):
    """Zero the discarded pixels on the original coordinate grid, preserving NaNs."""
    vertices = np.asarray(vertices, dtype=float)
    if vertices.ndim != 2 or vertices.shape[1] != 2 or len(vertices) < 3 or not np.isfinite(vertices).all():
        raise ValueError("Draw a rectangle or polygon with at least three finite vertices.")
    area = abs(np.sum(vertices[:, 0] * np.roll(vertices[:, 1], -1)
                      - vertices[:, 1] * np.roll(vertices[:, 0], -1)))
    if area <= 0:
        raise ValueError("The selected area is empty.")
    shape = (len(source.x), len(source.y))
    if source.magnitude.shape != shape:
        raise ValueError("This plot has no native mask source. Replot it before using Mask Cut.")
    # Positive radius includes the drawn boundary for counter-clockwise paths.
    signed = np.sum(vertices[:, 0] * np.roll(vertices[:, 1], -1)
                    - vertices[:, 1] * np.roll(vertices[:, 0], -1))
    if signed < 0:
        vertices = vertices[::-1]
    polygon = Path(np.vstack((vertices, vertices[0])))
    from .mask_geometry import BOUNDARY_TOLERANCE_METRES
    # Fixed physical tolerance keeps boundary membership invariant under units.
    scale = float(source.metadata.get('contract', {}).get('distance_scale_per_metre', 1.))
    tolerance = 2 * BOUNDARY_TOLERANCE_METRES * scale
    with operation_result_storage("auto"):
        magnitude = allocate_result(shape, source.magnitude.dtype)
        complex_image = (allocate_result(shape, source.complex_image.dtype)
                         if source.complex_image is not None else None)
        selected = 0
        for block in result_blocks(shape):
            _check(cancel)
            x, y = np.meshgrid(source.x[block[0]], source.y[block[1]], indexing="ij")
            inside = polygon.contains_points(np.column_stack((x.ravel(), y.ravel())),
                                             radius=tolerance).reshape(x.shape)
            selected += int(inside.sum())
            keep = inside if include else ~inside
            tile = source.magnitude[block]
            magnitude[block] = np.where(keep | ~np.isfinite(tile), tile, 0.)
            if complex_image is not None:
                tile = source.complex_image[block]
                complex_image[block] = np.where(keep | ~np.isfinite(tile), tile, 0.)
        if not selected:
            raise ValueError("The selected area contains no image pixel centers. Draw a larger area.")
    metadata = deepcopy(source.metadata)
    metadata.setdefault("cuts", []).append({"operation": "include" if include else "exclude",
                                            "vertices": vertices.tolist(), "source": source.name})
    return MaskImage(name.strip() or "Mask cut", source.x.copy(), source.y.copy(), magnitude,
                     complex_image, source.inverse, metadata)


def _unflipped_image(image, cancel=None):
    """Restore a cropped/mirrored image to its full shifted Fourier grid."""
    context = image.inverse
    full = allocate_result((len(context["full_x"]), len(context["full_y"])), np.complex64, fill=0)
    signs = image.metadata.get("contract", {}).get("coordinate_signs", [1, 1])
    x, y = signs[0] * image.x, signs[1] * image.y
    ix = np.rint((x - context["full_x"][0]) / np.diff(context["full_x"])[0]).astype(int)
    iy = np.rint((y - context["full_y"][0]) / np.diff(context["full_y"])[0]).astype(int)
    if np.any(ix < 0) or np.any(ix >= full.shape[0]) or np.any(iy < 0) or np.any(iy >= full.shape[1]):
        raise ValueError("Image coordinates do not match the retained reconstruction grid.")
    for block in result_blocks(image.complex_image.shape):
        _check(cancel)
        full[np.ix_(ix[block[0]], iy[block[1]])] = image.complex_image[block]
    return full


def inverse_spectrum(image, cancel=None):
    """Invert the normalized, padded IFFT in bounded batches of 1-D FFTs."""
    if image.inverse is None or image.complex_image is None:
        raise ValueError(image.reconstruction_reason)
    try:
        from scipy.fft import fft
    except ImportError:  # NumPy may use complex128 scratch; batches stay bounded.
        from numpy.fft import fft
    context = image.inverse
    full = _unflipped_image(image, cancel)
    na, nf = full.shape
    nv, nu = len(context["v"]), len(context["u"])
    stage = allocate_result((na, nv), np.complex64)
    rows = max(1, 8 * 1024**2 // (32 * nf))
    for start in range(0, na, rows):
        _check(cancel)
        tile = full[start:start + rows]
        stage[start:start + rows] = fft(np.fft.ifftshift(tile, axes=1), axis=1)[:, :nv]
    del full
    result = allocate_result((nu, nv), np.complex64)
    cols = max(1, 8 * 1024**2 // (32 * na))
    scale = context["gain"] / (na * nf)
    for start in range(0, nv, cols):
        _check(cancel)
        tile = stage[:, start:start + cols]
        result[:, start:start + cols] = fft(np.fft.ifftshift(tile, axes=0), axis=0)[:nu] * scale
    return result


def _project_backprojection_image(image, cancel=None):
    """Exact separable point sum on native looks, with bounded phase matrices."""
    from .quality import C0
    context = image.inverse
    signs = image.metadata.get('contract', {}).get('coordinate_signs', [1, 1])
    scale = context['distance_scale_per_metre']
    x, y = signs[0] * image.x / scale, signs[1] * image.y / scale
    theta = np.deg2rad(context['native_azimuth_degrees'] - context['look_center_degrees'])
    frequencies = context['native_frequency_hz']
    factor = 4 * np.pi / C0 * abs(np.cos(np.deg2rad(context['elevation_degrees'])))
    output = allocate_result((len(theta), len(frequencies)), np.complex128)
    for row, angle in enumerate(theta):
        for start in range(0, len(frequencies), 64):
            _check(cancel)
            k = factor * frequencies[start:start + 64]
            summed = np.zeros((len(k), len(y)), np.complex128)
            for first in range(0, len(x), 128):
                _check(cancel)
                phase_x = np.exp(-1j * np.sin(angle) * k[:, None] * x[None, first:first + 128])
                summed += phase_x @ image.complex_image[first:first + 128]
            phase_y = np.exp(-1j * np.cos(angle) * k[:, None] * y[None, :])
            values = np.sum(summed * phase_y, axis=1)
            output[row, start:start + len(k)] = values.conj() if context['conjugate_time'] else values
    return output


def _prepare_native_cut(source, cancel=None):
    """Invert the native forward model, then independently verify every sample."""
    from scipy.sparse.linalg import lsmr
    from .mask_operator import MaskPolarOperator
    context = source.inverse
    if context is None or 'native_samples' not in context:
        raise ValueError('Replot ISAR to retain the native samples required for a reproducible cut.')
    if 'native_fit' in context:
        return source
    measured = np.asarray(context['native_samples'], complex)
    valid = np.isfinite(measured)
    signs = source.metadata.get('contract', {}).get('coordinate_signs', [1, 1])
    scale = context['distance_scale_per_metre']
    x, y = signs[0]*source.x/scale, signs[1]*source.y/scale
    if measured.size * len(x) * len(y) > 2_000_000_000:
        raise ValueError('Mask reconstruction exceeds its verification budget. Select a smaller aperture, band, or scene.')
    op = MaskPolarOperator(x, y, context['native_azimuth_degrees']-context['look_center_degrees'],
                          context['native_frequency_hz'], context['elevation_degrees'], valid=valid, cancel=cancel)
    observed = measured[valid]
    amplitude = float(np.max(abs(observed), initial=0.))
    if context['conjugate_time']:
        observed = observed.conj()
    rhs = observed / (amplitude or 1.) / op.normalization
    solved = lsmr(op, rhs, atol=1e-8, btol=1e-4, conlim=1e6, maxiter=300)
    _check(cancel)
    coefficients = solved[0].reshape(source.magnitude.shape)*amplitude
    metadata = deepcopy(source.metadata)
    metadata['contract'].update(spatial_frequency_origin_hz=[0., 0.], spatial_frequency_step_hz=None,
                                amplitude_normalization='native_sample_fit_point_coefficients')
    fitted = replace(source, complex_image=coefficients, magnitude=abs(coefficients), metadata=metadata)
    predicted = _project_backprojection_image(fitted, cancel)
    residual = float(np.linalg.norm((predicted[valid]-measured[valid])/(amplitude or 1.)) /
                     max(np.linalg.norm(measured[valid]/(amplitude or 1.)), np.finfo(float).tiny))
    if not np.isfinite(residual) or residual > .001 or solved[1] in (3, 6):
        raise ValueError(f'The selected scene cannot reproduce the native data accurately (complex fit error {residual:.3%}). '
                         'No cut dataset was created. Try a narrower aperture or adjust the scene bounds.')
    diagnostic = dict(method='native_polar_lsmr_v1', operator='oversampled_kaiser_sinc_16_tap',
                      verification='exact_native_polar_sum_all_finite_samples',
                      relative_complex_residual=residual, maximum_relative_residual=.001,
                      iterations=int(solved[2]), stop_code=int(solved[1]), condition_estimate=float(solved[6]),
                      atol=1e-8, btol=1e-4, maxiter=300,
                      residual_policy='original_minus_include_retained_in_exclude')
    return replace(fitted, inverse=dict(context, native_fit=diagnostic))


def _radar_dataset(image, output):
    from GRIM_Backend.datasets.grid import RcsGrid
    context = image.inverse
    note = image.reconstruction_reason
    extra = dict(context['extra'])
    details = dict(version=1, image=image.metadata, reconstruction=note, window=context['window'])
    if 'native_fit' in context:
        details['native_fit'] = context['native_fit']
    extra['mask_cut_json'] = json.dumps(details, default=_json_value)
    extra['data_origin'] = 'ISAR mask reconstructed phase history'
    extra.pop('source_path', None)
    history = context['history'] + '\nMask Cut: ' + image.name + '. ' + note
    return RcsGrid(context['azimuths'], context['elevations'], context['frequencies'],
                   context['polarizations'], rcs=output[:, None, :, None],
                   units=dict(context['units']), extra=extra, history=history)


def reconstruct_radar(image, cancel=None):
    """Return an RcsGrid on the selected original azimuth/frequency axes."""
    with operation_result_storage("auto"):
        if image.inverse is not None and 'native_samples' in image.inverse:
            _check(cancel)
            context = image.inverse
            if image.metadata.get('cuts') and 'native_cut_samples' not in context:
                raise ValueError('Use Finish Cut to fit the native data before exporting a spatial mask.')
            output = np.array(context.get('native_cut_samples', context['native_samples']), copy=True)
            _check(cancel)
            return _radar_dataset(image, output)
        spectrum = inverse_spectrum(image, cancel)
        context = image.inverse
        az = np.asarray(context["native_azimuth_degrees"])
        frequencies = np.asarray(context["native_frequency_hz"])
        u, v = np.asarray(context["u"]), np.asarray(context["v"])
        output = allocate_result((len(az), len(frequencies)), np.complex64)
        for block in result_blocks(output.shape):
            _check(cancel)
            psi = np.deg2rad(az[block[0]] - context["look_center_degrees"])[:, None]
            f = frequencies[block[1]][None, :]
            qu = f * np.sin(psi)
            qv = f * (np.cos(psi) if context["accurate"] else np.ones_like(psi))
            a = (qu - u[0]) / (u[-1] - u[0]) * (len(u) - 1)
            b = (qv - v[0]) / (v[-1] - v[0]) * (len(v) - 1)
            valid = ((a >= -1e-7) & (a <= len(u)-1+1e-7)
                     & (b >= -1e-7) & (b <= len(v)-1+1e-7)
                     & context["native_valid"][block])
            a, b = np.clip(a, 0, len(u)-1), np.clip(b, 0, len(v)-1)
            ia = np.minimum(a.astype(int), len(u)-2)
            ib = np.minimum(b.astype(int), len(v)-2)
            da, db = a - ia, b - ib
            tile = ((1-da) * ((1-db)*spectrum[ia, ib] + db*spectrum[ia, ib+1])
                    + da * ((1-db)*spectrum[ia+1, ib] + db*spectrum[ia+1, ib+1]))
            if context["conjugate_time"]:
                tile = tile.conj()
            output[block] = np.where(valid, tile, complex(np.nan, np.nan))
        return _radar_dataset(image, output)


def _cut_results(source, vertices, *, name, cancel=None, mask_definition=None):
    from .mask_geometry import define_mask, with_operation, validate_definition
    if source.inverse is None or 'native_samples' not in source.inverse:
        raise ValueError('Replot ISAR before cutting: the original native samples are required.')
    if source.metadata.get('cuts'):
        raise ValueError('Plot the cut dataset in ISAR before applying another mask.')
    # A unity image tests geometry independently of zero or missing image values.
    membership = cut_image(replace(source, magnitude=np.ones(source.magnitude.shape), complex_image=None),
                           vertices, include=True, name=name, cancel=cancel).magnitude != 0
    measured = np.asarray(source.inverse['native_samples'], complex)
    valid = np.isfinite(measured)
    definition = (validate_definition(source, vertices, mask_definition) if mask_definition is not None
                  else define_mask(source, vertices))
    if membership.all():
        fitted = replace(source, inverse=dict(source.inverse, native_fit=dict(
            method='full_scene_identity_v1', relative_complex_residual=0.,
            residual_policy='original_minus_include_retained_in_exclude')))
        inside = np.array(measured, copy=True)
    else:
        fitted = _prepare_native_cut(source, cancel)
        selected = replace(fitted, complex_image=np.where(membership, fitted.complex_image, 0.))
        inside = _project_backprojection_image(selected, cancel)
    inside[~valid] = complex(np.nan, np.nan)
    outside = measured - inside
    outside[~valid] = complex(np.nan, np.nan)
    results = []
    for include, suffix, samples in ((True, 'include', inside), (False, 'exclude', outside)):
        _check(cancel)
        result = cut_image(fitted, vertices, include=include, name=f'{name}_{suffix}', cancel=cancel)
        result.metadata['mask_definition'] = deepcopy(with_operation(definition, include))
        result.inverse = dict(result.inverse, native_cut_samples=samples)
        results.append((result, _radar_dataset(result, samples)))
    return results


def finish_cut(source, vertices, *, include, name, cancel=None, mask_definition=None):
    results = _cut_results(source, vertices, name=name, cancel=cancel, mask_definition=mask_definition)
    image, radar = results[0 if include else 1]
    image.name = name
    # Use the final name in provenance as well as in the dataset table.
    return image, _radar_dataset(image, image.inverse['native_cut_samples'])


def comparison_cuts(source, vertices, *, name, cancel=None, mask_definition=None):
    """One fitted image and boundary; Include + Exclude equals the native data."""
    return _cut_results(source, vertices, name=name, cancel=cancel, mask_definition=mask_definition)


def _json_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Unsupported image metadata: {type(value).__name__}")


def save_image(path, image):
    """Atomic NPZ artifact; numeric arrays and JSON only, no pickle."""
    arrays = {"x": image.x, "y": image.y, "magnitude": image.magnitude}
    if image.complex_image is not None:
        arrays["complex_image"] = image.complex_image
    context = None
    if image.inverse is not None:
        context = {}
        for key, value in image.inverse.items():
            if isinstance(value, np.ndarray):
                arrays["inverse_" + key] = value
            else:
                context[key] = value
    arrays["manifest"] = np.asarray(json.dumps({"format": "GRIM-mask-image", "version": 1,
        "name": image.name, "metadata": image.metadata, "inverse": context}, default=_json_value))
    folder = os.path.dirname(os.path.abspath(path))
    descriptor, temporary = tempfile.mkstemp(prefix=".grim-mask-", suffix=".npz", dir=folder)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            np.savez_compressed(stream, **arrays)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def load_image(path, cancel=None, *, maximum_bytes=16 * 1024**3):
    from .mask_storage import read_mask_archive
    manifest, arrays = read_mask_archive(path, cancel=cancel, maximum_bytes=maximum_bytes)
    with operation_result_storage('auto'):
        x, y, magnitude = (arrays[key] for key in ("x", "y", "magnitude"))
        if (x.ndim != 1 or y.ndim != 1 or min(len(x), len(y)) < 2
                or not np.isfinite(x).all() or not np.isfinite(y).all()
                or not np.all(np.diff(x) > 0) or not np.all(np.diff(y) > 0)
                or magnitude.shape != (len(x), len(y)) or np.iscomplexobj(magnitude)):
            raise ValueError("Invalid image axes or magnitude shape.")
        for block in result_blocks(magnitude.shape):
            _check(cancel)
            if np.isinf(magnitude[block]).any() or np.any(magnitude[block] < 0):
                raise ValueError("Image magnitude must be nonnegative and cannot contain infinity.")
        complex_image = arrays.get("complex_image")
        if complex_image is not None and (complex_image.shape != magnitude.shape or not np.iscomplexobj(complex_image)):
            raise ValueError("Invalid complex image shape or type.")
        inverse = manifest["inverse"]
        if inverse is not None:
            inverse.update({key[len("inverse_"):]: value for key, value in arrays.items() if key.startswith("inverse_")})
        return MaskImage(manifest["name"], x, y, magnitude, complex_image, inverse, manifest["metadata"])
