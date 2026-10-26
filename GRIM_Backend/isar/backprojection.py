"""Coherent far-field backprojection on native polar measurements.

No angular interpolation or small-angle approximation is used. Uniform
frequency rows use oversampled range profiles; irregular rows use a bounded
direct Fourier sum. Both retain the carrier phase and missing-sample mask.
"""
from __future__ import annotations

import numpy as np
from scipy.fft import ifft, next_fast_len

from .quality import C0


def geometry(azimuth_degrees, frequency_hz, elevation_degrees, half=None, center=None):
    az, freq = np.asarray(azimuth_degrees, float), np.asarray(frequency_hz, float)
    az = np.unique(az)
    if len(az) > 1 and np.isclose(np.ptp(az), 360., atol=1e-8, rtol=0):
        az = az[:-1]
    if len(az) < 2:
        raise ValueError('Backprojection requires at least two distinct look directions; 0° and 360° are the same look.')
    if np.ptp(az) > 360. + 1e-8:
        raise ValueError('Select at most one azimuth revolution for coherent imaging.')
    steps = np.diff(az)
    full_circle = len(az) >= 3 and 360 - np.ptp(az) <= 1.5 * np.median(steps) + 1e-8
    center = float(center if center is not None else (0. if full_circle else np.mean(az)))
    theta = np.deg2rad(az - center)
    projection = abs(float(np.cos(np.deg2rad(elevation_degrees))))
    factor = 4 * np.pi * projection / C0
    # Conservative square contained in the frequency/angle sampling bounds.
    angle_step = np.max(np.hypot(np.diff(np.sin(theta)), np.diff(np.cos(theta))))
    radius = min(np.pi / (factor * np.max(np.diff(freq))),
                 np.pi / (factor * freq[-1] * angle_step))
    half = np.asarray(half if half is not None else [radius / np.sqrt(2)] * 2, float)
    u = np.outer(np.sin(theta), freq[[0, -1]])
    v = np.outer(np.cos(theta), freq[[0, -1]])
    resolution = 2 * np.pi / (factor * np.array([np.ptp(u), np.ptp(v)]))
    # Resolve the physical complex carrier as well as the magnitude envelope.
    spacing = np.minimum(resolution / 3, np.pi / (factor * np.array([abs(u).max(), abs(v).max()])))
    ratios = half / spacing
    if not np.all(np.isfinite(ratios)) or np.any(ratios <= 0):
        raise ValueError('Backprojection requires finite positive scene extents and pixel spacing.')
    counts = [max(1, int(np.ceil(value))) for value in ratios]
    from GRIM_Backend.plotting.modes.isar_mode import _validate_isar_working_set
    _validate_isar_working_set(120 * (2 * int(counts[0]) + 1) * (2 * int(counts[1]) + 1)
        + 32 * next_fast_len(32 * len(freq)) + 8*1024**2, operation='Coherent backprojection')
    axes = [np.linspace(-h, h, 2 * int(n) + 1) for h, n in zip(half, counts)]
    return theta, center, full_circle, factor, axes, resolution, radius


def working_set_bytes(axes, frequency_count):
    """Shared GUI/worker bound, including live complex output and row scratch."""
    return 120 * len(axes[0]) * len(axes[1]) + 32 * next_fast_len(32 * frequency_count) + 8*1024**2


def _cubic(values, coordinates):
    """Four-point periodic interpolation of an oversampled complex profile."""
    base = np.floor(coordinates).astype(np.int64)
    t = coordinates - base
    n = len(values)
    return (values[(base - 1) % n] * (-t * (t-1) * (t-2) / 6)
            + values[base % n] * ((t+1) * (t-1) * (t-2) / 2)
            - values[(base+1) % n] * ((t+1) * t * (t-2) / 2)
            + values[(base+2) % n] * ((t+1) * t * (t-1) / 6))


def project_row(samples, frequency_hz, ranges, factor, cancel_check=None):
    """Evaluate sum S(f) exp(+j factor f r), preserving absolute phase."""
    freq = np.asarray(frequency_hz, float)
    ranges = np.asarray(ranges, float)
    differences = np.diff(freq)
    output = np.empty(ranges.shape, complex)
    uniform = np.allclose(differences, differences[0], rtol=1e-10, atol=1e-5)
    if uniform:
        count = next_fast_len(32 * len(freq))
        profile = ifft(samples, n=count) * count
        for first in range(0, ranges.size, 16384):
            if cancel_check is not None and cancel_check():
                raise InterruptedError('ISAR computation superseded.')
            part = ranges[first:first+16384]
            coordinates = np.remainder(factor * differences[0] * part / (2*np.pi), 1.) * count
            output[first:first+16384] = _cubic(profile, coordinates) * np.exp(1j * factor * freq[0] * part)
    else:
        # Exact native-frequency evaluation, bounded to ~4 MiB phase blocks.
        output.fill(0)
        for first in range(0, ranges.size, 512):
            part = ranges[first:first+512]
            for start in range(0, len(freq), 256):
                if cancel_check is not None and cancel_check():
                    raise InterruptedError('ISAR computation superseded.')
                phase = np.exp(1j * factor * part[:, None] * freq[None, start:start+256])
                output[first:first+512] += phase @ samples[start:start+256]
    return output


def _cut_metrics(response, axis):
    power = abs(response)**2
    center = len(axis)//2
    if power[center] <= 0:
        return {'status': 'not_computed', 'reason': 'zero weighted support'}
    right = power[center:] / power[center]
    below = np.flatnonzero(right <= .5)
    fwhm = None
    if len(below) and below[0] > 0:
        i = below[0]
        fraction = (right[i-1]-.5)/(right[i-1]-right[i])
        fwhm = float(2*(axis[center+i-1] + fraction*(axis[center+i]-axis[center+i-1])))
    return {'status': 'partial', 'power_fwhm': fwhm,
            'reason': 'Origin cut on output grid; sidelobe metrics not evaluated'}


def form_band(params, indices, plan, resident_bytes=0):
    # Import lazily to share established source, units, storage and budget rules.
    from GRIM_Backend.plotting.modes import isar_mode as core
    dataset, freq = params['dataset'], np.asarray(params['freq_hz'])
    cancel, progress = params.get('cancel_check'), params.get('progress')
    az = core._angle_values_to_degrees(dataset, 'azimuth', np.asarray(dataset.azimuths)[indices])
    if params.get('az_center_deg') is not None:
        az = core._unwrap_degrees(az, params['az_center_deg'])
    order = np.argsort(az, kind='stable')
    az, ids = az[order], np.asarray(indices)[order]
    # A stored closing revolution endpoint is the same look, counted once.
    keep = np.r_[True, np.diff(az) > 1e-8]
    if len(az) > 2 and np.isclose(az[-1] - az[0], 360., atol=1e-8, rtol=0):
        keep[-1] = False
    az, ids = az[keep], ids[keep]
    theta, center, circular, factor, axes, resolution, radius = geometry(
        az, freq, params.get('elevation_deg', 0.), params.get('scene_half_extent_m'), params.get('az_center_deg'))
    x, y = axes
    cells = len(x) * len(y)
    native_bytes = 16 * len(ids) * len(freq) if params.get('retain_complex') else 0
    core._validate_isar_working_set(working_set_bytes(axes, len(freq)) + native_bytes,
        operation='Coherent backprojection', resident_bytes=resident_bytes)
    digest = core._selected_data_token(dataset, ids.tolist(), params['elev_idx'],
        params['freq_indices_sorted'], params['pol_idx'], cancel_check=cancel)
    if digest is None:
        raise InterruptedError('ISAR computation superseded.')
    # A full circle has no aperture endpoints to taper. Preserve all looks.
    wa = np.ones(len(az)) if circular else core._window_array(params['window_name'], len(az))
    wf = core._window_array(params['window_name'], len(freq))
    field = np.zeros(cells, complex)
    px, py = np.repeat(x, len(y)), np.tile(y, len(x))
    psf = np.zeros(len(x)+len(y), complex)
    gain, valid_count = 0., 0
    conjugate = core.plot_phase_conventions(dataset)[0] == '-jwt'
    native_samples = (core.allocate_result((len(ids), len(freq)), np.complex128)
                      if params.get('retain_complex') else None)
    for row, index in enumerate(ids):
        if cancel is not None and cancel():
            raise InterruptedError('ISAR computation superseded.')
        if progress is not None and row % 16 == 0:
            progress(f'Coherent backprojection: look {row+1} of {len(ids)}')
        samples = np.asarray(dataset.rcs_slice(np.ix_([index], [params['elev_idx']],
            params['freq_indices_sorted'], [params['pol_idx']])))[0, 0, :, 0]
        valid = np.isfinite(samples)
        if native_samples is not None:
            native_samples[row] = samples
        valid_count += int(valid.sum())
        weights = valid * wa[row] * wf
        gain += weights.sum()
        samples = np.where(valid, samples, 0.)
        if conjugate:
            samples = samples.conj()
        ranges = px * np.sin(theta[row]) + py * np.cos(theta[row])
        field += project_row(samples * weights, freq, ranges, factor, cancel)
        cut_ranges = np.r_[x*np.sin(theta[row]), y*np.cos(theta[row])]
        psf += project_row(weights, freq, cut_ranges, factor, cancel)
    if gain <= 0:
        raise ValueError('Backprojection requires finite complex samples with nonzero window weight.')
    image = core.allocate_result((len(x), len(y)), np.complex64)
    image[:] = (field / gain).reshape(image.shape)
    scale = params['unit_scale']
    contract = core.image_contract(theta, freq, center, params.get('elevation_deg', 0.), scale)
    contract.update(spatial_frequency_origin_hz=[0., 0.], spatial_frequency_step_hz=None,
                    reconstruction='backprojection', scene_frame='body_at_zero_azimuth' if center == 0 else 'look_frame')
    plan.update(scene_half_extent_m=[float(abs(x).max()), float(abs(y).max())],
                scene_frame=contract['scene_frame'], center_degrees=center)
    plan['warnings'].append('Backprojection coherently sums native looks; aspect-dependent scattering can interfere across a wide aperture.')
    if params.get('az_target_deg') is not None:
        plan['warnings'].append('Backprojection uses the selected native azimuth samples; the azimuth interpolation grid is not used.')
    magnitude = core.allocate_result(image.shape, np.float32)
    np.abs(image, out=magnitude)
    inverse = None
    if native_samples is not None:
        from .masking import reconstruction_extra
        inverse = dict(
            kind='backprojection', native_samples=native_samples,
            native_azimuth_degrees=az.copy(), native_frequency_hz=freq.copy(),
            look_center_degrees=center, elevation_degrees=params.get('elevation_deg', 0.),
            distance_scale_per_metre=scale, conjugate_time=conjugate,
            azimuths=np.asarray(dataset.azimuths)[ids].copy(),
            frequencies=np.asarray(dataset.frequencies)[params['freq_indices_sorted']].copy(),
            elevations=np.asarray(dataset.elevations)[[params['elev_idx']]].copy(),
            polarizations=[str(dataset.polarizations[params['pol_idx']])],
            units=dict(dataset.units or {}), extra=reconstruction_extra(dataset),
            history=str(dataset.history or ''), window=params['window_name'],
        )
    return dict(az_values=az, x_range=x*scale, y_range=y*scale, magnitude=magnitude,
        **({'complex_image': image} if params.get('retain_complex') else {}),
        **({'mask_inverse': inverse} if inverse is not None else {}),
        image_contract=contract, phase_coverage=valid_count/(len(ids)*len(freq)),
        source_selection_digest=digest.hex(), accurate_pfa=False,
        sampling=dict(cross_resolution=float(resolution[0]*scale), range_resolution=float(resolution[1]*scale),
            cross_half_extent=float(radius/np.sqrt(2)*scale), range_half_extent=float(radius/np.sqrt(2)*scale),
            definition='Inverse projected spatial-frequency spans; sampling extents are conservative square bounds, not a guarantee against ambiguity.'),
        psf={'definition': 'Origin cuts from the actual tapered native polar support, sampled on the output grid',
             'cross_range': _cut_metrics(psf[:len(x)], x*scale),
             'range': _cut_metrics(psf[len(x):], y*scale)},
        spatial_frequency_support={'coverage_interpolation': 'native angles and frequencies; 32x range-profile cubic interpolation for uniform frequencies',
            'angular_window': 'uniform over full circle' if circular else params['window_name'],
            'sample_weighting': 'equal acquired samples with aperture taper; missing samples zero weighted'})
