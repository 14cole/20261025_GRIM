"""Explicit adapters from producer coordinates to the imaging phase law.

GHOST 2-D B uses +j 2k (X cos(a) + Y sin(a)). The canonical imager
uses -j 2k (x sin(b) + y cos(b)); b=270-a preserves native body X/Y.
Only the identified GHOST field domain (or explicit versioned declaration)
enables this adapter. Unspecified third-party data is never guessed.
"""
from __future__ import annotations

import numpy as np

GHOST_FRAME = "ghost_2d_positive_cos_sin_v1"


def is_ghost_frame(dataset):
    extra = dataset.extra or {}
    explicit = extra.get("imaging_coordinate_convention")
    if explicit is not None:
        return str(explicit) == GHOST_FRAME
    from GRIM_Backend.datasets.phase_conventions import declared_time_sign
    try:
        return (str(extra.get("complex_field_domain", "")) ==
                "2d_layer_potential_bare_integral_amplitude_B"
                and declared_time_sign(dataset) == "+jwt")
    except ValueError:
        return False


class _ImagingView:
    def __init__(self, source):
        from GRIM_Backend.plotting.modes.isar_mode import _angle_values_to_degrees
        self.source = source
        self.azimuths = 270. - _angle_values_to_degrees(source, "azimuth", source.azimuths)
        self.units = dict(source.units, azimuth="deg")
        self.extra = dict(source.extra, imaging_coordinate_convention=GHOST_FRAME)

    def __getattr__(self, name):
        return getattr(self.source, name)


def processing_params(params):
    if not is_ghost_frame(params["dataset"]):
        return params
    result = dict(params, dataset=_ImagingView(params["dataset"]))
    if params.get("az_center_deg") is not None:
        result["az_center_deg"] = 270. - params["az_center_deg"]
    if params.get("az_target_deg") is not None:
        result["az_target_deg"] = np.sort(270. - np.asarray(params["az_target_deg"]))
    return result


def restore_source_contract(result, source):
    if isinstance(result, str) or not is_ghost_frame(source):
        return result
    from .masking import reconstruction_extra
    for band in result[0]:
        band["az_values"] = np.sort(270. - band["az_values"])
        band["image_contract"].update(source_coordinate_convention=GHOST_FRAME,
            source_to_processing_azimuth="270 - native azimuth (degrees)",
            body_coordinate_frame="native_geometry_xy")
        band.setdefault('isar_metadata_notes', []).append(
            'GHOST native body X/Y preserved using processing azimuth 270° − native azimuth; coherent axes use the declared look-frame basis.')
        inverse = band.get("mask_inverse")
        if inverse is not None:
            native = 270. - np.asarray(inverse["azimuths"])
            if str(source.units.get("azimuth", "deg")).lower() in {"rad", "radian", "radians"}:
                native = np.deg2rad(native)
            inverse["azimuths"] = native
            inverse["units"] = dict(source.units)
            inverse["extra"] = reconstruction_extra(source)
            inverse["extra"]["imaging_coordinate_convention"] = GHOST_FRAME
    return result
