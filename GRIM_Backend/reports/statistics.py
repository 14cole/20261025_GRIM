"""Statistical report panels on exact common azimuth samples, without Qt.

Every panel describes one frequency, elevation, and polarization. Overlay
populations use the same directions, count a periodic seam once, and preserve
physical zero power separately from invalid samples.
"""
from __future__ import annotations

import math
import numpy as np

from GRIM_Backend.plotting import distribution_analysis as distribution
from GRIM_Backend.plotting.sweep_analysis import reduce_power
from GRIM_Backend.reports import plot_data as cuts
from GRIM_Backend.reports.report import PlotSeries, PlotSpec


STATISTICAL_KINDS = frozenset({"cdf", "exceedance", "sector_stats"})
DISTRIBUTION_KINDS = frozenset({"cdf", "exceedance"})


def build_statistical_specs(datasets, *, frequencies, elevation, polarization,
                            kind="cdf", sector_width_deg=30., statistic="mean",
                            percentile=90., show_legend=True, tol=1e-6):
    """Build per-frequency distributions or equal-width azimuth sectors.

    CDFs rank displayed levels with equal sample weights. Sector statistics
    reduce linear power before converting to the dataset's displayed dB unit;
    percentile uses nearest rank. Empty sectors remain gaps.
    """
    if kind not in STATISTICAL_KINDS:
        raise ValueError(f"Unsupported statistical plot kind: {kind!r}")
    selected = cuts._coerce_named_grids(datasets)
    available = cuts.get_plot_availability(selected, tol=tol, evaluate_phase=False)
    frequencies = tuple(float(f) for f in frequencies)
    if not frequencies:
        raise ValueError("Select at least one frequency for statistical plots.")
    cuts._exact_numeric_indices(available.frequencies, frequencies, tol=tol,
                               axis_name="common frequency", dataset_name="selected overlays")
    cuts._exact_index(available.elevations, elevation, tol=tol,
                      axis_name="common elevation", dataset_name="selected overlays")
    polarizations = cuts._requested_polarizations(polarization, available.polarizations)
    if len(available.azimuths) < 2:
        raise ValueError("Statistical reports require at least two common azimuth samples.")
    band = (min(available.azimuths), max(available.azimuths))
    azimuths = cuts._common_azimuths_in_band(available.azimuths, band,
                                          angle_unit=available.azimuth_unit, tol=tol)
    if len(azimuths) < 2:
        raise ValueError("Statistical reports require at least two distinct common azimuth directions.")
    degrees = cuts._convert_angles(azimuths, available.azimuth_unit, "deg")
    low, high = cuts._convert_angles(band, available.azimuth_unit, "deg")
    if kind == "sector_stats":
        width = float(sector_width_deg)
        if not math.isfinite(width) or not 1 <= width <= 360:
            raise ValueError("Sector width must be between 1 and 360 degrees.")
        if statistic not in {"mean", "median", "max", "percentile"}:
            raise ValueError("Choose mean, median, maximum, or percentile sector statistics.")
        if not math.isfinite(float(percentile)) or not 0 <= percentile <= 100:
            raise ValueError("Sector percentile must be between 0 and 100.")
        edges = np.r_[np.arange(low, high, width), high]
        config = dict(mode=statistic, percentile=float(percentile), interpolate=False)
        statistic_label = {"mean": "Mean power", "median": "Median power", "max": "Maximum",
                           "percentile": f"P{percentile:g} nearest rank"}[statistic]
    else:
        config = distribution.options(dict(exceedance=kind == "exceedance"))
    prepared = []
    for item in selected:
        grid = item.grid
        ai = cuts._exact_numeric_indices(grid.azimuths, azimuths, tol=tol,
                                        axis_name="azimuth", dataset_name=item.name)
        ei = cuts._exact_index(grid.elevations, elevation, tol=tol,
                              axis_name="elevation", dataset_name=item.name)
        fi = cuts._exact_numeric_indices(grid.frequencies, frequencies, tol=tol,
                                        axis_name="frequency", dataset_name=item.name)
        prepared.append((item, ai, ei, fi))
    result = []
    elevation_deg = cuts._convert_angles([elevation], available.elevation_unit, "deg")[0]
    for pol in polarizations:
        for panel, frequency in enumerate(frequencies):
            series, curves = [], []
            has_zero = has_invalid = False
            for item, ai, ei, fi in prepared:
                grid = item.grid
                pi = cuts._exact_index(grid.polarizations, pol, tol=0.,
                                      axis_name="polarization", dataset_name=item.name)
                selection = (np.asarray(ai), np.array([ei]), np.array([fi[panel]]), np.array([pi]))
                if kind in DISTRIBUTION_KINDS:
                    curve = distribution.prepare_curve(item.name, grid, selection,
                        "All selected coordinates", config, scale="dbsm")
                    if curve is None:
                        raise ValueError(f"Dataset {item.name!r} has no valid statistical samples on this cut.")
                    curves.append(curve)
                    has_zero |= curve["population"]["counts"]["zero"] > 0
                    has_invalid |= curve["population"]["counts"]["excluded"] > 0
                else:
                    power = np.asarray(grid.power[np.ix_(*selection)], dtype=float).reshape(-1)
                    valid = np.isfinite(power) & (power >= 0)
                    if not np.any(valid):
                        raise ValueError(f"Dataset {item.name!r} has no valid statistical samples on this cut.")
                    has_invalid |= not np.all(valid)
                    x, y = [], []
                    for index, (start, stop) in enumerate(zip(edges[:-1], edges[1:])):
                        members = (degrees >= start) & ((degrees <= stop) if index == len(edges)-2 else (degrees < stop))
                        values = power[members]
                        if values.size:
                            reduced, _count = reduce_power(values, config)
                            level = float(grid.linear_to_default_db(np.asarray(reduced),
                                frequency_value=float(grid.frequencies[fi[panel]]), eps=0.))
                        else:
                            level = math.nan
                        has_zero |= bool(np.isneginf(level))
                        x.extend((float(start), float(stop), float(stop)))
                        y.extend((level, level, math.nan))
                    series.append(PlotSeries(tuple(x), tuple(y),
                        label=cuts._series_label(item, "magnitude", available)))
            if kind in DISTRIBUTION_KINDS:
                limits = distribution.bounds(curves)
                # A constant distribution still needs a readable level axis;
                # nanodecibel padding also collapses in the UI's spin boxes.
                if limits[1] - limits[0] < 1e-5:
                    center = (limits[0] + limits[1]) / 2
                    limits = (center - 1., center + 1.)
                for item, curve in zip(selected, curves):
                    x, y = distribution.draw_values(curve, limits, exceedance=kind == "exceedance")
                    series.append(PlotSeries(tuple(x), tuple(y),
                        label=cuts._series_label(item, "magnitude", available),
                        draw_style="steps-pre" if kind == "exceedance" else "steps-post"))
            frequency_ghz = cuts._convert_frequencies([frequency], available.frequency_unit, "GHz")[0]
            title = f"{frequency_ghz:g} GHz | Elevation {elevation_deg:g} deg | {pol}"
            note = []
            if has_zero:
                note.append("Zero power retained at −∞ dB." if kind in DISTRIBUTION_KINDS else "Zero-valued sectors: −∞ dB (not drawn).")
            if has_invalid:
                note.append("Invalid samples excluded.")
            distribution_kind = kind in DISTRIBUTION_KINDS
            result.append(PlotSpec(
                plot_id=f"{kind}_{cuts._safe_id(pol)}_{panel:03d}", kind=kind,
                title=title if distribution_kind else f"{title}\n{statistic_label} | {sector_width_deg:g}° sectors",
                x_label=cuts._quantity_axis_label(selected, "magnitude", available) if distribution_kind else "Azimuth (deg)",
                y_label=distribution.y_label(config) if distribution_kind else cuts._quantity_axis_label(selected, "magnitude", available),
                series=tuple(series), x_limits=limits if distribution_kind else (float(low), float(high)),
                y_limits=(0., 100.) if distribution_kind else None,
                show_legend=bool(show_legend), note="\n".join(note)))
    return tuple(result)
