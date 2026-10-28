MATLAB GRIM plotting tool -- R2024b
=================================

START
Open this folder in MATLAB and run:

    app = run_grim;

GrimPlotApp also starts the complete app directly. GrimPlotApp('gui') retains
the original empty GUI option. The complete folder must be on the MATLAB path.
Base MATLAB is sufficient. There is no Python call, interpreter, package,
bridge, or Python installation required by the app or the MATLAB tests.

WORKFLOW
Load one or multiple files. Select dataset rows, then choose polarization,
frequency, elevation and azimuth values. Ctrl/Shift selects multiple values.
Click Plot, or enable Auto Plot. The first selected table row supplies the
parameter coordinates; other datasets are matched by values, not indices.
All internal angular coordinates are degrees; frequencies are GHz.

Operations create new dataset rows and leave their inputs intact. Crop uses
the current parameter selection. Other transforms operate on the complete
selected datasets. Statistics are calculated on linear power and broadcast
back over the original axes. Percentile defaults to nearest rank across
azimuth. Align matches the first selected dataset's existing coordinates;
Regrid interpolates without extrapolation. Join rejects conflicting overlaps.
Operand order is table-row order, not selection-click order.

PLOTS
Azimuth rectangular/polar, frequency, elevation sweep, waterfall, RF Compare,
Delta Map, and sector statistics are implemented. Sweep plots draw individual
cuts. Waterfall and Delta Map use azimuth versus frequency; select exactly one
elevation and polarization. Delta Map needs two compatible datasets. RF Compare
needs two datasets, one varying axis and one polarization; it uses the Python
GRIM concordance/error/sector-weighted score. Sector statistics shows the mean
linear power over the selected azimuth population for each cut. Advanced Python
plot-specific dialogs and interactive analysis overlays are outside this
simplified GUI. Zoom/pan/data tips use MATLAB's built-in axes interactions.

INPUT / OUTPUT
Input: native .grim (NPY-in-ZIP, numeric/string payloads, versions 1 and 2),
MATLAB datasets saved by this app, GRIM flat CSV/TSV, Pioneer PIO/CMPLX_DI,
PTM, OUT, Xpatch SS, and the supported CST/SENTRi table schemas.
Unknown arbitrary tables are rejected rather than guessing column physics.
Legacy object/pickle-based .grim archives are not supported.

Save: .grim or .mat. Export: CSV, TSV, PIO, PTM; plot export: PNG, PDF, JPEG,
TIFF. PIO/PTM export uses the full azimuth/frequency cut at the selected single
elevation/polarization. Both require 3-D RCS and known phase. PTM additionally
requires finite float32-representable samples, uniform axes, no repeated full-
turn endpoint, and at least 37 frequencies. Native/CSV preserve unknown phase.
Load uses a file chooser; no third-party drag/drop helper is included.

NUMERICAL CONVENTIONS
Array order: [azimuth, elevation, frequency, polarization]. Power is physical
sigma_3d (m^2), sigma_2d scattering width (m), or a dimensionless power ratio.
Phase is radians. NaN power means absent data; NaN phase means unknown phase.
Coherent operations refuse positive-power samples with unknown phase.
Incoherent/statistical operations do not invent phase.
GRIM dBke = 10*log10(k0*scattering_width), k0=2*pi*f/c. It is not simply
10*log10(scattering_width). Raw native solver fields are normalized on import.
CST theta is converted to elevation=90-theta. SENTRi native theta is retained
and tagged; this tool does not include the excluded beta coordinate converter.

FILES / PYTHON SOURCE MAP
GrimPlotApp.m    MATLAB UI (original Python ui/app.py and dataset_sidebar.py)
GrimBackend.m    UI/controller dispatch (ui/dataset_actions.py and dialogs)
GrimData.m       Array model and math (datasets/grid, axes, arithmetic,
                combine, transforms; MATLAB in-memory implementation)
GrimIO.m         Native/flat table dispatch (io/native.py, csv.py, loaders.py)
GrimNpy.m        MATLAB implementation of the NPY/NPZ serialization used by GRIM
GrimVendorIO.m   io/pioneer.py, ptm.py, out.py, xpatch.py, cst.py, sentri.py
GrimPlots.m      plotting/modes and RF comparison calculations
run_grim.m      Launcher

The port uses ordinary in-memory MATLAB arrays. Python's disk-backed storage,
background worker/cancellation infrastructure, caches and assembly/HPC features
are not part of this plotting tool. Large files must fit available MATLAB RAM.

VALIDATION
From this folder:

    results = runtests('tests');
    assertSuccess(results);

Tests cover native/raw normalization, reference results for arithmetic and
transforms, file readers, interchange round trips, missing-phase rejection,
RF agreement, and all eight GUI plot modes. Reference .mat/.grim/vendor files
were generated with the repository's Python implementation during development;
they are fixed data and do not execute Python. Numeric tests can run in batch;
the GUI smoke test is skipped when the MATLAB desktop is unavailable.

Development status: all .m files pass the MISS_HIT MATLAB syntax checker.
MATLAB R2024b was not available in the development environment, so the MATLAB
tests and visual/runtime verification must still be run in MATLAB. This is a
source port, not a claim of complete parity with every Python edge case or
vendor variation. No ISAR, beta buttons, provenance UI or axis-units UI is added.
