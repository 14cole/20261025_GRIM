# GRIM

Desktop electromagnetic data analysis with integrated GHOST solvers, FREDDY
material-stack tools, assembly workflows, plotting, ISAR, and PowerPoint reporting.

## Install and run

Use 64-bit Python 3.11 or later. From this repository's root on Windows:

```powershell
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe GRIM_Backend\run_gui.py
```

After installation, `Launch_GRIM.bat` opens the application using that environment.
Keep `GRIM_Backend` and `tools` together: the GHOST and FREDDY tabs load their
implementations from the bundled subfolders.

On other desktop platforms, create a virtual environment with `python3 -m venv
.venv`, install the same requirements with `.venv/bin/python -m pip install -r
requirements.txt`, and run `.venv/bin/python GRIM_Backend/run_gui.py`.

PowerPoint automation requires Windows and an installed desktop copy of
Microsoft PowerPoint. Other analysis and plotting features do not require it.

## Included runtime

- `GRIM_Backend`: desktop application, dataset operations, imaging, and reports.
- `tools/GHOST`: 2D and body-of-revolution solvers and their interfaces.
- `tools/FREDDY`: material stacks, impedance calculations, and material exploration.
- Application assets, report templates, built-in material presets, and native
  accelerator source/build helpers.

The solvers have NumPy fallbacks; machine-specific compiled accelerator binaries
are not included. Optional accelerators can be built with the helpers under
`tools/GHOST/ghost_backend/twod/assembly/native` and
`tools/GHOST/ghost_backend/bor/native`.

Tests, audit archives, standalone validation cases, recovery files, and generated
datasets are omitted. The small validation service modules required by the
existing GHOST integration remain included. Third-party license notices are
retained with their runtime code.

The runtime includes the box-selection correction: increasing azimuth bounds
retain literal endpoints, so 90–360 does not include 0 and 0–30 does not include
360. A stop below the start still selects a range across the angular seam.
