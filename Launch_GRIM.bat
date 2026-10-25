@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
    echo Install the environment using the commands in README.md first.
    pause
    exit /b 1
)
".venv\Scripts\python.exe" "GRIM_Backend\run_gui.py"
if errorlevel 1 pause
