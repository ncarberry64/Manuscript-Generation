@echo off
setlocal
pushd "%~dp0"
if not exist ".venv\Scripts\python.exe" (
    python -m venv .venv
    if errorlevel 1 goto error
    ".venv\Scripts\python.exe" -m pip install -e ".[dev]"
    if errorlevel 1 goto error
)
".venv\Scripts\python.exe" -m bhtracker serve --open
popd
exit /b 0
:error
echo Setup failed. See README.md for the explicit installation commands.
pause
popd
exit /b 1
