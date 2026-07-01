@echo off
cd /d "%~dp0"
echo Installing dependencies...
python -m pip install -r requirements.txt --quiet
echo.
echo Starting HelixZero-CMS API on http://localhost:8000
echo The browser will open automatically.
echo Press Ctrl+C to stop.
echo.
start "HelixZero-CMS API" cmd /k "uvicorn api.main:app --reload --port 8000"
timeout /t 2 /nobreak >nul
start "" http://localhost:8000
pause
