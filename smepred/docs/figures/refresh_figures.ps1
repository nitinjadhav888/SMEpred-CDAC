# Refresh all figures by running the standalone figure generation scripts.
# Run from the HelixZero-CMS project root (Helixx/helixzero_cms).

$ProjectRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$FigsDir = $PSScriptRoot
$Scripts = @(
    "schematic_diagram.py",
    "venn_diagrams.py",
    "workflow_diagram.py"
)

# Try to activate venv if it exists
$VenvActivate = Join-Path $ProjectRoot ".venv\Scripts\Activate.ps1"
if (Test-Path $VenvActivate) {
    . $VenvActivate
} else {
    # Check parent directory venv
    $ParentVenv = Join-Path (Split-Path -Parent $ProjectRoot) ".venv\Scripts\Activate.ps1"
    if (Test-Path $ParentVenv) {
        . $ParentVenv
    }
}

Write-Host "`nRefreshing figures in $FigsDir ...`n" -ForegroundColor Cyan

Push-Location $FigsDir
try {
    foreach ($Script in $Scripts) {
        $ScriptPath = Join-Path $FigsDir $Script
        if (Test-Path $ScriptPath) {
            Write-Host "Running $Script ..." -NoNewline
            python $ScriptPath 2>&1 | ForEach-Object { "$_" }
            if ($LASTEXITCODE -eq 0) {
                Write-Host "  OK" -ForegroundColor Green
            } else {
                Write-Host "  FAILED (exit code: $LASTEXITCODE)" -ForegroundColor Red
            }
        } else {
            Write-Host "  $Script not found, skipping" -ForegroundColor Yellow
        }
    }
} finally {
    Pop-Location
}

Write-Host "`nDone.`n" -ForegroundColor Cyan
