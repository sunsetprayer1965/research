$ErrorActionPreference = "Stop"

$ROOT = "C:\Users\weilu\Desktop\research william luo\python_script"
$DATASET_ROOT = Join-Path $ROOT "dataset_linking_number"
$RUNS_ROOT = Join-Path $ROOT "runs"

# ===== Process only one dataset (for an initial test) =====
#$datasets = @("HindIII+EcoRI")

# ===== To process all datasets, uncomment the line above and use this line instead. =====
$datasets = Get-ChildItem -Path $DATASET_ROOT -Directory | Select-Object -ExpandProperty Name

New-Item -ItemType Directory -Force -Path $RUNS_ROOT | Out-Null

# ===== Compile Java once before processing datasets. =====
Push-Location (Join-Path $ROOT "java\computewritec")
javac -encoding UTF-8 -d . ComputeWriteC3.java
Pop-Location

foreach ($ds in $datasets) {
    $runDir = Join-Path $RUNS_ROOT $ds
    $dsDir  = Join-Path $DATASET_ROOT $ds

    Write-Host ""
    Write-Host "==============================="
    Write-Host "Dataset: $ds"
    Write-Host "Input : $dsDir"
    Write-Host "Output: $runDir"
    Write-Host "==============================="

    if (-not (Test-Path $dsDir)) {
        Write-Host "❌ Dataset folder not found: $dsDir"
        continue
    }

    # ===== Remove old results to avoid mixing them with this run. =====
    if (Test-Path $runDir) {
        Remove-Item $runDir -Recurse -Force
    }

    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    New-Item -ItemType Directory -Force -Path (Join-Path $runDir "logs") | Out-Null
    New-Item -ItemType Directory -Force -Path (Join-Path $runDir "linking") | Out-Null

    # ===== Copy config.ini and the BED file into runDir. =====
    Copy-Item (Join-Path $ROOT "config.ini") $runDir -Force
    Copy-Item (Join-Path $ROOT "duan.SC.10000.raw_sub.bed") $runDir -Force

    # ===== Run the full pipeline in Docker, in order. =====
    # ===== Write the container script into runDir, which is mounted as /app. =====
$shPath = Join-Path $runDir "pipeline.sh"

$cmd = @'
set -euo pipefail
cd /app
mkdir -p /app/logs

echo "[INFO] locus..."
python3 -u /repo/locus.py --dataset_dir /data --bed /app/duan.SC.10000.raw_sub.bed --out /app/duan.SC.10000.raw_sub.matrix > /app/logs/locus.log 2>&1

echo "[INFO] pastis-pm2..."
pastis-pm2 . > /app/logs/pastis_pm2.log 2>&1

echo "[INFO] pastis-pm1..."
pastis-pm1 . > /app/logs/pastis_pm1.log 2>&1

echo "[INFO] pastis-mds..."
pastis-mds . > /app/logs/pastis_mds.log 2>&1

echo "[INFO] pastis-nmds..."
pastis-nmds . > /app/logs/pastis_nmds.log 2>&1

echo "[INFO] coordinates..."
python3 -u /repo/coordinates_blank.py > /app/logs/coordinates.log 2>&1

echo "[INFO] split..."
python3 -u /repo/java/chromosome/split.py --base_path . > /app/logs/split.log 2>&1

echo "[INFO] container pipeline done."
'@

# Force LF line endings to avoid carriage returns.
$cmd = $cmd -replace "`r`n", "`n"

# Write the script into runDir (/app) as UTF-8 without a BOM.
$utf8NoBom = New-Object System.Text.UTF8Encoding($false)
[System.IO.File]::WriteAllText($shPath, $cmd, $utf8NoBom)

# ===== Run the script file in Docker instead of passing multiline text to -lc. =====
docker run --rm `
  -v "${runDir}:/app" `
  -v "${dsDir}:/data:ro" `
  -v "${ROOT}:/repo:ro" `
  pastis bash /app/pipeline.sh

# Check Docker's exit code so Java does not run after a failed pipeline.
if ($LASTEXITCODE -ne 0) {
    throw "Docker pipeline failed for dataset [$ds], exit code = $LASTEXITCODE"
}
    # ===== Run Java locally to read the split results and write linking frequencies. =====
    $resultDir = Join-Path $runDir "linking"
    $mds  = Join-Path $runDir "MDS_split"
    $nmds = Join-Path $runDir "NMDS_split"
    $pm1  = Join-Path $runDir "PM1_split"
    $pm2  = Join-Path $runDir "PM2_split"

    java -cp (Join-Path $ROOT "java\computewritec") computewritec.ComputeWriteC3 `
        $resultDir $mds $nmds $pm1 $pm2 *> (Join-Path $runDir "logs\java_linking.log")

    Write-Host "✅ Finished: $ds"
}

Write-Host ""
Write-Host "🎉 All datasets done."
