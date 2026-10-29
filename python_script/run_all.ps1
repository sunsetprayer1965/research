$ErrorActionPreference = "Stop"

$ROOT = "C:\Users\weilu\Desktop\research william luo\python_script"
$DATASET_ROOT = Join-Path $ROOT "dataset_linking_number"
$RUNS_ROOT = Join-Path $ROOT "runs"

# ===== 只跑一个 dataset（先测试这个）=====
#$datasets = @("HindIII+EcoRI")

# ===== 如果以后想跑全部，把上面注释掉，改成下面这行 =====
$datasets = Get-ChildItem -Path $DATASET_ROOT -Directory | Select-Object -ExpandProperty Name

New-Item -ItemType Directory -Force -Path $RUNS_ROOT | Out-Null

# ===== 先编译 Java（只一次）=====
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

    # ===== 清理旧结果，避免干扰 =====
    if (Test-Path $runDir) {
        Remove-Item $runDir -Recurse -Force
    }

    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    New-Item -ItemType Directory -Force -Path (Join-Path $runDir "logs") | Out-Null
    New-Item -ItemType Directory -Force -Path (Join-Path $runDir "linking") | Out-Null

    # ===== 把 config.ini 和 bed 复制到 runDir =====
    Copy-Item (Join-Path $ROOT "config.ini") $runDir -Force
    Copy-Item (Join-Path $ROOT "duan.SC.10000.raw_sub.bed") $runDir -Force

    # ===== Docker 内按顺序跑完整链路 =====
    # ===== 生成容器内要跑的脚本文件（写到 runDir，也就是 /app）=====
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

# ✅ 关键：强制 LF（避免 \r）
$cmd = $cmd -replace "`r`n", "`n"

# ✅ 写脚本到 runDir（/app），用 UTF8 不带 BOM 更稳
$utf8NoBom = New-Object System.Text.UTF8Encoding($false)
[System.IO.File]::WriteAllText($shPath, $cmd, $utf8NoBom)

# ===== Docker 执行脚本文件（不再把多行塞进 -lc 参数）=====
docker run --rm `
  -v "${runDir}:/app" `
  -v "${dsDir}:/data:ro" `
  -v "${ROOT}:/repo:ro" `
  pastis bash /app/pipeline.sh

# ✅ 必须检查 docker 成败，否则你会继续跑 Java 然后假装 Finished
if ($LASTEXITCODE -ne 0) {
    throw "Docker pipeline failed for dataset [$ds], exit code = $LASTEXITCODE"
}
    # ===== Java 本机：读取 split 结果，输出 linking =====
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