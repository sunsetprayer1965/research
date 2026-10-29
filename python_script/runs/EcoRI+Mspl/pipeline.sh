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