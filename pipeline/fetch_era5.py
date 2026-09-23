"""
Fetch ERA5 hourly total precipitation over Spain (1976 -> latest) and store daily sums.

Source: Copernicus Climate Change Service (C3S) ERA5 reanalysis, as mirrored
(hourly, 0.25 deg, unchanged values) in the public Open-Meteo S3 bucket
s3://openmeteo/data/copernicus_era5/precipitation/ . No API key needed.

If you have a CDS account you can instead use fetch_cds.py (same output format).

Output: pipeline/raw/<tag>.npz with
    daily  float32 (lat, lon, days)   mm/day
    start  ISO date of first day
"""
import os, sys, json, datetime as dt
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

BUCKET = "openmeteo/data/copernicus_era5/precipitation"
# Global grid: lat -90..90 (721), lon -180..179.75 (1440), 0.25 deg
BOXES = {  # south, north, west, east
    "mainland": (35.5, 44.0, -9.5, 4.5),    # peninsula + Balearic Islands
    "canarias": (27.5, 29.5, -18.25, -13.25),
}
BOX = os.environ.get("BOX", "mainland")
LAT0, LAT1, LON0, LON1 = BOXES[BOX]
R0, R1 = int((LAT0 + 90) / 0.25), int((LAT1 + 90) / 0.25) + 1
C0, C1 = int((LON0 + 180) / 0.25), int((LON1 + 180) / 0.25) + 1
ROW_BLOCK = 7
OUT = os.path.join(os.path.dirname(__file__), "raw" if BOX == "mainland" else f"raw_{BOX}")
FIRST_YEAR = 1976
LAST_YEAR_FILE = 2021       # later data lives in 504-hour chunk files
CHUNK_H = 504


def reader(key):
    import s3fs, omfiles
    fs = s3fs.S3FileSystem(anon=True, default_block_size=256 * 1024, default_cache_type="none")
    return omfiles.OmFileReader.from_fsspec(fs, f"{BUCKET}/{key}")


def task_year(year, r0, r1):
    tag = f"y{year}_r{r0}"
    path = os.path.join(OUT, tag + ".npz")
    if os.path.exists(path):
        return tag
    r = reader(f"year_{year}.om")
    nt = r.shape[2]
    a = r.read_array((slice(r0, r1), slice(C0, C1), slice(0, nt)))
    a = np.nan_to_num(a)
    # hour index 0 = YYYY-01-01T00:00 (accumulation of preceding hour); sum per UTC day
    days = nt // 24
    daily = a[:, :, : days * 24].reshape(a.shape[0], a.shape[1], days, 24).sum(-1)
    np.savez_compressed(path, daily=daily.astype(np.float32), start=f"{year}-01-01", r0=r0)
    return tag


def task_chunk(ci):
    tag = f"c{ci}"
    path = os.path.join(OUT, tag + ".npz")
    if os.path.exists(path):
        return tag
    r = reader(f"chunk_{ci}.om")
    nt = r.shape[2]
    a = np.nan_to_num(r.read_array((slice(R0, R1), slice(C0, C1), slice(0, nt))))
    np.savez_compressed(path, hourly=a.astype(np.float32), start_hour=ci * CHUNK_H)
    return tag


def main():
    os.makedirs(OUT, exist_ok=True)
    import s3fs
    fs = s3fs.S3FileSystem(anon=True)
    meta = json.loads(fs.cat(f"openmeteo/data/copernicus_era5/static/meta.json"))
    end_h = meta["data_end_time"] // 3600
    first_chunk = int(dt.datetime(LAST_YEAR_FILE + 1, 1, 1, tzinfo=dt.timezone.utc).timestamp() // 3600) // CHUNK_H
    last_chunk = end_h // CHUNK_H
    jobs = []
    with ProcessPoolExecutor(max_workers=int(sys.argv[1]) if len(sys.argv) > 1 else 24) as ex:
        for ci in range(first_chunk, last_chunk + 1):
            jobs.append(ex.submit(task_chunk, ci))
        for y in range(FIRST_YEAR, LAST_YEAR_FILE + 1):
            for r0 in range(R0, R1, ROW_BLOCK):
                jobs.append(ex.submit(task_year, y, r0, min(r0 + ROW_BLOCK, R1)))
        done = 0
        for f in as_completed(jobs):
            done += 1
            try:
                print(f"[{done}/{len(jobs)}] {f.result()}", flush=True)
            except Exception as e:
                print(f"[{done}/{len(jobs)}] ERROR {e!r}", flush=True)
    json.dump({"end_hour": end_h, "first_chunk": first_chunk, "last_chunk": last_chunk},
              open(os.path.join(OUT, "meta.json"), "w"))


if __name__ == "__main__":
    main()
