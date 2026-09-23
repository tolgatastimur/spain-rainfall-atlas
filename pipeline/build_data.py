"""
Turn the raw ERA5 daily sums (pipeline/raw, pipeline/raw_canarias) into compact files for the web page.

Outputs in site/data/:
  meta.json   grid cells, dates, yearly totals per cell, national series, headline stats
  daily.b64.txt  uint8 [day][cell]  daily precipitation, q = round(sqrt(mm) * Q), base64
  clim.b64.txt   uint16 [365][cell] 1991-2020 daily normal (31-day smoothed), mm * 100, base64
"""
import os, json, glob, datetime as dt
import numpy as np
from shapely.geometry import shape, Point
from shapely.prepared import prep

HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "..", "site", "data")
RES = 0.25
FIRST_YEAR, LAST_YEAR_FILE, CHUNK_H = int(os.environ.get("FIRST_YEAR", 1976)), 2021, 504
Q = 14.0                     # daily encoding: max ~331 mm
NORMAL = tuple(int(v) for v in os.environ.get("NORMAL", "1991,2020").split(","))  # WMO climate normal
# raw download boxes: directory, south lat, west lon, n lat, n lon (see fetch_era5.py)
BOXES = {"mainland": ("raw", 35.5, -9.5, 35, 57), "canarias": ("raw_canarias", 27.5, -18.25, 9, 21)}
ISLANDS = {"Baleares": "baleares", "Canarias": "canarias"}
# map groups: peninsula as 0.5 deg cells (2x2 ERA5 points), small islands at native 0.25 deg
GROUPS = [("mainland", "mainland", 2), ("baleares", "mainland", 1), ("canarias", "canarias", 1)]


def load_daily(raw, nlat, nlon):
    raw = os.path.join(HERE, raw)
    parts = []
    for y in range(FIRST_YEAR, LAST_YEAR_FILE + 1):
        files = sorted(glob.glob(os.path.join(raw, f"y{y}_r*.npz")), key=lambda f: int(np.load(f)["r0"]))
        assert files, f"missing year {y} in {raw}"
        parts.append(np.concatenate([np.load(f)["daily"] for f in files], axis=0))
        assert parts[-1].shape[:2] == (nlat, nlon), (y, parts[-1].shape)
    meta = json.load(open(os.path.join(raw, "meta.json")))
    hourly = np.concatenate([np.load(os.path.join(raw, f"c{c}.npz"))["hourly"]
                             for c in range(meta["first_chunk"], meta["last_chunk"] + 1)], axis=2)
    h0 = meta["first_chunk"] * CHUNK_H
    start_2022 = int(dt.datetime(LAST_YEAR_FILE + 1, 1, 1, tzinfo=dt.timezone.utc).timestamp() // 3600)
    hourly = hourly[:, :, start_2022 - h0: meta["end_hour"] + 1 - h0]
    ndays = hourly.shape[2] // 24          # complete UTC days only
    parts.append(hourly[:, :, : ndays * 24].reshape(nlat, nlon, ndays, 24).sum(-1))
    return np.clip(np.concatenate(parts, axis=2), 0, None)


def make_cells(group, daily, lat0, lon0, block, preg):
    """Aggregate ERA5 points into block x block cells that overlap the group's land."""
    names = [n for n in preg if ISLANDS.get(n, "mainland") == group]
    half = block * RES / 2
    s = np.linspace(-half * .92, half * .92, 7)
    cells, series = [], []
    for by in range(daily.shape[0] // block):
        for bx in range(daily.shape[1] // block):
            lat = lat0 + (by * block + (block - 1) / 2) * RES
            lon = lon0 + (bx * block + (block - 1) / 2) * RES
            hits = {}
            for dy in s:
                for dx in s:
                    p = Point(lon + dx, lat + dy)
                    for n in names:
                        if preg[n].contains(p):
                            hits[n] = hits.get(n, 0) + 1
                            break
            frac = sum(hits.values()) / len(s) ** 2
            if frac < 0.06:
                continue
            cells.append({"lat": round(lat, 3), "lon": round(lon, 3), "frac": round(frac, 3),
                          "region": max(hits, key=hits.get), "group": group, "size": block * RES})
            series.append(daily[by * block:(by + 1) * block, bx * block:(bx + 1) * block].mean(axis=(0, 1)))
    return cells, series


def doy365(d):
    """Day index in a 365-day calendar; Feb 29 shares Feb 28's slot."""
    n = d.timetuple().tm_yday - 1
    leap = d.year % 4 == 0 and (d.year % 100 != 0 or d.year % 400 == 0)
    if leap and n >= 59:
        n -= 1
    return n


def main():
    os.makedirs(OUT, exist_ok=True)
    raw = {k: load_daily(*v[:1], *v[3:]) for k, v in BOXES.items()}
    nd = min(a.shape[2] for a in raw.values())
    start = dt.date(min(FIRST_YEAR, LAST_YEAR_FILE + 1), 1, 1)
    dates = [start + dt.timedelta(i) for i in range(nd)]
    print("days", nd, dates[0], "->", dates[-1])

    # --- cells + land masks --------------------------------------------------------
    regs = json.load(open(os.path.join(HERE, "geo", "regions_mask.json")))["regions"]
    preg = {r["name"]: prep(shape(r["geom"])) for r in regs}
    cells, series = [], []
    for group, box, block in GROUPS:
        _, lat0, lon0, _, _ = BOXES[box]
        c, sr = make_cells(group, raw[box][:, :, :nd], lat0, lon0, block, preg)
        cells += c; series += sr
        print(group, len(c), "cells")
    X = np.stack(series, axis=1).astype(np.float64)           # (days, cells)
    nc = X.shape[1]
    # area weight: land fraction x cell area (cos lat x size^2)
    w = np.array([c["frac"] * np.cos(np.radians(c["lat"])) * c["size"] ** 2 for c in cells])
    w /= w.sum()
    print("cells", nc)

    # --- daily binary ------------------------------------------------------------
    q = np.clip(np.rint(np.sqrt(X) * Q), 0, 255).astype(np.uint8)
    q.tofile(os.path.join(OUT, "daily.bin"))

    # --- climatology (1991-2020), 365-day calendar, 31-day circular smoothing -----------
    years = np.array([d.year for d in dates])
    doys = np.array([doy365(d) for d in dates])
    sel = (years >= NORMAL[0]) & (years <= NORMAL[1])
    clim = np.zeros((365, nc))
    cnt = np.zeros(365)
    np.add.at(clim, doys[sel], X[sel])
    np.add.at(cnt, doys[sel], 1)
    clim /= cnt[:, None]
    k = np.ones(31) / 31
    clim = np.stack([np.convolve(np.r_[c[-15:], c, c[:15]], k, mode="valid") for c in clim.T], axis=1)
    np.rint(clim * 100).astype(np.uint16).tofile(os.path.join(OUT, "clim.bin"))

    # --- yearly totals (exact) -----------------------------------------------------
    ylist = sorted(set(years.tolist()))
    yearly = np.stack([X[years == y].sum(0) for y in ylist])  # (years, cells)
    normal_cell = yearly[[ylist.index(y) for y in range(NORMAL[0], NORMAL[1] + 1)]].mean(0)
    nat = yearly @ w
    nat_normal = float(normal_cell @ w)
    last_complete = ylist[-1] if dates[-1] == dt.date(ylist[-1], 12, 31) else ylist[-2]
    full = [y for y in ylist if y <= last_complete]
    fy = np.array(full); fv = nat[:len(full)]
    slope = np.polyfit(fy, fv, 1)[0] * 10
    order = np.argsort(fv)
    half = len(full) // 2
    # same-period comparison for the partial current year
    partial = None
    if last_complete != ylist[-1]:
        cy = ylist[-1]
        m = (years == cy)
        doy_end = doys[m].max()
        ytd = float((X[m].sum(0)) @ w)
        norm_ytd = float(clim[:doy_end + 1].sum(0) @ w)
        partial = {"year": cy, "through": dates[-1].isoformat(), "mm": round(ytd, 1),
                   "normal_mm": round(norm_ytd, 1)}

    nat_daily = X @ w
    wettest_day = int(np.argmax(nat_daily))
    stats = {
        "first_year": full[0], "last_complete_year": last_complete,
        "normal_mm": round(nat_normal, 1),
        "trend_mm_per_decade": round(float(slope), 1),
        "driest": [{"year": int(fy[i]), "mm": round(float(fv[i]), 1)} for i in order[:5]],
        "wettest": [{"year": int(fy[i]), "mm": round(float(fv[i]), 1)} for i in order[::-1][:5]],
        "first_half": {"from": full[0], "to": full[half - 1], "mm": round(float(fv[:half].mean()), 1)},
        "second_half": {"from": full[half], "to": full[-1], "mm": round(float(fv[half:].mean()), 1)},
        "dry_years_since_2000": int(sum(1 for i in order[:10] if fy[i] >= 2000)),
        "partial": partial,
        "wettest_day": {"date": dates[wettest_day].isoformat(), "mm": round(float(nat_daily[wettest_day]), 1)},
        "max_cell_day_mm": round(float(X.max()), 1),
    }
    # regional (autonomous community) yearly totals
    rnames = sorted(set(c["region"] for c in cells))
    reg_yearly = {}
    for rn in rnames:
        idx = [i for i, c in enumerate(cells) if c["region"] == rn]
        rw = w[idx] / w[idx].sum()
        reg_yearly[rn] = {"yearly": [round(float(v), 1) for v in yearly[:, idx] @ rw],
                          "normal": round(float(normal_cell[idx] @ rw), 1)}

    meta = {
        "source": "ERA5 hourly total precipitation (Copernicus C3S / ECMWF), 0.25 deg, "
                  "aggregated to daily UTC sums and 0.5 deg cells",
        "start": dates[0].isoformat(), "end": dates[-1].isoformat(), "days": nd,
        "q": Q, "normal_period": NORMAL,
        "cells": cells, "weights": [round(float(v), 6) for v in w],
        "years": ylist,
        "yearly": [[round(float(v)) for v in row] for row in yearly],
        "normal_cell": [round(float(v), 1) for v in normal_cell],
        "national_yearly": [round(float(v), 1) for v in nat],
        "regions": reg_yearly,
        "stats": stats,
    }
    json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), separators=(",", ":"), ensure_ascii=False)
    import base64  # artifact hosting only serves text types: ship binaries as base64 .txt
    for b in ("daily", "clim"):
        raw = open(os.path.join(OUT, b + ".bin"), "rb").read()
        open(os.path.join(OUT, b + ".b64.txt"), "w").write(base64.b64encode(raw).decode())
        os.remove(os.path.join(OUT, b + ".bin"))
    for f in ("meta.json", "daily.b64.txt", "clim.b64.txt"):
        print(f, os.path.getsize(os.path.join(OUT, f)) // 1024, "KB")
    print(json.dumps(stats, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
