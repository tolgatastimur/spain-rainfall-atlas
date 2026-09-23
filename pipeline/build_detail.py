"""
Per-community detail files for the drill-down views (site/data/detail/<slug>.json).

Each file holds the community outline, its sub-units (provinces; islands for the
Balearics) and the native 0.25 deg ERA5 points that overlap it:
  points   lat, lon, f (land fraction in the community), u (sub-unit fractions)
  yearly   [year][point] totals (mm), normal [point] 1991-2020 annual mean
  daily    base64 uint8 [day][point], q = round(sqrt(mm) * Q)   (same encoding as daily.b64.txt)
  clim     base64 uint16 [365][point] daily normal (31-day smoothed) * 100
  munis    municipalities of each sub-unit (or of the community when it has a single
           province): name, outline, and gaussian weights to the points at its centre

Run after fetch_era5.py (both boxes) and build_geo.py.
"""
import os, json, base64, re, unicodedata, datetime as dt
import numpy as np
from shapely.geometry import shape, Point, Polygon
from shapely.prepared import prep
from shapely.ops import unary_union
import build_data as bd
from build_geo import rings, label_point

HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "..", "site", "data", "detail")
SKIP = {"Ceuta", "Melilla"}
# the Balearics are one province; use the islands as sub-units instead
BALEARIC_ISLANDS = [("Mallorca", 2.95, 39.6), ("Menorca", 4.1, 39.95), ("Eivissa", 1.4, 38.98), ("Formentera", 1.47, 38.69)]


def slug(name):
    s = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")


SIGMA = 0.16      # same smoothing the site uses for 0.25 deg points


def decode_topo(topo, obj):
    """TopoJSON (quantized, delta-encoded arcs) -> list of (properties, shapely geometry)."""
    sx, sy = topo["transform"]["scale"]; tx, ty = topo["transform"]["translate"]
    arcs = []
    for arc in topo["arcs"]:
        x = y = 0; pts = []
        for dx, dy in arc:
            x += dx; y += dy; pts.append((x * sx + tx, y * sy + ty))
        arcs.append(pts)

    def ring(idx):
        pts = []
        for i in idx:
            a = arcs[i] if i >= 0 else arcs[~i][::-1]
            pts.extend(a if not pts else a[1:])
        return pts
    out = []
    for g in topo["objects"][obj]["geometries"]:
        polys = g["arcs"] if g["type"] == "MultiPolygon" else [g["arcs"]] if g["type"] == "Polygon" else []
        ps = [Polygon(ring(p[0]), [r for r in map(ring, p[1:]) if len(r) > 3]) for p in polys if len(ring(p[0])) > 3]
        if ps:
            out.append((dict(g.get("properties", {}), id=g.get("id")), unary_union([q.buffer(0) for q in ps])))
    return out


def muni_entries(munis, container, pts):
    """Municipalities inside `container` with outline + weights to the community's points."""
    cos = np.cos(np.radians(container.centroid.y))
    P = np.array([[p["lon"], p["lat"]] for p in pts])
    res = []
    for props, g in munis:
        rp = g.representative_point()
        if not container.contains(rp):
            continue
        d2 = ((P[:, 0] - rp.x) * cos) ** 2 + (P[:, 1] - rp.y) ** 2
        w = np.exp(-d2 / (2 * SIGMA ** 2))
        keep = np.where(w > .01)[0]
        if not len(keep):
            keep = np.array([int(np.argmin(d2))]); w = np.ones(len(P))
        ws = w[keep] / w[keep].sum()
        res.append({"n": props.get("name", ""), "r": rings(g, 0.0015, 0.0), "w": [[int(c), round(float(x), 3)] for c, x in zip(keep, ws)]})
    return res


def clim_of(X, years, doys):
    sel = (years >= bd.NORMAL[0]) & (years <= bd.NORMAL[1])
    clim = np.zeros((365, X.shape[1])); cnt = np.zeros(365)
    np.add.at(clim, doys[sel], X[sel]); np.add.at(cnt, doys[sel], 1)
    clim /= cnt[:, None]
    k = np.ones(31) / 31
    return np.stack([np.convolve(np.r_[c[-15:], c, c[:15]], k, mode="valid") for c in clim.T], axis=1)


def main():
    os.makedirs(OUT, exist_ok=True)
    comm = json.load(open(os.path.join(HERE, "geo", "communities.geojson")))["features"]
    prov = json.load(open(os.path.join(HERE, "geo", "provinces.geojson")))["features"]
    munis = decode_topo(json.load(open(os.path.join(HERE, "geo", "municipalities.topo.json"))), "municipalities")
    print("municipalities", len(munis))
    raw = {k: bd.load_daily(*v[:1], *v[3:]) for k, v in bd.BOXES.items()}
    nd = min(a.shape[2] for a in raw.values())
    start = dt.date(bd.FIRST_YEAR, 1, 1)
    dates = [start + dt.timedelta(i) for i in range(nd)]
    years = np.array([d.year for d in dates]); doys = np.array([bd.doy365(d) for d in dates])
    ylist = sorted(set(years.tolist()))
    s = np.linspace(-0.115, 0.115, 7)
    index = []

    for f in comm:
        name = f["properties"]["name"]
        if name in SKIP:
            continue
        cgeom = shape(f["geometry"]).buffer(0)
        # sub-units
        if name == "Baleares":
            polys = list(cgeom.geoms) if hasattr(cgeom, "geoms") else [cgeom]
            units = []
            for uname, lon, lat in BALEARIC_ISLANDS:
                g = unary_union([p for p in polys if p.distance(Point(lon, lat)) < 0.12])
                units.append({"name": uname, "geom": g})
        else:
            # the two source files number communities differently: match provinces spatially
            units = [{"name": p["properties"]["name"], "geom": g} for p in prov
                     for g in [shape(p["geometry"]).buffer(0)] if cgeom.contains(g.representative_point())]
        if len(units) < 2:
            units = []
        box = "canarias" if name == "Canarias" else "mainland"
        _, lat0, lon0, nlat, nlon = bd.BOXES[box]
        data = raw[box]
        pc, pu = prep(cgeom), [prep(u["geom"]) for u in units]
        minx, miny, maxx, maxy = cgeom.bounds
        pts, series = [], []
        for j in range(nlat):
            lat = lat0 + j * bd.RES
            if lat < miny - .2 or lat > maxy + .2:
                continue
            for i in range(nlon):
                lon = lon0 + i * bd.RES
                if lon < minx - .2 or lon > maxx + .2:
                    continue
                hit, uh = 0, [0] * len(units)
                for dy in s:
                    for dx in s:
                        p = Point(lon + dx, lat + dy)
                        if pc.contains(p):
                            hit += 1
                            for k, u in enumerate(pu):
                                if u.contains(p):
                                    uh[k] += 1
                                    break
                frac = hit / len(s) ** 2
                if frac < 0.06:
                    continue
                pts.append({"lat": round(lat, 3), "lon": round(lon, 3), "f": round(frac, 3),
                            "u": [[k, round(n / len(s) ** 2, 3)] for k, n in enumerate(uh) if n]})
                series.append(data[j, i, :nd])
        X = np.stack(series, axis=1).astype(np.float64)
        q = np.clip(np.rint(np.sqrt(X) * bd.Q), 0, 255).astype(np.uint8)
        clim = clim_of(X, years, doys)
        yearly = np.stack([X[years == y].sum(0) for y in ylist])
        normal = yearly[[ylist.index(y) for y in range(bd.NORMAL[0], bd.NORMAL[1] + 1)]].mean(0)
        out = {
            "name": name, "q": bd.Q, "start": dates[0].isoformat(), "days": nd, "years": ylist,
            "rings": rings(cgeom, 0.004, 0.0004),
            "units": [{"name": u["name"], "rings": rings(u["geom"], 0.004, 0.0004), "label": label_point(u["geom"]),
                       "munis": muni_entries(munis, u["geom"], pts)} for u in units],
            "munis": [] if units else muni_entries(munis, cgeom, pts),
            "points": pts,
            "yearly": [[round(float(v)) for v in row] for row in yearly],
            "normal": [round(float(v), 1) for v in normal],
            "daily": base64.b64encode(q.tobytes()).decode(),
            "clim": base64.b64encode(np.rint(clim * 100).astype(np.uint16).tobytes()).decode(),
        }
        fn = slug(name) + ".json"
        json.dump(out, open(os.path.join(OUT, fn), "w"), separators=(",", ":"), ensure_ascii=False)
        nm = len(out["munis"]) + sum(len(u["munis"]) for u in out["units"])
        index.append({"name": name, "file": fn, "units": len(units), "points": len(pts)})
        print(f"{name:20s} {len(pts):4d} pts {len(units):2d} units {nm:4d} munis {os.path.getsize(os.path.join(OUT, fn)) // 1024:6d} KB")
    json.dump(index, open(os.path.join(OUT, "index.json"), "w"), ensure_ascii=False, indent=0)


if __name__ == "__main__":
    main()
