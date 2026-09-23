"""Build simplified map geometry for the site: mainland Spain, its autonomous
communities (borders + label points) and the Balearic and Canary Islands."""
import json, os
from shapely.geometry import shape, mapping, Polygon
from shapely.ops import unary_union

HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "..", "site", "data")


def rings(geom, tol, min_area):
    geom = geom.simplify(tol, preserve_topology=True)
    polys = [geom] if isinstance(geom, Polygon) else list(getattr(geom, "geoms", []))
    res = []
    for p in polys:
        if p.area < min_area or p.is_empty:
            continue
        res.append([[round(x, 3), round(y, 3)] for x, y in p.exterior.coords])
    return res


ISLANDS = {"Baleares": "baleares", "Canarias": "canarias"}


def label_point(g):
    """Point for the region's name: inside its largest polygon."""
    polys = [g] if isinstance(g, Polygon) else list(g.geoms)
    big = max(polys, key=lambda p: p.area)
    c = big.centroid
    p = c if big.contains(c) else big.representative_point()
    return [round(p.x, 3), round(p.y, 3)]


def main():
    os.makedirs(OUT, exist_ok=True)
    comm = json.load(open(os.path.join(HERE, "geo", "communities.geojson")))
    regions, islands = [], {}
    for f in comm["features"]:
        name = f["properties"]["name"]
        if name in ("Ceuta", "Melilla"):
            continue
        g = shape(f["geometry"]).buffer(0)
        if name in ISLANDS:
            islands[ISLANDS[name]] = {"name": name, "geom": g}
        else:
            regions.append({"name": name, "geom": g})
    mainland = unary_union([r["geom"] for r in regions])

    geo = {
        "spain": rings(mainland, 0.012, 0.004),
        "regions": [{"name": r["name"], "rings": rings(r["geom"], 0.012, 0.004), "label": label_point(r["geom"])}
                    for r in regions],
        "islands": {k: {"name": v["name"], "rings": rings(v["geom"], 0.004, 0.0015), "label": label_point(v["geom"])}
                    for k, v in islands.items()},
    }
    json.dump(geo, open(os.path.join(OUT, "geo.json"), "w"), separators=(",", ":"))
    # shapely-friendly copy (all communities incl. islands) for the grid masks
    allr = regions + list(islands.values())
    json.dump({"regions": [{"name": r["name"], "geom": mapping(r["geom"])} for r in allr]},
              open(os.path.join(HERE, "geo", "regions_mask.json"), "w"))
    print("geo.json", os.path.getsize(os.path.join(OUT, "geo.json")) // 1024, "KB",
          {k: len(v["rings"]) for k, v in geo["islands"].items()})


if __name__ == "__main__":
    main()
