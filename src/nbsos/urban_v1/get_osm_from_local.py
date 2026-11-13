import osmium
import shapely.geometry as geom
import geopandas as gpd

def tag_match(tags, filter_dict):
    """Return True if OSM tags match the filter dictionary."""
    for key, rule in filter_dict.items():
        if key not in tags:
            continue

        val = tags[key]

        if rule is True:
            return True

        if isinstance(rule, list):
            if val in rule:
                return True

        if isinstance(rule, str):
            if val == rule:
                return True

    return False


class GenericExtractor(osmium.SimpleHandler):
    def __init__(self, tag_filter):
        super().__init__()
        self.tag_filter = tag_filter
        self.points = []
        self.lines = []
        self.polygons = []

    def node(self, n):
        if tag_match(n.tags, self.tag_filter):
            point = geom.Point(n.lon, n.lat)
            self.points.append({
                "id": n.id,
                "geometry": point,
                **dict(n.tags)
            })

    def way(self, w):
        if not tag_match(w.tags, self.tag_filter):
            return

        coords = [(n.lon, n.lat) for n in w.nodes]
        if len(coords) < 2:
            return

        # Polygon if closed
        if coords[0] == coords[-1]:
            geometry = geom.Polygon(coords)
            self.polygons.append({
                "id": w.id,
                "geometry": geometry,
                **dict(w.tags)
            })
        else:
            geometry = geom.LineString(coords)
            self.lines.append({
                "id": w.id,
                "geometry": geometry,
                **dict(w.tags)
            })

    def relation(self, r):
        if not tag_match(r.tags, self.tag_filter):
            return

        # Only multipolygon relations are supported here
        if r.tags.get("type") != "multipolygon":
            return

        outers = []
        inners = []

        for member in r.members:
            if member.role == "outer" and isinstance(member, osmium.osm.Way):
                coords = [(n.lon, n.lat) for n in member.nodes]
                outers.append(coords)
            elif member.role == "inner" and isinstance(member, osmium.osm.Way):
                coords = [(n.lon, n.lat) for n in member.nodes]
                inners.append(coords)

        if outers:
            geometry = geom.MultiPolygon(
                [geom.Polygon(outer, inners)] if inners else
                [geom.Polygon(outer) for outer in outers]
            )

            self.polygons.append({
                "id": r.id,
                "geometry": geometry,
                **dict(r.tags)
            })


# -----------------------------
# Run extractor
# -----------------------------

tags_water = {
    'waterway': True,
    'water': True,
    'natural': ['water', 'wetland', 'bay']
}

handler = GenericExtractor(tags_water)
handler.apply_file("yourfile.osm.pbf")

# Convert to GeoDataFrames
gdf_points = gpd.GeoDataFrame(handler.points, geometry="geometry", crs="EPSG:4326")
gdf_lines = gpd.GeoDataFrame(handler.lines, geometry="geometry", crs="EPSG:4326")
gdf_polygons = gpd.GeoDataFrame(handler.polygons, geometry="geometry", crs="EPSG:4326")

# Export shapefiles
gdf_points.to_file("points.shp")
gdf_lines.to_file("lines.shp")
gdf_polygons.to_file("polygons.shp")