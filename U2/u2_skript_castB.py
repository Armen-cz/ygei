import rasterio
from rasterio.transform import from_bounds
from rasterio.features import shapes
from pyproj import Transformer
import geopandas as gpd
from shapely.geometry import shape

"""Part 1 - creates georeferenced raster"""

# reads to image from matlab script
with rasterio.open("lesy.tif") as pic:
    arr = pic.read()  # shape = (bands, rows, cols)
    rows_num, cols_num = arr.shape[1], arr.shape[2]

# defines bounding box in geographic coordinates
xmin_deg = 14 + 30/60       
ymin_deg = 50 + 35/60           
xmax_deg = 14 + 37.5/60       
ymax_deg = 50 + 40/60         

# converts to cartesian coordinates (EPSG:28403)
transformer = Transformer.from_crs("EPSG:4326", "EPSG:28403", always_xy=True)
xmin_m, ymin_m = transformer.transform(xmin_deg, ymin_deg)
xmax_m, ymax_m = transformer.transform(xmax_deg, ymax_deg)

# return affine transform for the raster
affine_transform = from_bounds(xmin_m, ymin_m, xmax_m, ymax_m, cols_num, rows_num)

# writes proper georeferenced raster
output_path = "lesy_georef.tif"
with rasterio.open(
    output_path,
    "w",
    driver="GTiff",
    height=rows_num,
    width=cols_num,
    count=arr.shape[0],    
    dtype=arr.dtype,
    crs="EPSG:28403",
    transform=affine_transform,
) as gtif:
    gtif.write(arr)

print(f"georeferenced raster saved: {output_path}")


"""Part 2 - creates polygons from raster and saves it as Geopackage"""

# opens the created raster
with rasterio.open("lesy_georef.tif") as gtif:
    image = gtif.read(1)
    mask = image == 0 # creates mask for value 0
    transform = gtif.transform
    crs = gtif.crs

# creates polygons and saves their area to a new attribute
geoms = []
areas = []
for geom, val in shapes(image, mask=mask, transform=transform):
    geoms.append(shape(geom))
    areas.append(shape(geom).area)  

# creates GeoDataFrame
gdf = gpd.GeoDataFrame({"area": areas}, geometry=gpd.GeoSeries(geoms, crs=crs))

# saves to Geopackage
output_gpkg_file = "lesy_polygons.gpkg"
gdf.to_file(output_gpkg_file, layer="lesy_TM25_sk2", driver="GPKG")

print(f"geopackage succesfully saves: {output_gpkg_file}")

