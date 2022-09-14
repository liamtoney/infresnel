from pyproj import CRS
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info


# Find UTM CRS of a (lat, lon) point (see https://gis.stackexchange.com/a/423614)
def _estimate_utm_crs(lat, lon, datum_name='WGS 84'):
    utm_crs_list = query_utm_crs_info(
        datum_name=datum_name,
        area_of_interest=AreaOfInterest(
            west_lon_degree=lon,
            south_lat_degree=lat,
            east_lon_degree=lon,
            north_lat_degree=lat,
        ),
    )
    return CRS.from_epsg(utm_crs_list[0].code)  # Taking first entry of list here!
