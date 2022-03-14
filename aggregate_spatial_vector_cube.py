import xarray as xr
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.multipoint import MultiPoint
import json

input_raster_cube = xr.open_dataarray('input.nc')
input_raster_cube = input_raster_cube.rio.set_crs('EPSG:32632')
input_raster_cube

# Read the input geoJSON
input_vector_cube = gpd.read_file('urban_forest_points.geojson')
input_vector_cube

# Sample transformation of vector data from EPSG:4326 to EQUI7
# input_vector_cube = input_vector_cube.set_crs(4326)
# input_vector_cube_equi7 = input_vector_cube.to_crs("PROJCRS[\"Azimuthal_Equidistant\",BASEGEOGCRS[\"WGS 84\",DATUM[\"World Geodetic System 1984\",ELLIPSOID[\"WGS 84\",6378137,298.257223563,LENGTHUNIT[\"metre\",1]]],PRIMEM[\"Greenwich\",0,ANGLEUNIT[\"degree\",0.0174532925199433]],ID[\"EPSG\",4326]],CONVERSION[\"Modified Azimuthal Equidistant\",METHOD[\"Modified Azimuthal Equidistant\",ID[\"EPSG\",9832]],PARAMETER[\"Latitude of natural origin\",53,ANGLEUNIT[\"degree\",0.0174532925199433],ID[\"EPSG\",8801]],PARAMETER[\"Longitude of natural origin\",24,ANGLEUNIT[\"degree\",0.0174532925199433],ID[\"EPSG\",8802]],PARAMETER[\"False easting\",5837287.81977,LENGTHUNIT[\"metre\",1],ID[\"EPSG\",8806]],PARAMETER[\"False northing\",2121415.69617,LENGTHUNIT[\"metre\",1],ID[\"EPSG\",8807]]],CS[Cartesian,2],AXIS[\"easting\",east,ORDER[1],LENGTHUNIT[\"metre\",1,ID[\"EPSG\",9001]]],AXIS[\"northing\",north,ORDER[2],LENGTHUNIT[\"metre\",1,ID[\"EPSG\",9001]]]]")
# input_vector_cube_equi7

# The computation also stores information about the total count of pixels (valid + invalid pixels) and
# the number of valid pixels (see is_valid) for each geometry. These values are added as a new dimension with
# a dimension name derived from target_dimension by adding the suffix _meta. The new dimension has the
# dimension labels total_count and valid_count.


def aggregate_spatial(raster_cube,vector_cube,target_dimension='result'):
    input_raster_cube_dims = list(raster_cube.dims)
    if len(input_raster_cube_dims)>3:
        raise Exception('TooManyDimensions - The number of dimensions must be reduced to three for aggregate_spatial. Input raster-cube dimensions: {}'.format(input_raster_cube_dims))
    if 'x' in input_raster_cube_dims:
        input_raster_cube_dims.remove('x')
    if 'y' in input_raster_cube_dims:
        input_raster_cube_dims.remove('y')
    if len(input_raster_cube_dims) == 0:
        input_raster_cube_dims = ['result']
    bands_or_timesteps = None
    if input_raster_cube_dims[0] in list(raster_cube.dims):
        bands_or_timesteps = raster_cube[input_raster_cube_dims[0]].values
        
    # Case when geoJSON is provided
    if type(vector_cube) == dict:
        vector_cube = gpd.GeoDataFrame.from_features(vector_cube)
    # Case when a vector-cube (result of vector_to_random_points for instance) is provided
    elif type(vector_cube) == gpd.geodataframe.GeoDataFrame:
        pass
    else:
        raise Exception('[!] No compatible vector input data has been provided.')
    
    input_vector_cube_columns = list(vector_cube.columns)

        
    output_raster_cube_columns = input_vector_cube_columns + [target_dimension,target_dimension + '_meta']
    output_vector_cube = gpd.GeoDataFrame(columns = output_raster_cube_columns)
    # print(output_vector_cube)
    
    ## Input geometries are in EPSG:4326 and the data has a different projection. We reproject the vector-cube
    vector_cube = vector_cube.set_crs(4326)
    vector_cube_utm = vector_cube.to_crs(32632) #TODO: read the crs from raster-cube. For EODC this will likely be almost always EQUI7
    
    ## First clip the data keeping only the data within the polygons
    crop = raster_cube.rio.clip(vector_cube_utm.geometry, drop=True)
    reducer = 'mean' #TODO: handle multiple reducers
    
    ## Loop over the geometries in the FeatureCollection and apply the reducer
    geom_crop_list = []
    for i in range(len(vector_cube_utm)):
        if reducer == 'mean':
            geom_crop = crop.rio.clip(vector_cube_utm.loc[[i]].geometry)
            valid_data_dict = {}
            if bands_or_timesteps is not None:
                total_count = len(geom_crop.x) * len(geom_crop.y) * len(geom_crop[input_raster_cube_dims[0]])
            else:
                total_count = len(geom_crop.x) * len(geom_crop.y)
            invalid_count = np.isnan(geom_crop).sum().values
            
            valid_count = total_count - invalid_count
            valid_data_dict['total_count'] = float(total_count)
            valid_data_dict['valid_count'] = float(valid_count)
            
            reduced_value = geom_crop.mean(dim=['x','y'])
            # print(geom_crop)
            
            raster_data_dict = {}
            if bands_or_timesteps is not None:
                for b_or_t in bands_or_timesteps:
                    raster_data_dict[str(b_or_t)] = reduced_value.loc[{input_raster_cube_dims[0]:b_or_t}].item()
            else:
                raster_data_dict = reduced_value.item()

            vector_data_dict = {}
            vector_data_dict[target_dimension] = raster_data_dict
            vector_data_dict[target_dimension + '_meta'] = valid_data_dict
            for ic in input_vector_cube_columns:
                # print('ic',ic)
                # print(vector_cube.loc[[i],ic].item())
                # print(vector_cube.loc[[i]])
                vector_data_dict[ic] = vector_cube.loc[[i],ic].item()
            # print(vector_data_dict)
            output_vector_cube = output_vector_cube.append(vector_data_dict,ignore_index=True)

    return output_vector_cube

result = aggregate_spatial(input_raster_cube.mean(dim='variable').rio.set_crs('EPSG:32632').mean(dim='time'),input_vector_cube)
result.to_file('sample_output.geojson', driver='GeoJSON') 

result = aggregate_spatial(input_raster_cube.mean(dim='time').rio.set_crs('EPSG:32632'),input_vector_cube)
result.to_file('sample_output_with_bands.geojson', driver='GeoJSON')

result = aggregate_spatial(input_raster_cube.mean(dim='variable').rio.set_crs('EPSG:32632'),input_vector_cube)
result.to_file('sample_output_with_time.geojson', driver='GeoJSON') 
