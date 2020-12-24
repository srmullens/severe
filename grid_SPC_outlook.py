#####################################################################
# Plots SPC Outlooks for the area you are concerned about.          #
#   Considering outlooks represent a forecast of probabilities,     #
#   outlooks are smoothed by 25mi at their borders. The goal is     #
#   for people to not worry about what color they are in when       #
#   they are near the border. Focus on preparing for storms!        #
#                                                                   #
# Code by: Stephen Mullens                                          #
# June 2020                                                         #
#                                                                   #
# Code inspired by:                                                 #
#   Communications:                                                 #
#       https://twitter.com/ounwcm/status/1130197459618156545?s=20  #
#       https://twitter.com/ounwcm/status/712703585218203648?s=20   #
#       https://twitter.com/ounwcm/status/591582647123480577?s=20   #
#       https://twitter.com/ounwcm/status/524932920037752834?s=20   #
#       https://twitter.com/ounwcm/status/451364902050213888?s=20   #
#   Twitter capability:                                             #
#       @betelbot and https://github.com/hippke/betelbot            #
#   Data plotting:                                                  #
#       https://github.com/frontogenesis/metpy/.ipynb_checkpoints/  #
#           SPC Outlooks-checkpoint.ipynb                           #
#           download_shapes-checkpoint.ipynb                        #
#                                                                   #
# See guiding map making principles at:                             #
#   https://github.com/srmullens/map_making_principles              #
#####################################################################

import requests
from contextlib import closing
import zipfile, tarfile
import os, glob, shutil
import numpy as np
import pandas as pd
import pyproj
import itertools
from functools import partial
from math import degrees, radians, cos, sin, asin, sqrt

import geopandas
import shapely.vectorized
import shapely.ops as sops
import shapely.geometry as sg
from shapely.geometry import Point,Polygon

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from cartopy import crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

# This is needed because of some error with cartopy and matplotlib axes 
# See also [Axes issue](https://github.com/SciTools/cartopy/issues/1120)
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

from datetime import datetime as dt, timedelta
from dateutil import tz
import time as t
from timezonefinder import TimezoneFinder

from twython import Twython



# What area do you want to plot? 'data', 'CONUS', 'Southeast', 'Florida'
where='data'

# What resolution do you want? 'low', 'mid', 'high'
grid_res = 'low'

# What SPC day do you want to plot?
#   Current Operational Data (int between 1-8)
plot_day = 1
#   Archived date and time: ['YYYYMMDD','HHMM',int between 1-8]
#   Archived date must be 2020, due to SPC shapefile formatting.
#plot_day = ['20200514','1730',2]

# What type of plot: 'exact' or 'smooth'
plot_type = 'smooth'

# Send tweet?
send_tweet = True

# Need a plot_day, smoothing, and plot_type override?
override = False                # Master 'override' flag
plot_type_override = False      # Independent of master 'override' flag.
get_average_override = False    # Smoothing

if override:
    send_tweet = False
    plot_type_override = False
    get_average_override = False



########################
# Make some functions. #
########################

# Clears all contents of folder specified in the argument
def clear_folder_contents(folder):

    folder = f'{folder}/'
    extensions = ['.zip','.dbf','.prj','.shp','.shx']
    last = ['-shp']
    image = ['.png']
    map_features = ['ne_','cou']
    delete = ['.DS_Store']

    for file in os.listdir(folder):
        file_path = os.path.join(folder, file)
        try:
            if os.path.isfile(file_path) and file[-4:] in image:
                shutil.move(file_path,f'images/{file}')
            elif os.path.isdir(file_path) and file[:3] in map_features:
                shutil.move(file_path,f'map_features/{file}')
            elif os.path.isfile(file_path) and file[-4:] in extensions:
                os.unlink(file_path)
            elif os.path.isdir(file_path) and file[-4:] in last:
                shutil.rmtree(file_path)
            elif os.path.isfile(file_path) and file in delete:
                os.unlink(file_path)
        except Exception as e:
            print(e)


# Downloads a zip file and extracts it in a target directory of choice
def download_zip_files(issuing_center,product,shapefiles):
    for issuing_center in shapefiles:
        for product in shapefiles[issuing_center]:
            download_zip_file(file_url = shapefiles[issuing_center][product], root_folder = issuing_center)

    # Move the zip file to an appropriate folder.
    for zip_file in glob.glob(f'{os.getcwd()}/*.zip'):
        head, tail = os.path.split(zip_file)
        shutil.move(zip_file, f'zips/{tail}')


# Download the specified ZIP file, extract everything, rename it.
def download_zip_file(file_url, root_folder):
    # Check the inputs are the appropriate data type
    if not isinstance(file_url,str):
        raise TypeError(f'file URL must by type string, not {type(file_url)}')
    if not isinstance(root_folder,str):
        raise TypeError(f'folder must by type string, not {type(root_folder)}')

    # Get the bits needed for processing.
    file = file_url.split('/')[-1]
    folder = file.split('.')[0]
    remove_chars = folder[8:-4]
    folder = folder.replace(remove_chars,'')

    # Download the ZIP file.
    with requests.get(file_url, stream=True) as r:
        with open(file, 'wb') as f:
            for chunk in r.iter_content(chunk_size=128):
                f.write(chunk)

    # Extract all the files from the ZIP file.
    with zipfile.ZipFile(file, 'r') as zip_ref:
        zip_ref.extractall(f'{root_folder}/{folder}')
        zip_ref.close()

    # Rename all the files to uniform naming system.
    mypath = f'{root_folder}/{folder}/'
    files = [os.path.join(mypath,f) for f in os.listdir(mypath)
                if os.path.isfile(os.path.join(mypath, f))]
    for j,file in enumerate(files):
        new_fname = file.replace(remove_chars,'')
        shutil.move(file,new_fname)


# Get the polygon for a categorical risk.
#   If no polygon exists, create one near the equator 
#   that won't overlap with anything else.
def get_polygons(gdf,loc):
    if len(gdf)>=loc+1:
        gdf_loc = gdf.iloc[loc]
        gdf_geom = gdf_loc['geometry']
        if gdf_geom.geom_type == 'MultiPolygon':
            polygons = [polygon for polygon in gdf_geom]
        else:
            polygons = [gdf_geom]
        return polygons
    else:
        polygon =  Polygon([ [0+2*loc, 0+2*loc],
                             [1+2*loc, 0+2*loc],
                             [1+2*loc, 1+2*loc],
                             [0+2*loc, 1+2*loc] ])
        polygons = [polygon]
        return polygons


# Determine if one polygon is inside the other.
#   Probably helps if polygons listed in order of increasing risk,
#   but that may not be strictly necessary.
def this_contains_that(*args):
    # Make one polygon out of the 2, 3, or 4 polygons in the list.
    if len(args)==2:
        union_boundary = args[1].union(args[0]).boundary
    elif len(args)==3:
        union_01 = args[1].union(args[0])
        union_boundary = args[2].union(union_01).boundary
    elif len(args)==4:
        union_01 = args[1].union(args[0])
        union_012 = union_01.union(args[2])
        union_boundary = args[3].union(union_012).boundary
    else:
        raise RuntimeWarning(f'Must contain 2 to 4 polygons. You gave {len(args)} polygons.')

    # If that produces one boundary...
    if union_boundary.geom_type == 'LineString':
        # ...is it in a simple circle?
        inside_upper = union_boundary.is_ring

    # If that produces two boundaries...
    elif union_boundary.geom_type == 'MultiLineString':
        ubx = []

        # ...make polygons of them
        for line in union_boundary:
            coordinates=[]
            for index, point in enumerate(line.coords):
                if index == 0: first_pt = point
                coordinates.append(point)
            coordinates.append(first_pt)

            if len(coordinates) >= 3:
                polygon = Polygon(coordinates)
                ubx.append(polygon)

        # ...and see if one is within the other.
        inside_or_outside = [ubx[-1].within(ubx[0]),
                             ubx[0].within(ubx[-1])]

        inside_upper = any(inside_or_outside)

    return inside_upper


# Return the extent of several polygons, making sure they are all included.
def combined_extent(*args):
    arg_bounds = [arg.bounds for arg in args]
    west = min([float(ab[0])-1 for ab in arg_bounds])
    south = min([float(ab[1])-1 for ab in arg_bounds])
    east = max([float(ab[2])+1 for ab in arg_bounds])
    north = max([float(ab[3])+1 for ab in arg_bounds])
    combined_bounds = [west,east,south,north]
    return combined_bounds


# Check to see if area of lower risk is way bigger than area of higher risk.
def size_check(l_risk,u_risk):
    # Want size of l_risk + u_risk.
    lu_risk = l_risk.union(u_risk)

    # Convert u_risk polygon to projected equal area coordinates (eac).
    u_eac = sops.transform(
                partial(
                    pyproj.transform,
                    pyproj.Proj(init='EPSG:4326'),
                    pyproj.Proj(proj='moll')),
                u_risk)

    # Compute area
    u_area = u_eac.area

    # Convert lu_risk polygon to projected equal area coordinates (eac).
    #   Albers Equal Area = aea
    lu_eac = sops.transform(
                partial(
                    pyproj.transform,
                    pyproj.Proj(init='EPSG:4326'),
                    pyproj.Proj(proj='moll')),
                lu_risk)
    # Compute area
    lu_area = lu_eac.area

    # Return true if lu_risk is 5x bigger area than u_risk.
    return True if lu_area/u_area*100 > 500 else False


# Take the list of extents, order them, and keep the unique ones.
def keep_unique_extents(extents_list,order):
    # Copy the original list, so it doesn't change.
    extents_copy = extents_list.copy()

    # Order the original list in highest to lowest risk.
    extents_ordered = [x for _,x in sorted(zip(order,extents_copy))]

    # Keep only the first instance of each item.
    unique_extents=[]
    for item in extents_ordered:
        count=0
        if item not in unique_extents:
            unique_extents.append(item)
            count+=1
        elif item in unique_extents: continue
        else: print('problem')

    return unique_extents


# Gets shapefile records from natural earth website
#   https://www.naturalearthdata.com/downloads/
def get_shapes_list(res,cat,name):
    # Look for file here.
    file_location = f'map_features/ne_{res}_{name}/ne_{res}_{name}.shp'

    # If it's not there, get it.
    if not os.path.isfile(file_location):
        shpfilename = shpreader.natural_earth(resolution=res,
                                        category=cat,
                                        name=name)

        orig_folder = shpfilename.split('/')[:-1]
        file = shpfilename.split('/')[-1]
        folder = file.split('.')[0]

        for shp in glob.glob(f'{"/".join(orig_folder)}/{folder}.*'):
            if not os.path.isdir(f'map_features/{folder}'):
                os.mkdir(f'map_features/{folder}')
            shutil.copy2(shp,f'map_features/{folder}/{shp.split("/")[-1]}')

    # With data in place, read the file...
    reader = shpreader.Reader(file_location)

    # ...and get the records in a list.
    list_to_return = list(reader.records())

    return list_to_return


# Trim datasets to right around the extent of the map for easier plotting.
def trim_to_extent(x,y,data,mask,extent):
    start_timer = dt.now()
    print('--> Trimming data')

    # How many grids do you include for averaging?
    grids = dist_in_grids(x,extent)

    # Find the location of the extent lats/lons
    extent[0] = min(x[0], key=lambda x:abs(x-extent[0]))
    extent[1] = min(x[0], key=lambda x:abs(x-extent[1]))
    extent[2] = min(y[:,0], key=lambda x:abs(x-extent[2]))
    extent[3] = min(y[:,0], key=lambda x:abs(x-extent[3]))

    # Get the grid location of the extent plus the extra needed for averaging.
    west = int(np.where(x[0]==extent[0])[0][0]-grids)
    east = int(np.where(x[0]==extent[1])[0][0]+grids)
    south = int(np.where(y[:,0]==extent[2])[0][0]-grids)
    north = int(np.where(y[:,0]==extent[3])[0][0]+grids)

    # Make sure the numbers stay within the original grid.
    if west<0: west=0
    if south<0: south=0
    if east>x.shape[1]: east = x.shape[1]-1
    if north>x.shape[0]: north = x.shape[0]-1

    # Trim the data.
    x = x[south:north+1,west:east+1]
    y = y[south:north+1,west:east+1]
    data = data[south:north+1,west:east+1]
    mask = mask[south:north+1,west:east+1]

    time_elapsed = dt.now() - start_timer
    tsec = round(time_elapsed.total_seconds(),2)
    print(f'  --> Trimmed shape: {data.shape} ({tsec:.2f} seconds)')

    return x,y,data,mask,grids


# Calculate distance between lat/lon points.
def distance(lat1, lon1, lat2, lon2):

    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))

    # Radius of earth in kilometers is 6371
    km = 6371* c

    return km


# Get average value of places within 25 miles.
def get_average_values(x,y,data,mask,grids):

    start_timer = dt.now()
    print('--> Smoothing')

    i,j = x.shape
    smooth_data = np.zeros((i,j))
    last_percent = ''

    # For each grid cell in the data...
    for lat in range(i):
        for lon in range(j):

            # ...get its characteristics...
            sum_cats = count = 0
            orig_mask = mask[lat,lon]
            orig_value = data[lat,lon]
            orig_x = x[lat,lon]
            orig_y = y[lat,lon]

            # ...search through its neighboring cells...
            for look_x in range(i):
              if abs(look_x-lat)<=grids:
                for look_y in range(j):
                  if abs(look_y-lon)<=grids:

                    # ...and their characteristics.
                    include_mask = mask[look_x,look_y]
                    include_value = data[look_x,look_y]
                    include_x = x[look_x,look_y]
                    include_y = y[look_x,look_y]

                    # If the neighbor is close enough...
                    if include_mask:
                        dist = distance(include_y,include_x,orig_y,orig_x)
                        if dist < 40.2:     # 25 miles = 40.23 km

                            # ...include it in the smoothing around this location.
                            sum_cats += include_value
                            count += 1

            # Record the average category of all the nearby locations.
            if count>0: smooth_data[lat,lon] = sum_cats/count
            else: smooth_data[lat,lon] = orig_value

            # Show loop progress
            percent = round(lat/i*100,0)
            if percent%10==0 and lon==0 and percent!=last_percent:
                print(f'  --> {percent:.0f}%')
                last_percent = percent

    # Output results
    """
    for i in range(7):
        if i==0: print(f'  {i}: {(smooth_data == i).sum()}')
        else:
            print(f'{i-1}->{i}: {(smooth_data < i).sum() - (smooth_data <= i-1).sum()}')
            print(f'  {i}: {(smooth_data == i).sum()}')
    """

    time_elapsed = dt.now() - start_timer
    tsec = round(time_elapsed.total_seconds(),2)
    print(f'--> Smoothed ({tsec:.2f} seconds)')

    return smooth_data


# Determine an ideal legend location.
def legend_location(category,mask):
    print('--> Determining Legend Location')

    # Get items for slicing the array
    height,width = category.shape
    dh = int(round(height*0.33,0))
    dw = int(round(width*0.33,0))

    # Mask array to US coast
    ma_category = np.ma.masked_where(np.logical_not(mask),category)

    # Max category within the slice.
    upper_right_max = np.amax(ma_category[height-dh:height,width-dw:width])
    upper_left_max = np.amax(ma_category[height-dh:height,0:dw])
    lower_right_max = np.amax(ma_category[0:dh,width-dw:width])
    lower_left_max = np.amax(ma_category[0:dh,0:dw])
    max_corners = [lower_right_max,lower_left_max,upper_right_max,upper_left_max]

    # Sum of category values across the slice.
    upper_right_sum = np.sum(ma_category[height-dh:height,width-dw:width])
    upper_left_sum = np.sum(ma_category[height-dh:height,0:dw])
    lower_right_sum = np.sum(ma_category[0:dh,width-dw:width])
    lower_left_sum = np.sum(ma_category[0:dh,0:dw])
    sum_corners = [lower_right_sum,lower_left_sum,upper_right_sum,upper_left_sum]

    print(f'    --> lower right: Max: {lower_right_max:.1f}  Sum: {lower_right_sum:.2f}')
    print(f'    --> lower left: Max: {lower_left_max:.1f}  Sum: {lower_left_sum:.2f}')
    print(f'    --> upper right: Max: {upper_right_max:.1f}  Sum: {upper_right_sum:.2f}')
    print(f'    --> upper left: Max: {upper_left_max:.1f}  Sum: {upper_left_sum:.2f}')

    # List of legend locations corresponding to the corners in max_corners and sum_corners.
    corners = [4,3,1,2]

    # Sort legend locations by max category at any pixel. Lowest to highest.
    #   Do the same with sum_corners.
    corners_sort_by_low_max = [c for _,c in sorted(zip(max_corners,corners), key=lambda pair: pair[0])]
    sum_corners = [c for _,c in sorted(zip(max_corners,sum_corners), key=lambda pair: pair[0])]

    # Trim lists to locations with equivalently lowest max category at each pixel.
    #   Maybe there are two corners that max out at the second category.
    corners_sort_by_low_max = corners_sort_by_low_max[:max_corners.count(max_corners[0])]
    sum_corners = sum_corners[:max_corners.count(max_corners[0])]

    # Sort remaining corner locations by sum of categories across all pixels.
    leg_loc = [ll for _,ll in sorted(zip(sum_corners,corners_sort_by_low_max), key=lambda pair: pair[0])][0]

    print(f' --> Legend location: {leg_loc}')

    return leg_loc


# Convert time zones for printing
#   start_time gets used in the title on the map.
#   end_time isn't used right now, but it's included here just in case.
#   issue_time gets used in the tweet.
def convert_datetime_from_spc_to_local(polygon,start_time,end_time,issue_time,where):

    # Convert all times to UTC.
    from_zone = tz.gettz('UTC')
    start_utc_time = dt.strptime(start_time, '%Y%m%d%H%M').replace(tzinfo=from_zone)
    end_utc_time = dt.strptime(end_time, '%Y%m%d%H%M').replace(tzinfo=from_zone)
    issue_utc_time = dt.strptime(issue_time, '%Y%m%d%H%M').replace(tzinfo=from_zone)

    # Create some lists.
    new_zones_list = []
    start_zones_list = []
    end_zones_list = []
    issue_zones_list = []

    # Determine time zones by where the data is...
    if where in ['data','CONUS']:
        # Set up getting the time zone.
        tf = TimezoneFinder()

        # Extract a list of coordinates from the polygon(s).
        try:
            if polygon.geom_type == 'MultiPolygon':
                multipolygon = polygon
                coords = [point for polygon in multipolygon for point in polygon.exterior.coords[:-1]]
            elif polygon.geom_type == 'Polygon':
                coords = list(polygon.exterior.coords)
        except:
            coords = None
            new_zones_list=['UTC']

        # List the timezone at each polygon vertex.
        if coords is not None:
            for coord in coords:
                # What timezone is each coordinate of the highest risk in?
                new_zones_list.append(tf.timezone_at(lng=coord[0], lat=coord[1]))

        # Remove None values.
        new_zones_list = [i for i in new_zones_list if i]

        # Standardize time zones
        for i,item in enumerate(new_zones_list):
            zone = f'{start_utc_time.astimezone(tz.gettz(item)):%Z}'
            print(f'  --> {item}, {zone}')
            if zone in ['EST','EDT','-05']:
                new_zones_list[i] = 'America/New_York'
            elif zone in ['CST','CDT']:
                new_zones_list[i] = 'America/Chicago'
            elif zone in ['MST','MDT']:
                new_zones_list[i] = 'America/Denver'
            elif zone in ['PST','PDT']:
                new_zones_list[i] = 'America/Los_Angeles'

        # Reduce list to unique time zones.
        print(f'  --> All time zones: {new_zones_list}')
        new_zones_list = list(set(new_zones_list))

        # Sort the resulting list from west to east.
        west_to_east = ['America/Los_Angeles','America/Denver','America/Chicago','America/New_York','UTC']
        sort_by = []
        print(f'  --> Unique time zones: {new_zones_list}')
        for item in new_zones_list: sort_by.append(west_to_east.index(item))
        new_zones_list = [nzl for _,nzl in sorted(zip(sort_by,new_zones_list), key=lambda pair: pair[0])]

        print(f'  --> Sorted time zones: {new_zones_list}')

        # Use str timezone names to modify datetime objects
        for i,item in enumerate(new_zones_list):
            start_zones_list.append(start_utc_time.astimezone(tz.gettz(item)))
            end_zones_list.append(end_utc_time.astimezone(tz.gettz(item)))
            issue_zones_list.append(issue_utc_time.astimezone(tz.gettz(item)))

    # Use Eastern time zone for Florida or Southeast.
    else:
        to_zone = tz.gettz('America/New_York')
        start_zones_list = [start_utc_time.astimezone(tz.gettz(item))]
        end_zones_list = [end_utc_time.astimezone(tz.gettz(item))]
        issue_zones_list = [issue_utc_time.astimezone(tz.gettz(item))]


    # Generate string outputs based on how many time zones the highest risk covers.
    s1    = f'{start_zones_list[0]:%a, %b %d, %Y %-I:%M}'
    s1_pm = f'{start_zones_list[0].strftime("%p").lower()}'
    s1_tz = f'{start_zones_list[0]:%Z}'
    s2    = f'{start_zones_list[-1]:%-I:%M}'
    s2_pm = f'{start_zones_list[-1].strftime("%p").lower()}'
    s2_tz = f'{start_zones_list[-1]:%Z}'
    start_time = f'{s1}{s1_pm} {s1_tz}, {s2}{s2_pm} {s2_tz}'

    e1    = f'{end_zones_list[0]:%-I:%M}'
    e1_pm = f'{end_zones_list[0].strftime("%p").lower()}'
    e1_tz = f'{end_zones_list[0]:%Z}'
    e2    = f'{end_zones_list[-1]:%-I:%M}'
    e2_pm = f'{end_zones_list[-1].strftime("%p").lower()}'
    e2_tz = f'{end_zones_list[-1]:%Z}'

    i1    = f'{issue_zones_list[0]:%-I:%M}'
    i1_pm = f'{issue_zones_list[0].strftime("%p").lower()}'
    i1_tz = f'{issue_zones_list[0]:%Z}'
    i2    = f'{issue_zones_list[-1]:%-I:%M}'
    i2_pm = f'{issue_zones_list[-1].strftime("%p").lower()}'
    i2_tz = f'{issue_zones_list[-1]:%Z}'

    if len(start_zones_list)>=2:
        start_time = f'{s1}{s1_pm} {s1_tz}, {s2}{s2_pm} {s2_tz}'
        end_time = f'{e1}{e1_pm} {e1_tz}, {e2}{e2_pm} {e2_tz}'
        issue_time = f'{i1}{i1_pm} {i1_tz}, or {i2}{i2_pm} {i2_tz}'

    elif len(start_zones_list)==1:
        start_time = f'{s1}{s1_pm} {s1_tz}'
        end_time = f'{e1}{e1_pm} {e1_tz}'
        issue_time = f'{i1}{i1_pm} {i1_tz}'

    print(f'  --> start_time: {start_time}')
    print(f'  --> end_time: {end_time}')
    print(f'  --> issue_time: {issue_time}')

    return start_time, end_time, issue_time


# Calculate how many grids 25km is
def dist_in_grids(x,extent):

    resolution = round(x[0,1]-x[0,0],2)

    # Largest distance is determined by northernmost latitude
    lat = extent[3]

    # Calculate the distance in degrees.
    #   Source: https://stackoverflow.com/questions/3182260/python-geocode-filtering-by-distance
    #   Derived from: http://janmatuschek.de/LatitudeLongitudeBoundingCoordinates
    #   WARNING: problems if North/South Pole is in circle of interest
    #   WARNING: problems if longitude meridian +/-180 degrees intersects circle of interest
    distance_km = 25
    Earth_radius_km = 6371.0
    RADIUS = Earth_radius_km
    dlat = distance_km / RADIUS
    dlon = asin(sin(dlat) / cos(radians(lat)))

    # Conver the degrees to grid points.
    dy = round(degrees(dlat) / resolution * 10,0)
    dx = round(degrees(dlon) / resolution * 10,0)
    grids = max(dx,dy)
    print(f'  --> Res: {resolution}')
    print(f'  --> Grids in 25mi: {int(grids)} ({int(dx)} vs {int(dy)})')

    return grids


# Get the shapes of US counties.
def get_counties(data_crs):
    file_location = './map_features/countyp010g/countyp010g.shp'
    data_file = 'https://prd-tnm.s3.amazonaws.com/StagedProducts/Small-scale/data/Boundaries/countyp010g.shp_nt00934.tar.gz'
    # If there is no local file, download it.
    if not os.path.isfile(file_location):
        file = data_file.split('/')[-1]
        folder = file.split('.')[0]
        with requests.get(data_file, stream=True) as r:
            with open(f'map_features/{folder}/{file}', 'wb') as f:
                for chunk in r.iter_content(chunk_size=128):
                    f.write(chunk)

        with tarfile.open(file, mode='r:gz') as tar_ref:
            tar_ref.extractall(f'map_features/{folder}')
            tar_ref.close()

        if os.path.isfile(f'map_features/{folder}/{file}'):
            os.unlink(f'map_features/{folder}/{file}')
    # Read local file and get county geometries.
    reader = shpreader.Reader(file_location)
    counties = list(reader.geometries())
    # Make a shapely feature to plot.
    COUNTIES = cfeature.ShapelyFeature(counties, data_crs)

    return COUNTIES


# Get the country outlines of Mexico and Canada for the map.
def Mexico_Canada():
    # Get country shapes.
    countries = get_shapes_list('50m','cultural','admin_0_countries')
    # Get Mexico and Canada's shapes.
    country_list = ['Mexico','Canada']
    MC = [country for country in countries if country.attributes['NAME_EN'] in country_list]
    # Make a list of the geometries for these countries.
    Mex_Can = [mc.geometry for mc in MC]

    return Mex_Can


# Get the outlines of the Great Lakes for the map.
def Great_Lakes():
    # Get lake shapes.
    all_lakes = get_shapes_list('50m','physical','lakes')
    # Get Great Lake's shapes.
    list_of_lakes = ['Lake Superior','Lake Huron','Lake Michigan','Lake Erie','Lake Ontario']
    GL = [lake for lake in all_lakes if lake.attributes['name'] in list_of_lakes]
    # Make a list of teh geometries for these lakes.
    lakes = [gl.geometry for gl in GL]

    return lakes


# Tweet the results.
def tweet(text, image, send_tweet, reply):
    if send_tweet:
        consumer_key = os.environ.get('consumer_key')
        consumer_secret = os.environ.get('consumer_secret')
        access_token = os.environ.get('access_token')
        access_token_secret = os.environ.get('access_token_secret')

        print('  --> Tweeting...')
        twitter = Twython(consumer_key, consumer_secret, access_token, access_token_secret)

        # Tweet status.
        if not reply:
            response = twitter.upload_media(media=open(image, 'rb'))
            twitter.update_status(status=text, media_ids=[response['media_id']])

        # If there is another image, make or continue a thread.
        elif reply:
            # ...Get most recent tweet's ID from the timeline...
            timeline = twitter.get_user_timeline(screen_name='SmoothedPC',count=5)
            tweet_list = []
            for tweet in timeline:
                created = dt.strptime(tweet['created_at'],'%a %b %d %H:%M:%S %z %Y')
                tweet_list.append({'created_at':created,'id':tweet['id']})
                tweet_list = sorted(tweet_list, key = lambda i: i['created_at'],reverse=True)
                tweet_id = tweet_list[0]['id']

            response = twitter.upload_media(media=open(image, 'rb'))
            twitter.update_status(status=text,
                                media_ids=[response['media_id']],
                                in_reply_to_status_id=tweet_id,
                                auto_populate_reply_metadata=True)
        print('  --> Tweeted.')
    else:
        if not reply: print('    --> TEST Original tweet')
        elif reply: print('    --> TEST Reply to tweet')



########################################
# The main function to make the plots. #
########################################
def get_SPC_data(where,plot_type,plot_type_override,plot_day,grid_res,override):
    #
    # STEP 1: Get the files
    #

    print('--> Get the files')

    if isinstance(plot_day,int):
        print('--> Getting operational data')
        if plot_day<4:
            shapefiles = {
                'spc': {f'categorical_day{plot_day}': f'https://www.spc.noaa.gov/products/outlook/day{plot_day}otlk-shp.zip'}
                        }
        else:
            shapefiles = {
                'spc': {f'categorical_day{plot_day}': f'https://www.spc.noaa.gov/products/exper/day4-8/day{plot_day}prob-shp.zip'}
                        }
    elif isinstance(plot_day,list):
        archive_ymd = plot_day[0]
        archive_hm = plot_day[1]
        plot_day = plot_day[2]
        print('--> Getting archive: {archive_ymd}_{archive_hm}')
        shapefiles = {
            'spc': {f'categorical_day{plot_day}': f'https://www.spc.noaa.gov/products/outlook/archive/{archive_ymd[:4]}/day{plot_day}otlk_{archive_ymd}_{archive_hm}-shp.zip'} }

    """
    shapefiles = {
        'spc': {
            'categorical_day1': 'https://www.spc.noaa.gov/products/outlook/day1otlk-shp.zip',
            'categorical_day2': 'https://www.spc.noaa.gov/products/outlook/day2otlk-shp.zip',
            'categorical_day3': 'https://www.spc.noaa.gov/products/outlook/day3otlk-shp.zip',
            'categorical_day4': 'https://www.spc.noaa.gov/products/exper/day4-8/day4prob-shp.zip',
            'categorical_day5': 'https://www.spc.noaa.gov/products/exper/day4-8/day5prob-shp.zip',
            'categorical_day6': 'https://www.spc.noaa.gov/products/exper/day4-8/day6prob-shp.zip',
            'categorical_day7': 'https://www.spc.noaa.gov/products/exper/day4-8/day7prob-shp.zip',
            'categorical_day8': 'https://www.spc.noaa.gov/products/exper/day4-8/day8prob-shp.zip'
        },
        'wpc': {
            'excessive_rain_day1': 'https://ftp.wpc.ncep.noaa.gov/shapefiles/qpf/excessive/EXCESSIVERAIN_Day1_latest.zip',
            'excessive_rain_day2' : 'https://ftp.wpc.ncep.noaa.gov/shapefiles/qpf/excessive/EXCESSIVERAIN_Day2_latest.zip',
            'excessive_rain_day3' : 'https://ftp.wpc.ncep.noaa.gov/shapefiles/qpf/excessive/EXCESSIVERAIN_Day3_latest.zip'
        },
        'nhc': {
            'tropical_wx_outlook': 'https://www.nhc.noaa.gov/xgtwo/gtwo_shapefiles.zip'
        },
        'drought': {
            'latest': 'https://droughtmonitor.unl.edu/data/shapefiles_m/USDM_current_M.zip',
        }
    }
    """

    # Clear directories before getting data
    target_folders = ['spc', 'zips', 'images', 'map_features']
    #target_folders = ['nhc', 'spc', 'wpc', 'drought', 'zips']
    for folder in target_folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
        else:
            clear_folder_contents(folder)

    # Get the shapefiles
    print('--> Getting ZIPs')
    start_timer = dt.now()

    # Read in Shapefile
    #   If it's not available, wait and try again.
    tries = 1

    # For Day 1 and Day 3 forecasts...
    if plot_day<4:
        while tries<30:
            # Download ZIP files
            for issuing_center in shapefiles:
                for product in shapefiles[issuing_center]:
                    download_zip_files(issuing_center,product,shapefiles)
                # Read file you want
                cat_gdf = geopandas.read_file(f'{issuing_center}/day{plot_day}otlk-shp/day{plot_day}otlk_cat.shp')
                time_since_issued = dt.utcnow()-dt.strptime(cat_gdf['ISSUE'][0],'%Y%m%d%H%M')
                time_since_issued = time_since_issued.total_seconds()

                if time_since_issued > 9000 and override==False:  # 2.5 hours
                    print(f"  --> Not available yet. {time_since_issued:.0f}={dt.utcnow():%H%M} UTC - {dt.strptime(cat_gdf['ISSUE'][0],'%Y%m%d%H%M').strftime('%H%M')}")
                    if tries==29: print(f'Could not find day{plot_day}otlk_cat.shp after 15 minutes.'); return
                    else:
                        t.sleep(120)
                        tries += 1
                else:
                    print(f"  --> Got it! {time_since_issued:.0f}={dt.utcnow():%H%M} UTC - {dt.strptime(cat_gdf['ISSUE'][0],'%Y%m%d%H%M').strftime('%H%M')}")
                    tries+=30

    # For Day 4 through Day 8 forecasts...
    else:
        while tries<30:
            try:
                # Download ZIP files
                for issuing_center in shapefiles:
                    for product in shapefiles[issuing_center]:
                        download_zip_files(issuing_center,product,shapefiles)
                    # Read file you want
                    if plot_day==3:
                        cat_gdf = geopandas.read_file(f'{issuing_center}/day{plot_day}otlk-shp/day{plot_day}otlk_cat.shp')
                    else:
                        cat_gdf = geopandas.read_file(f'{issuing_center}/day{plot_day}prob-shp/day{plot_day}otlk_{start_timer:%Y%m%d}_prob.shp')
                    print(f"  --> Got it! {dt.utcnow():%H%M} UTC - {dt.strptime(cat_gdf['ISSUE'][0],'%Y%m%d%H%M').strftime('%H%M')}")
                    tries+=30
            except:
                print(f'  --> Not available yet. {dt.utcnow():%H%M} UTC')
                if tries==29: print(f'Could not find day{plot_day}otlk_{start_timer:%Y%m%d}_prob.shp after 15 minutes.'); return
                else:
                    t.sleep(120)
                    tries+=1

    # Delete data.
    for folder in target_folders[:-2]:
        print(folder)
        if os.path.exists(folder):
            clear_folder_contents(folder)
            # If no more files in the folder, delete the folder.
            print(folder,os.listdir(folder))
            if not len(os.listdir(folder)):
                shutil.rmtree(folder)


    # If there is no polygon, plot_nothing=True
    plot_nothing = len(cat_gdf)==1 and cat_gdf.loc[0].geometry is None

    # Plot low categories as exact polygons,
    #   regardless of the original setting.
    if plot_type_override==False:
        if plot_day<4 and cat_gdf.iloc[2:5].empty:
            plot_type = 'exact'
        elif plot_nothing:
            plot_type = 'exact'
        else: plot_type = 'smooth'


    time_elapsed = dt.now() - start_timer
    tsec = round(time_elapsed.total_seconds(),2)
    print(f'  --> Got ZIPs ({tsec:.2f} seconds)')


    #
    # STEP 2: Set the map extent, legend location, and coordinate reference
    #   system for plotting the map.
    #

    ### MAP EXTENT, LEGEND LOCATION, AND COORDINATE REFERENCE SYSTEM ###
    extents = False

    # Set Coordinate Reference System, extent, and legend location for the map
    if plot_nothing or where=='CONUS':
        extent = [-125,-66,24,50]
        leg_loc = 4
    elif where == 'Southeast':
        extent = [-95, -75, 24, 35]
        leg_loc = 4
    elif where=='Florida':
        extent = [-88,-79,24,32]
        leg_loc = 3
    elif where=='data':
        # Days 4-8, if there is any risk, calculate map extent.
        if len(cat_gdf)>=1 and plot_day>3:
            # Get polygon for lowest risk category.
            lowest_category = cat_gdf.iloc[0]['geometry']

            # Determine the extent.
            cat_extent = lowest_category.bounds
            extent = [float(cat_extent[0])-1,
                        float(cat_extent[2])+1,
                        float(cat_extent[1])-1,
                        float(cat_extent[3])+1 ]

            # NOTE: legend location is determined later.
            leg_loc = 0

        # Days 1-3, if there is at least a slight risk, calculate map extent.
        elif plot_day<4 and len(cat_gdf)>2:

            #0=TSTM; 1=MRGL; 2=SLGT; 3=ENH; 4=MDT; 5=HIGH

            new_extents=[]
            order=[]
            reply_polygons=[]
            reply_categories=[]

            # Get a list of polygons at each category from the shapefiles.
            HIGHs = get_polygons(cat_gdf,5)
            MDTs = get_polygons(cat_gdf,4)
            ENHs = get_polygons(cat_gdf,3)
            SLGTs = get_polygons(cat_gdf,2)
            MRGLs = get_polygons(cat_gdf,1)

            # Report at which levels polygons exist, and how many.
            no_poly_zone = Polygon([ [0, 0],[12, 0],[12, 12],[0, 12] ])
            nmrgl = f'{len(MRGLs)}' if any([not this_contains_that(MRGL,no_poly_zone) for MRGL in MRGLs]) else '-'
            nslgt = f'{len(SLGTs)}' if any([not this_contains_that(SLGT,no_poly_zone) for SLGT in SLGTs]) else '-'
            nenh = f'{len(ENHs)}' if any([not this_contains_that(ENH,no_poly_zone) for ENH in ENHs]) else '-'
            nmdt = f'{len(MDTs)}' if any([not this_contains_that(MDT,no_poly_zone) for MDT in MDTs]) else '-'
            nhigh = f'{len(HIGHs)}' if any([not this_contains_that(HIGH,no_poly_zone) for HIGH in HIGHs]) else '-'
            print(f'\nlengths: M, S, E, M, H')
            print(f' exists: {nmrgl}, {nslgt}, {nenh}, {nmdt}, {nhigh}')

            top_risk=None
            if nhigh!='-': top_risk='high'
            elif nmdt!='-': top_risk='mdt'
            elif nenh!='-': top_risk='enh'


            # Cycle through all the polygons, finding out what to plot.
            for r,MRGL in enumerate(MRGLs):
                for s,SLGT in enumerate(SLGTs):
                    for e,ENH in enumerate(ENHs):
                        for m,MDT in enumerate(MDTs):
                            for h,HIGH in enumerate(HIGHs):
                                one = this_contains_that(MDT,HIGH)
                                two = this_contains_that(ENH,MDT,HIGH)
                                three = this_contains_that(ENH,MDT)
                                four = this_contains_that(SLGT,ENH,MDT,HIGH)
                                five = this_contains_that(SLGT,ENH,MDT)
                                six = this_contains_that(SLGT,ENH)
                                seven = this_contains_that(MRGL,SLGT,ENH)
                                eight = this_contains_that(MRGL,SLGT)

                                high_exists = any([not this_contains_that(H,no_poly_zone) for H in HIGHs])
                                mdt_exists = any([not this_contains_that(M,no_poly_zone) for M in MDTs])
                                enh_exists = any([not this_contains_that(E,no_poly_zone) for E in ENHs])

                                ### ZOOMED IN ###

                                # if MDT contains HIGHs...
                                if one:
                                    # ...add MDT+HIGH extent
                                    print(f'    MDT {m+1} contains HIGH {h+1}')
                                    add_extent = combined_extent(MDT,HIGH)
                                    new_extents.append(add_extent)
                                    reply_polygons.append(HIGH)
                                    reply_categories.append('HIGH')
                                    order.append(0)
                                    # Check if MDT is way bigger than HIGH.
                                    if top_risk=='high':
                                        size_small = size_check(MDT,HIGH)
                                        if size_small and plot_day<3:
                                            # ...add HIGH extent
                                            print(f'    HIGH {h+1} much smaller than MDT {m+1}')
                                            add_extent = combined_extent(HIGH)
                                            new_extents.append(add_extent)
                                            reply_polygons.append(HIGH)
                                            reply_categories.append('HIGH')
                                            order.append(0.5)
                                # if ENH contains MDT & HIGH, but MDT doesn’t contain HIGH.
                                if not one and two and three:
                                    # ...add MDT extent
                                    print(f'    MDT {m+1} has no HIGH')
                                    add_extent = combined_extent(MDT)
                                    new_extents.append(add_extent)
                                    reply_polygons.append(MDT)
                                    reply_categories.append('MDT')
                                    order.append(1)
                                # if ENH doesn’t contain HIGH, and ENH contains MDT.
                                if not two and three:
                                    # ...add ENH+MDT extent
                                    print(f'    ENH {e+1} contains MDT {m+1}')
                                    add_extent = combined_extent(ENH,MDT)
                                    new_extents.append(add_extent)
                                    reply_polygons.append(MDT)
                                    reply_categories.append('MDT')
                                    order.append(2)
                                    # Check if ENH is way bigger than MDT.
                                    if top_risk=='mdt':
                                        size_small = size_check(ENH,MDT)
                                        if size_small and plot_day<3:
                                            # ...add MDT extent
                                            print(f'    MDT {m+1} much smaller than ENH {e+1}')
                                            add_extent = combined_extent(MDT)
                                            new_extents.append(add_extent)
                                            reply_polygons.append(MDT)
                                            reply_categories.append('MDT')
                                            order.append(2.5)
                                # if SLGT contains ENH & MDT, but ENH doesn’t contain MDT.
                                if not three and five and six:
                                    # ...add ENH extent
                                    print(f'    ENH {e+1} has no MDT')
                                    add_extent = combined_extent(ENH)
                                    new_extents.append(add_extent)
                                    reply_polygons.append(ENH)
                                    reply_categories.append('ENH')
                                    order.append(3)
                                # if SLGT contains ENH, and neither contain MDT.
                                if not three and not five and six:
                                    # ...add SLGT extent
                                    print(f'    SLGT {s+1} contains ENH {e+1}')
                                    add_extent = combined_extent(SLGT,ENH)
                                    new_extents.append(add_extent)
                                    reply_polygons.append(ENH)
                                    reply_categories.append('ENH')
                                    order.append(4)
                                    # Check if SLGT is way bigger than ENH.
                                    if top_risk=='enh':
                                        size_small = size_check(SLGT,ENH)
                                        if size_small and plot_day<3:
                                            # ...add MDT extent
                                            print(f'    ENH {e+1} much smaller than SLGT {s+1}')
                                            add_extent = combined_extent(ENH)
                                            new_extents.append(add_extent)
                                            reply_polygons.append(ENH)
                                            reply_categories.append('ENH')
                                            order.append(4.5)
                                # if MRGL contains SLGT & ENH, but SLGT doesn't contain ENH.
                                if not six and seven and eight:
                                    # ...add SLGT extent
                                    print(f'    SLGT {s+1} has no ENH')
                                    add_extent = combined_extent(SLGT)
                                    new_extents.append(add_extent)
                                    reply_polygons.append(SLGT)
                                    reply_categories.append('SLGT')
                                    order.append(5)
                                # if MRGL contains SLGT, but neither contain ENH.
                                if not six and not seven and eight and enh_exists:
                                    # ... add SLGT extent
                                    print(f'    SLGT {s+1} not in ENH')
                                    add_extent = combined_extent(SLGT)
                                    new_extents.append(add_extent)
                                    reply_polygons.append(SLGT)
                                    reply_categories.append('SLGT')
                                    order.append(6)

                                ### ZOOMED OUT ###

                                # If SLGT contains HIGH.
                                if four and high_exists:
                                    print(f'    Zoomed out: SLGT {s+1} in HIGH')
                                    add_extent = combined_extent(SLGT,ENH,MDT,HIGH)
                                    new_extents.append(add_extent)
                                    reply_polygons.append(SLGT)
                                    reply_categories.append('HIGH')
                                    order.append(7)
                                # If SLGT contains MDT.
                                elif five and mdt_exists:
                                    # ...add SLGT extent
                                    print(f'    Zoomed out: SLGT {s+1} in MDT')
                                    add_extent = combined_extent(SLGT,ENH,MDT)
                                    new_extents.append(add_extent)
                                    reply_polygons.append(SLGT)
                                    reply_categories.append('MDT')
                                    order.append(8)

                                ### ONLY SLIGHT ###

                                # if MRGL contains SLGT, and there is no ENH+.
                                if eight and not enh_exists:
                                    # ... add MRGL+SLGT extent
                                    print(f'    MRGL {r+1} contains SLGT {s+1}')
                                    add_extent = combined_extent(MRGL,SLGT)
                                    new_extents.append(add_extent)
                                    reply_polygons.append(SLGT)
                                    reply_categories.append('SLGT')
                                    order.append(9)

            # Remove extent duplicates and order from highest to lowest risk.
            extents = keep_unique_extents(new_extents,order)

            # Use the list of unique extents to find the polygons to use
            # for timestamp formatting while keeping the order specified above.
            polygons_for_timestamps=[]
            categories_for_roads=[]
            for i,item in enumerate(extents):
                idx = new_extents.index(item)
                polygons_for_timestamps.append(reply_polygons[idx])
                categories_for_roads.append(reply_categories[idx])

            extent=False
            leg_loc = 0
        else:
            # Plot CONUS with no categories.
            extent = [-125,-66,24,50]
            leg_loc = 4

    #
    # On to steps 3-5, plotting the maps!
    #
    if not extents:
        reply = False
        plot_SPC_outlook(where,plot_type,plot_type_override,plot_day,grid_res,override,extent,leg_loc,cat_gdf,reply)
    elif not extent and len(extents)>0:
        for i,extent in enumerate(extents):
            reply = False if i==0 else i
            plot_SPC_outlook(where,plot_type,plot_type_override,plot_day,grid_res,override,extent,leg_loc,cat_gdf,reply,this_polygon=polygons_for_timestamps[i],this_category=categories_for_roads[i])



###################################################
# Makes the grid, smooths it, and plots the data. #
###################################################
def plot_SPC_outlook(where,plot_type,plot_type_override,plot_day,grid_res,override,extent,leg_loc,cat_gdf,reply,this_polygon=False,this_category=False):
    #
    # STEP 3: Make the lat/lon grid.
    #

    if plot_type == 'smooth':
        start_timer = dt.now()
        print('\n\n--> Make grid')

        # Make a grid of 0.1x0.1 degrees that covers the CONUS.
        #   CONUS covers N,S,W,E: 50,24,-125,-66

        if grid_res=='high':
            y, x = grid_points = np.mgrid[24:50:2601j, -125:-66:5901j]
            US_mask_count = 8357331     # 54.45%
            grid_size = 15348501
        elif grid_res=='mid':
            y, x = grid_points = np.mgrid[24:50:521j, -125:-66:1181j]
            US_mask_count = 334094      # 54.30%
            grid_size = 615301
        elif grid_res=='low':
            y, x = grid_points = np.mgrid[24:50:261j, -125:-66:591j]
            US_mask_count = 83477       # 54.11%
            grid_size = 154251


        #
        # Get the US's shape from natural_earth's countries.
        #   Helps speed up the smoothing process.
        #

        countries = get_shapes_list('50m','cultural','admin_0_countries')
        US, = [country for country in countries if country.attributes['NAME_EN'] == 'United States of America']
        US_geom = US.geometry

        # Mask of the United States.
        US_mask = shapely.vectorized.contains(US_geom, x, y)


        time_elapsed = dt.now() - start_timer
        tsec = round(time_elapsed.total_seconds(),2)
        print(f'  --> Grid made ({tsec:.2f} seconds)')



    #
    # STEP 4: Determine which grid points are in the outlook.
    #

    ### RISK LEVELS AND COLORS ###
    # Set a variable of SPC categories and their fill colors.
    cat_plot_colors = {'General Thunderstorms Risk':'#C1E9C1',
                       'Marginal Risk': '#66A366',
                       'Slight Risk': '#FFE066',
                       'Enhanced Risk': '#FFA366',
                       'Moderate Risk': '#E06666',
                       'High Risk': 'magenta'}
    prob_plot_colors = {'15% Any Severe Risk':'#FFE066',
                        '30% Any Severe Risk':'#FFA366'}

    if plot_day>3:
        cat_plot_colors = prob_plot_colors

    category_labels = [key for key in cat_plot_colors.keys()]

    if plot_type == 'smooth':
        print('--> Mask SPC categories')

        # Create grid value for each SPC category
        # 0=N/A, 1=TSTM, ... , 6=HIGH
        categories = []

        for key in category_labels:
            geometries = cat_gdf[cat_gdf['LABEL2'] == key]

            if len(geometries) > 0:
                row = geometries.index[0]
                print(f'  --> {key}: {geometries.at[row,"fill"]}')

                # Get the polygon or multipolygon for this category.
                cat_geom = geometries.at[row,'geometry']

                # Make TRUE/FALSE grid for whether location is in the polygon.
                cat_mask = shapely.vectorized.contains(cat_geom, x, y)

                # Turn False=0 and True=category value.
                cat_mask = cat_mask.astype(int)*(row+1)

                # Gather the 1D arrays for all the categories.
                categories.append(cat_mask)

        category = []
        for i in categories:
            if len(category)==0: category = i
            else: category = np.add(category,i)

        # Output results.
        if len(category)>0:
            for i in range(7):
                if i==0:
                    print(f'\nNONE: {(category == i).sum()} = {(category == i).sum()/grid_size*100:.2f}% of Grid')
                elif i<=len(cat_gdf):
                    label2 = category_labels[i-1]
                    label = list(cat_gdf[cat_gdf['LABEL2'] == label2]['LABEL'])[0]
                    print(f'{label}: {(category >= i).sum()} = {(category >= i).sum()/US_mask_count*100:.2f}% of US')
    print()



    #
    # STEP 5: For each grid point, find the average value of all points within 25 miles.
    #
    if plot_type=='smooth':
        # Trim the data to the area you want to plot.
        x,y,category,US_mask,grids = trim_to_extent(x,y,category,US_mask,extent)

        # Smooth the categories.
        if not get_average_override:
            category = get_average_values(x,y,category,US_mask,grids)

        if where=='data':
            leg_loc = legend_location(category,US_mask)


    #
    # STEP 6: Plot the maps.
    #
    start_timer = dt.now()
    print('--> Making maps')

    # Set Coordinate Reference System from the Shapefile Data
    data_crs = ccrs.PlateCarree()

    # Get time data and convert it to useful strings
    #   Used in tweet text and in the plot title
    issue_time = cat_gdf['ISSUE'][0]
    start_time = cat_gdf['VALID'][0]
    end_time = cat_gdf['EXPIRE'][0]

    start_time_dt = dt.strptime(start_time, '%Y%m%d%H%M').replace(tzinfo=tz.gettz('UTC'))
    issue_time_dt = dt.strptime(issue_time, '%Y%m%d%H%M').replace(tzinfo=tz.gettz('UTC'))

    # Use coords from polygon in middle of image to compute times to print.
    polygon = cat_gdf.iloc[-1]['geometry'] if not this_polygon else this_polygon
    start_time,end_time,issue_time = convert_datetime_from_spc_to_local(polygon,start_time,end_time,issue_time,where)

    # Generate legend patches
    legend_patches = []
    for i,risk in enumerate(reversed(category_labels)):
        if risk == 'General Thunderstorms Risk':
            patch = mpatches.Patch(color=cat_plot_colors[risk], label='Lightning Storms')
        elif plot_day<4:
            patch = mpatches.Patch(color=cat_plot_colors[risk], label=f'{5-i}: {risk}')
        else:
            patch = mpatches.Patch(color=cat_plot_colors[risk], label=f'{risk}')
        legend_patches.append(patch)

    # Setup matplotlib figure #
    ###########################
    # Make the map coordinate reference system (crs)
    if plot_type=='CONUS': map_crs = ccrs.Orthographic(central_latitude=39.833333, central_longitude=-98.583333)
    else:
        clon = (extent[0]+extent[1])/2
        clat = (extent[2]+extent[3])/2
        map_crs = ccrs.Orthographic(central_latitude=clat, central_longitude=clon)

    # Make the figure and axis.
    fig = plt.figure(1, figsize=(1024/48, 512/48))
    ax = plt.subplot(1, 1, 1, projection=map_crs)

    print('--> Plot SPC polygons')

    # Create colormap
    cat_fill_colors = [(1,1,1,0)]
    for color in cat_plot_colors: cat_fill_colors.append(cat_plot_colors[color])

    cmap_name='SPC'
    n_bin = 100
    cm = LinearSegmentedColormap.from_list(cmap_name, cat_fill_colors, N=n_bin)

    # Plot the categories! #
    ########################
    if plot_type=='smooth':
        kwargs = {'cmap':cm,
                    'vmin':0,
                    'vmax':6,
                    'rasterized':True,
                    'snap':True,
                    'transform':data_crs}
        if plot_day >3: kwargs['vmax']=2
        pm = ax.pcolormesh(x,y,category,**kwargs)

    elif plot_type=='exact':
        for i,key in enumerate(category_labels):
            geometries = cat_gdf[cat_gdf['LABEL2'] == key]

            # Check to see if there an area outlooked at all. 
            # If so, add the polygons to the map.
            if len(geometries) > 0:
                print('  --> Plot ax.add_geometries')
                ax.add_geometries(geometries['geometry'], crs=data_crs,
                              facecolor=cat_fill_colors[i+1],alpha=1)

    # Set plot extent #
    ###################
    ax.set_extent(extent, data_crs)

    # Adjust extent
    #   Find lat of bottom middle of graph.
    p_a = (0.5,0)
    p_a_disp = ax.transAxes.transform(p_a)
    p_a_data = ax.transData.inverted().transform(p_a_disp)
    p_a_cart = data_crs.transform_point(*p_a_data, src_crs=map_crs)

    # How much should the plot be moved south by?
    south_offset = extent[2]-p_a_cart[1]
    extent = [extent[0],extent[1],extent[2]+south_offset,extent[3]]

    # Set new extent
    ax.set_extent(extent, data_crs)

    # Add map features #
    ####################
    print('  --> Adding cfeatures')

    if plot_type=='exact':
        if grid_res=='high':
            y, x = grid_points = np.mgrid[24:25:101j, -125:-124:101j]
        elif grid_res=='mid':
            y, x = grid_points = np.mgrid[24:24:21j, -125:-124:21j]
        elif grid_res=='low':
            y, x = grid_points = np.mgrid[24:25:11j, -125:-124:11j]
        grids = dist_in_grids(x,extent)

    num_grids = int(round( (((extent[3]-extent[2])*10)+(2*grids)+1) *
                            (((extent[1]-extent[0])*10)+(2*grids)+1) ))

    # Ocean, Land, Coastline
    if plot_type=='exact': ax.add_feature(cfeature.OCEAN.with_scale('50m'))
    elif plot_type=='smooth': ax.add_feature(cfeature.OCEAN.with_scale('50m'),zorder=2,edgecolor='k')
    ax.add_feature(cfeature.LAND.with_scale('50m'),facecolor='w')
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'))


    # Add major roads
    cats = ['TSTM','MRGL','SLGT','ENH','MDT','HIGH']
    road_cat = len(cat_gdf)-1 if not this_category else cats.index(this_category)
    print(f' ****** road_cat: {road_cat} ****** ')
    if plot_day==1 and road_cat>2 and plot_type=='smooth' and num_grids<=(170*170):
        if num_grids<=(100*100): rank = 7
        elif num_grids<=(140*140): rank = 6
        elif num_grids<=(160*160): rank = 5
        elif num_grids<=(170*170): rank = 4

        roads = get_shapes_list('10m','cultural','roads')
        roads = [road for road in roads
                    if road.attributes['sov_a3']=='USA'
                    #and road.attributes['type'] in ['Major Highway','Beltway']
                    and road.attributes['scalerank']<=rank]

        #print(set([road.attributes['scalerank'] for road in roads]))
        #for i in range(3,8): print(i,len([road for road in roads if road.attributes['scalerank']==i]))

        roads = [rd.geometry for rd in roads]

        for rd in roads:
            rd_buffer = rd.buffer(0.015)
            ax.add_geometries([rd_buffer],crs=data_crs,facecolor=(1,1,1),edgecolor='none')
            ax.add_geometries([rd],crs=data_crs,facecolor='none',edgecolor=(0,0,0))


    # Show county borders
    if plot_day in [1,2] and len(cat_gdf)>2 and plot_type=='smooth' and num_grids<=(180*180):
        COUNTIES = get_counties(data_crs)
        ax.add_feature(COUNTIES, facecolor='none', edgecolor='dimgray', linewidth=0.25)


    # Show urban areas
    if plot_day in [1,2] and len(cat_gdf)>2 and plot_type=='smooth':
        if num_grids<=(160*160):
            urban = cfeature.NaturalEarthFeature('cultural','urban_areas','10m')
            ax.add_feature(urban,facecolor=(0,0,0,0.2))
        elif num_grids<=(275*275):
            urban = cfeature.NaturalEarthFeature('cultural','urban_areas','50m')
            ax.add_feature(urban,facecolor=(0,0,0,0.2))


    # Add text labels for select cities
    if plot_day in [1,2] and len(cat_gdf)>2 and plot_type=='smooth' and num_grids<=(275*275):
        if num_grids<=(110*110): rank = 7
        elif num_grids<=(140*140): rank = 6
        elif num_grids<=(160*160): rank = 4
        elif num_grids<=(275*275): rank = 3

        names_records = get_shapes_list('10m','cultural','populated_places_simple')
        names = [shp for shp in names_records
                    if shp.attributes['adm0name']=='United States of America'
                    and shp.attributes['scalerank']<=rank
                    and extent[0]+0.4<float(shp.attributes['longitude'])<extent[1]-0.4
                    and extent[2]+0.4<float(shp.attributes['latitude'])<extent[3]-0.4]
        name = [pt.attributes['name'] for pt in names]
        x = [pt.attributes['longitude'] for pt in names]
        y = [pt.attributes['latitude'] for pt in names]

        for i,_ in enumerate(x):
            kwargs = {'horizontalalignment':'center',
                        'verticalalignment':'top',
                        'fontsize':8,
                        'clip_on':True,
                        'bbox':{'facecolor':'white','edgecolor':'none','alpha':0.5,'pad':1},
                        'transform':data_crs }
            ax.text(x[i],y[i],name[i],**kwargs)

    # Add state labels to day 3-8 plots.
    elif plot_day in [3,4,5,6,7,8] and plot_type=='smooth' and num_grids<=(175*175):
        # Make polygon from extent. W E S N
        corner1 = Point(extent[0]+0.4,extent[2]+0.4)
        corner2 = Point(extent[1]-0.4,extent[2]+0.4)
        corner3 = Point(extent[1]-0.4,extent[3]-0.4)
        corner4 = Point(extent[0]+0.4,extent[3]-0.4)
        corner_list = [corner1,corner2,corner3,corner4]
        extent_poly = Polygon([[p.x, p.y] for p in corner_list])

        # Get geometries and attributes of US states.
        states_records = get_shapes_list('50m','cultural','admin_1_states_provinces_lakes')
        state_info = [shp for shp in states_records
                        if shp.attributes['admin']=='United States of America'
                        and shp.attributes["name"]!='District of Columbia']

        # Find states that lie within the extent
        for state in state_info:
            if extent_poly.intersection(state.geometry):
                # Get a point within the state to place the label.
                partial_state = extent_poly.intersection(state.geometry)
                label_point = partial_state.representative_point()

                # Process state name
                state_name = state.attributes["name"].upper()
                state_name = state_name.replace(' ','\n')

                # Add state labels to the map.
                kwargs = {'horizontalalignment':'center',
                        'verticalalignment':'top',
                        'fontsize':8,
                        'clip_on':True,
                        'bbox':{'facecolor':'white','edgecolor':'none','alpha':0.5,'pad':1},
                        'transform':data_crs }
                ax.text(label_point.x,label_point.y,state_name,**kwargs)
            else: continue


    # Show states
    ax.add_feature(cfeature.STATES.with_scale('50m'))


    # Add Mexico and Canada to mask output in those areas.
    mexico_canada = Mexico_Canada()
    for country in mexico_canada:
        ax.add_geometries( [country], crs=data_crs, facecolor='w', edgecolor='k')

    # Add Great Lakes to mask output in those areas.
    great_lakes = Great_Lakes()
    for lake in great_lakes:
        ax.add_geometries( [lake], crs=data_crs, facecolor=cfeature.COLORS['water'], edgecolor='k' )


    # Legend, Title, Attribution #
    ##############################
    print('  --> Legend, Title, Attribution')
    # Plot the legend
    kwargs = {'loc':leg_loc,
                'fontsize':'medium',
                'framealpha':0.9
             }
    plt.legend(handles=legend_patches,**kwargs).set_zorder(7)

    # Plot the titles
    mid = (fig.subplotpars.right + fig.subplotpars.left)/2
    y_pos = ax.get_position().y1+0.05
    kwargs = {'fontsize':18,
                'fontweight':'bold',
                'x':mid,
                'y':y_pos
            }
    plt.suptitle(f'SPC Day {plot_day} Severe Storm Outlook', **kwargs)
    if plot_day==1:
        plt.title(f'{start_time} through tomorrow morning.', fontsize=12, loc='center')
    else:
        plt.title(f'{start_time} through the next morning.', fontsize=12, loc='center')


    # Put attribution text in opposite corner as the legend.
    if leg_loc==4:
        att_x = 0.006
        att_y = 0.01
        att_ha='left'
    else:
        att_x = 1-0.006
        att_y = 0.01
        att_ha='right'

    # Add attribution
    text = '@SmoothedPC'
    kwargs = {'weight':'bold',
                    'va':'bottom',
                    'ha':att_ha,
                    'snap':True,
                    'zorder':7,
                    'transform':ax.transAxes
                }
    ax.text(att_x,att_y,text,**kwargs)


    time_elapsed = dt.now() - start_timer
    tsec = round(time_elapsed.total_seconds(),2)
    print(f'--> Map made ({tsec:.2f} seconds)')


    # Save #
    ########
    print('--> Save')
    start_timer = dt.now()

    save_location = f'images/day{plot_day}_categorical.png'

    if plot_type=='exact':
        save_location = f'images/day{plot_day}_categorical.png'

        # Save the figure.
        plt.savefig(save_location, dpi=96, bbox_inches='tight')

        # Copy file to places.
        shutil.copy2(save_location,f'images/latest_day{plot_day}_categorical.png')
        shutil.copy2(save_location,f'images/latest_exact.png')

    elif plot_type=='smooth':
        if not reply: save_location = f'images/day{plot_day}_grid_categorical.png'
        else: save_location = f'images/day{plot_day}_{reply}_grid_categorical.png'

        # Save the figure.
        plt.savefig(save_location, dpi=96, bbox_inches='tight')

        # Copy file to places.
        if not reply:
            shutil.copy2(save_location,f'images/latest_day{plot_day}_categorical.png')
            shutil.copy2(save_location,f'images/latest_smooth.png')

        # If high risk, keep it!
        if len(cat_gdf)==6:
            shutil.copy2(save_location,f'images/latest_high.png')

        # Tweet the result.
        print('  --> Smooth and tweet.')
        US_time_dt = start_time_dt.astimezone(tz.gettz('America/New_York'))

        if plot_day==1:
            if 3<int(issue_time_dt.strftime('%-H'))<10:
                tweet_text = f'Issued at {issue_time}: SPC forecast for the upcoming day, {US_time_dt:%A, %B %-d}.'
            elif int(issue_time_dt.strftime('%-H'))<2:
                tweet_text = f'Issued at {issue_time}: SPC forecast for TONIGHT, {US_time_dt:%A, %B %-d}.'
            else:
                tweet_text = f'Issued at {issue_time}: SPC forecast for TODAY, {US_time_dt:%A, %B %-d}.'
            print(f'    --> Tweet: {tweet_text}')
            tweet(tweet_text, save_location, send_tweet, reply)

        elif plot_day==2:
            if int(issue_time_dt.strftime('%-H'))<12:
                tweet_text = f'Issued at {issue_time}: SPC forecast for {start_time_dt:%A, %B %-d}.'
            else:
                tweet_text = f'Issued at {issue_time}: SPC forecast for tomorrow, {start_time_dt:%A, %B %-d}.'
            print(f'    --> Tweet: {tweet_text}')
            tweet(tweet_text, save_location, send_tweet, reply)

        else:
            tweet_text = f'SPC forecast for {start_time_dt:%A, %B %-d}.'
            print(f'    --> Tweet: {tweet_text}')
            tweet(tweet_text, save_location, send_tweet, reply)

        print('  --> Smooth: Tweeted.')



    # Clear figure.
    plt.clf()


    time_elapsed = dt.now() - start_timer
    tsec = round(time_elapsed.total_seconds(),2)
    print(f'--> Saved ({tsec:.2f} seconds)')








if __name__ == '__main__':
    time = dt.utcnow()
    if time.hour in [1,6,12,13,16,20]: h1=1; h2=2; st=1
    elif time.hour in [17]: h1=2; h2=3; st=1
    elif time.hour in [7]: h1=2; h2=4; st=1
    else: h1=8; h2=3; st=-1

    if override:
        if isinstance(plot_day,int):
            h1=plot_day; h2=plot_day+1; st=1

    # Time it
    big_start_timer = dt.now()

    # Run the code
    if isinstance(plot_day,int):
        for plot_day in range(h1,h2,st):
            print(f'\n*** Day {plot_day} ***')
            get_SPC_data(where,plot_type,plot_type_override,plot_day,grid_res,override)
    elif isinstance(plot_day,list):
        print(f'\n*** Date: {plot_day[0]}_{plot_day[1]} ***')
        get_SPC_data(where,plot_type,plot_type_override,plot_day,grid_res,override)

    # Report time
    big_time_elapsed = dt.now() - big_start_timer
    tsec = round(big_time_elapsed.total_seconds(),2)
    print(f'\n--> Done ({tsec:.2f} seconds)')
