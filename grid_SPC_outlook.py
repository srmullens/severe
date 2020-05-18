################################################################
# Plots SPC Outlooks for the area you are concerned about.     #
#   URLs exist for plotting WPC, NHC, drought data.            #
#                                                              #
# Code by: Stephen Mullens                                     #
# May 2020                                                     #
#                                                              #
# Code inspired by https://github.com/frontogenesis/metpy/     #
#   look at .ipynb_checkpoints/SPC Outlooks-checkpoint.ipynb   #
#   and at .ipynb_checkpoints/download_shapes-checkpoint.ipynb #
################################################################

import requests
from contextlib import closing
import zipfile
import os
import glob
import shutil
import numpy as np
import pandas as pd
from math import degrees, radians, cos, sin, asin, sqrt

import geopandas
import shapely.vectorized

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
setting = 'low'

# What type of plot: 'exact' or 'smooth'
plot_type_override = False
plot_type = 'smooth'

# What SPC day do you want to plot?
plot_day = 2

# Need a time, plot_type, and plot_day override?
override = False
if override: plot_type_override=override


#
# Make some functions.
#

# Downloads a zip file and extracts it in a target directory of choice
def download_zip_files(issuing_center,product):
    for issuing_center in shapefiles:
        for product in shapefiles[issuing_center]:
            download_zip_file(file_url = shapefiles[issuing_center][product], root_folder = issuing_center)

    for zip_file in glob.glob(f'{os.getcwd()}/*.zip'):
        head, tail = os.path.split(zip_file)
        shutil.move(zip_file, f'zips/{tail}')


def download_zip_file(file_url, root_folder):
    if not isinstance(file_url,str):
        raise TypeError(f'file URL must by type string, not {type(file_url)}')
    if not isinstance(root_folder,str):
        raise TypeError(f'folder must by type string, not {type(root_folder)}')

    file = file_url.split('/')[-1]
    folder = file.split('.')[0]

    with requests.get(file_url, stream=True) as r:
        with open(file, 'wb') as f:
            for chunk in r.iter_content(chunk_size=128):
                f.write(chunk)


    with zipfile.ZipFile(file, 'r') as zip_ref:
        zip_ref.extractall(f'{root_folder}/{folder}')
        zip_ref.close()


# Clears all contents of folder specified in the argument
def clear_folder_contents(folder):

    folder = f'{folder}/'
    extensions = ['.zip','.dbf','.prj','.shp','.shx']
    last = ['-shp']

    for file in os.listdir(folder):
        file_path = os.path.join(folder, file)
        try:
            if os.path.isfile(file_path) and file[-4] in extensions:
                os.unlink(file_path)
            elif os.path.isdir(file_path) and file[-4] in last:
                shutil.rmtree(file_path)
        except Exception as e:
            print(e)


# Calculate how many grids 25km is
def dist_in_grids(x,y,data,extent):

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



# Trim datasets to right around the extent of the map for easier plotting.
def trim_to_extent(x,y,data,mask,extent):
    start_timer = dt.now()
    print('--> Trimming data')

    # How many grids do you include for averaging?
    grids = dist_in_grids(x,y,data,extent)

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
                            #if lat==int(i/2) and lon==int(i/2):
                            #    print(f'  {orig_x}/{orig_y}->{include_x}/{include_y}->{dist:.1f}: {include_value}')

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


# Convert time zones for printing
def convert_datetime_from_spc_to_local(polygon,string,start_end,from_zone,to_zone='America/New_York'):
    utc_time = dt.strptime(string, '%Y%m%d%H%M').replace(tzinfo=from_zone)

    # Set up getting the time zone.
    tf = TimezoneFinder()
    new_zones_list = []

    # Extract a list of coordinates from the polygon(s).
    try:
        if polygon.geom_type == 'MultiPolygon':
            multipolygon = polygon
            coords = [point for polygon in multipolygon for point in polygon.exterior.coords[:-1]]
        elif polygon.geom_type == 'Polygon':
            coords = list(polygon.exterior.coords)
    except:
        coords = None

    if coords is not None:
        for coord in coords:
            new_zones_list.append(tf.timezone_at(lng=coord[0], lat=coord[1]))

        # Remove None values.
        new_zones_list = [i for i in new_zones_list if i]

        # Standardize time zones
        for i,item in enumerate(new_zones_list):
            zone = f'{utc_time.astimezone(tz.gettz(item)):%Z}'
            if zone in ['EST','EDT']:
                new_zones_list[i] = 'America/New_York'
            elif zone in ['CST','CDT']:
                new_zones_list[i] = 'America/Chicago'
            elif zone in ['MST','MDT']:
                new_zones_list[i] = 'America/Denver'
            elif zone in ['PST','PDT']:
                new_zones_list[i] = 'America/Los_Angeles'

        # Sort the resulting list from west to east.
        new_zones_list = list(set(new_zones_list))

        # Sort the resulting list from west to east.
        west_to_east = ['America/Los_Angeles','America/Denver','America/Chicago','America/New_York']
        sort_by = []
        print(new_zones_list)
        for item in new_zones_list: sort_by.append(west_to_east.index(item))
        new_zones_list = [nzl for _,nzl in sorted(zip(sort_by,new_zones_list), key=lambda pair: pair[0])]

        print(new_zones_list)

        # Use str time zone names to modify datetime objects
        for i,item in enumerate(new_zones_list):
            new_zones_list[i] = utc_time.astimezone(tz.gettz(item))

    else:
        new_zones_list=[utc_time]

    # Generate string outputs
    tweet_valid_time = ''
    if len(new_zones_list)>=3:
        if start_end=='start':
            date_time = f'{new_zones_list[0]:%a, %b %d, %Y %-I:%M}{new_zones_list[0].strftime("%p").lower()} {new_zones_list[0]:%Z}, {new_zones_list[-1]:%-I:%M}{new_zones_list[-1].strftime("%p").lower()} {new_zones_list[-1]:%Z}'
            tweet_valid_time = f'{new_zones_list[0]:%-I:%M}{new_zones_list[0].strftime("%p").lower()} {new_zones_list[0]:%Z}, {new_zones_list[-1]:%-I:%M}{new_zones_list[-1].strftime("%p").lower()} {new_zones_list[-1]:%Z}'
        elif start_end=='end':
            date_time = f'{new_zones_list[0]:%-I:%M}{new_zones_list[0].strftime("%p").lower()} {new_zones_list[0]:%Z}, {new_zones_list[-1]:%-I:%M}{new_zones_list[-1].strftime("%p").lower()} {new_zones_list[-1]:%Z}'

    elif len(new_zones_list)==2:
        if start_end=='start':
            date_time = f'{new_zones_list[0]:%a, %b %d, %Y %-I:%M}{new_zones_list[0].strftime("%p").lower()} {new_zones_list[0]:%Z}, {new_zones_list[1]:%-I:%M}{new_zones_list[1].strftime("%p").lower()} {new_zones_list[1]:%Z}'
            tweet_valid_time = f'{new_zones_list[0]:%-I:%M}{new_zones_list[0].strftime("%p").lower()} {new_zones_list[0]:%Z}, {new_zones_list[1]:%-I:%M}{new_zones_list[1].strftime("%p").lower()} {new_zones_list[1]:%Z}'
        elif start_end=='end':
            date_time = f'{new_zones_list[0]:%-I:%M}{new_zones_list[0].strftime("%p").lower()} {new_zones_list[0]:%Z}, {new_zones_list[1]:%-I:%M}{new_zones_list[1].strftime("%p").lower()} {new_zones_list[1]:%Z}'

    elif len(new_zones_list)==1:
        if start_end=='start':
            date_time = f'{new_zones_list[0]:%a, %b %d, %Y %-I:%M}{new_zones_list[0].strftime("%p").lower()}'
            tweet_valid_time = f'{new_zones_list[0]:%-I:%M}{new_zones_list[0].strftime("%p").lower()} {new_zones_list[0]:%Z}'
        elif start_end=='end':
            date_time = f'{new_zones_list[0]:%-I:%M}{new_zones_list[0].strftime("%p").lower()} {new_zones_list[0]:%Z}'

    print(f'  --> {start_end}: {date_time}')
    print(f'  --> {tweet_valid_time}')

    return date_time, tweet_valid_time


# Get the country outlines of Mexico and Canada for the map.
def Mexico_Canada():
    shpfilename = shpreader.natural_earth(resolution='50m',
                                    category='cultural',
                                    name='admin_0_countries')
    reader = shpreader.Reader(shpfilename)
    countries = list(reader.records())

    country_list = ['Mexico','Canada']

    MC = [country for country in countries if country.attributes['NAME_EN'] in country_list]
    Mex_Can = [mc.geometry for mc in MC]

    return Mex_Can


# Get the outlines of the Great Lakes for the map.
def Great_Lakes():
    shpfilename = shpreader.natural_earth(resolution='50m',
                                      category='physical',
                                      name='lakes')
    reader = shpreader.Reader(shpfilename)
    all_lakes = list(reader.records())

    list_of_lakes = ['Lake Superior','Lake Huron','Lake Michigan','Lake Erie','Lake Ontario']

    GL = [lake for lake in all_lakes if lake.attributes['name'] in list_of_lakes]
    lakes = [gl.geometry for gl in GL]

    return lakes


# Determine an ideal legend location.
def legend_location(category):
    print('--> Determining Legend Location')

    # Get items for slicing the array
    height,width = category.shape
    dh = int(round(height*0.25,0))
    dw = int(round(width*0.25,0))

    # Max category within the slice.
    upper_right_max = np.amax(category[height-dh:height,width-dw:width])
    upper_left_max = np.amax(category[height-dh:height,0:dw])
    lower_right_max = np.amax(category[0:dh,width-dw:width])
    lower_left_max = np.amax(category[0:dh,0:dw])
    max_corners = [lower_right_max,lower_left_max,upper_right_max,upper_left_max]

    # Sum of category values across the slice.
    upper_right_sum = np.sum(category[height-dh:height,width-dw:width])
    upper_left_sum = np.sum(category[height-dh:height,0:dw])
    lower_right_sum = np.sum(category[0:dh,width-dw:width])
    lower_left_sum = np.sum(category[0:dh,0:dw])
    sum_corners = [lower_right_sum,lower_left_sum,upper_right_sum,upper_left_sum]

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


# Tweet the results.
def tweet(text, image):
    consumer_key = os.environ.get('consumer_key')
    consumer_secret = os.environ.get('consumer_secret')
    access_token = os.environ.get('access_token')
    access_token_secret = os.environ.get('access_token_secret')

    print('--> Tweeting...')
    twitter = Twython(consumer_key, consumer_secret, access_token, access_token_secret)
    response = twitter.upload_media(media=open(image, 'rb'))
    twitter.update_status(status=text, media_ids=[response['media_id']])
    print('  --> Tweeted.')




# The main function to make the plots.
def grid_SPC_outlook(where,plot_type,plot_type_override,plot_day,setting,override):
    #
    # STEP 1: Get the files
    #

    big_start_timer = dt.now()
    print('--> Get the files')

    if plot_day<4:
        shapefiles = {
            'spc': {f'categorical_day{plot_day}': f'https://www.spc.noaa.gov/products/outlook/day{plot_day}otlk-shp.zip'}
                    }
    else:
        shapefiles = {
            'spc': {f'categorical_day{plot_day}': f'https://www.spc.noaa.gov/products/exper/day4-8/day{plot_day}prob-shp.zip'}
                    }
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
    target_folders = ['spc', 'zips']
    #target_folders = ['nhc', 'spc', 'wpc', 'drought', 'zips']
    for folder in target_folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
        else:
            clear_folder_contents(folder)

    # Get the shapefiles
    print('--> Getting ZIPs')
    start_timer = dt.now()

    """
    for issuing_center in shapefiles:
        for product in shapefiles[issuing_center]:
            download_zip_file(file_url = shapefiles[issuing_center][product], root_folder = issuing_center)

    for zip_file in glob.glob(f'{os.getcwd()}/*.zip'):
        head, tail = os.path.split(zip_file)
        shutil.move(zip_file, f'zips/{tail}')
    """

    # Read in Shapefile
    tries = 1
    if plot_day<3:
        while tries<30:
            # Download ZIP files
            download_zip_files(issuing_center,product)
            # Read file you want
            cat_gdf = geopandas.read_file(f'spc/day{plot_day}otlk-shp/day{plot_day}otlk_cat.shp')
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
    else:
        while tries<30:
            try:
                # Download ZIP files
                download_zip_files(issuing_center,product)
                # Read file you want
                cat_gdf = geopandas.read_file(f'spc/day{plot_day}prob-shp/day{plot_day}otlk_{big_start_timer:%Y%m%d}_prob.shp')
                print(f"  --> Got it! {dt.utcnow():%H%M} UTC - {dt.strptime(cat_gdf['ISSUE'][0],'%Y%m%d%H%M').strftime('%H%M')}")
                tries+=30
            except:
                print(f'  --> Not available yet. {dt.utcnow():%H%M} UTC')
                if tries==29: print(f'Could not find day{plot_day}otlk_{big_start_timer:%Y%m%d}_prob.shp after 15 minutes.'); return
                else:
                    t.sleep(120)
                    tries+=1

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
    # STEP 2: Set some variables we'll need for later.
    #

    # Set a variable of SPC categories and their fill colors.
    cat_plot_colors = {'General Thunderstorms Risk':'#C1E9C1',
                       'Marginal Risk': '#66A366',
                       'Slight Risk': '#FFE066',
                       'Enhanced Risk': '#FFA366',
                       'Moderate Risk': '#E06666',
                       'High Risk': '#e29ee9'}
    prob_plot_colors = {'15% Any Severe Risk':'#FFE066',
                        '30% Any Severe Risk':'#FFA366'}

    cat_fill_colors = [(1,1,1,0),'#C1E9C1','#66A366','#FFE066','#FFA366','#E06666','magenta']
    prob_fill_colors = [(1,1,1,0),'#FFE066','#FFA366']

    if plot_day>3:
        cat_plot_colors = prob_plot_colors
        cat_fill_colors = prob_fill_colors

    """
    # Set colors
    cat_plot_colors = {'Marginal Risk': 'green',
                       'Slight Risk': 'yellow',
                       'Enhanced Risk': 'orange',
                       'Moderate Risk': 'red',
                       'High Risk': 'magenta'}
    """


    # Set Coordinate Reference System, extent, and legend location for the map
    if plot_nothing or where=='CONUS':
        map_crs = ccrs.Orthographic(central_latitude=39.833333, central_longitude=-98.583333)
        extent = [-125,-66,24,50]
        leg_loc = 4
    elif where == 'Southeast':
        map_crs = ccrs.Orthographic(central_latitude=30.7, central_longitude=-88)
        extent = [-95, -75, 24, 35]
        leg_loc = 4
    elif where=='Florida':
        map_crs = ccrs.Orthographic(central_latitude=28, central_longitude=-83)
        extent = [-88,-79,24,32]
        leg_loc = 3
    elif where=='data':
        extents = []
        # Find the extent of each category in the outlook.
        for key in cat_plot_colors.keys():
            geometries = cat_gdf[cat_gdf['LABEL2'] == key]
            if len(geometries) > 0:
                geo_extent = geometries.bounds
                extents.append([float(geo_extent['minx'])-1,
                                float(geo_extent['maxx'])+1,
                                float(geo_extent['miny'])-1,
                                float(geo_extent['maxy'])+1 ])
        if len(extents)>1:
            # Plot extent of second highest category [-2].
            extent = extents[-2]
            clon = (extent[0]+extent[1])/2
            clat = (extent[2]+extent[3])/2
            map_crs = ccrs.Orthographic(central_latitude=clat, central_longitude=clon)
            leg_loc = 0
        elif len(extents)==1:
            # Plot extent of only category.
            extent = extents[0]
            clon = (extent[0]+extent[1])/2
            clat = (extent[2]+extent[3])/2
            map_crs = ccrs.Orthographic(central_latitude=clat, central_longitude=clon)
            leg_loc = 0
        else:
            # Plot CONUS with no categories.
            extent = [-125,-66,24,50]
            map_crs = ccrs.Orthographic(central_latitude=39.833333, central_longitude=-98.583333)
            leg_loc = 4

    # Put attribution text in opposite corner as the legend.
    if leg_loc==3:
        att_x = 1-0.006; att_y = 0.01; att_ha='right'
    if leg_loc==4:
        att_x = 0.006; att_y = 0.01; att_ha='left'



    #
    # STEP 3: Make the lat/lon grid.
    #

    if plot_type == 'smooth':
        start_timer = dt.now()
        print('--> Make grid')

        # Make a grid of 0.01x0.01 degrees that covers the CONUS.
        #   CONUS covers N,S,W,E: 50,24,-125,-66

        if setting=='high':
            y, x = grid_points = np.mgrid[24:50:2601j, -125:-66:5901j]
            US_mask_count = 8357331     # 54.45%
            grid_size = 15348501
        elif setting=='mid':
            y, x = grid_points = np.mgrid[24:50:521j, -125:-66:1181j]
            US_mask_count = 334094      # 54.30%
            grid_size = 615301
        elif setting=='low':
            y, x = grid_points = np.mgrid[24:50:261j, -125:-66:591j]
            US_mask_count = 83477       # 54.11%
            grid_size = 154251


        #
        # Get the US's shape from natural_earth's countries.
        #   Helps speed up the smoothing process.
        #

        shpfilename = shpreader.natural_earth(resolution='10m',
                                              category='cultural',
                                              name='admin_0_countries')
        reader = shpreader.Reader(shpfilename)
        countries = list(reader.records())

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
    if plot_type == 'smooth':
        print('--> Mask SPC categories')

        # Create grid value for each SPC category
        # 0=N/A, 1=TSTM, ... , 6=HIGH
        categories = []

        for key in cat_plot_colors.keys():
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
                    label2 = [key for key in cat_plot_colors.keys()][i-1]
                    label = list(cat_gdf[cat_gdf['LABEL2'] == label2]['LABEL'])[0]
                    print(f'{label}: {(category >= i).sum()} = {(category >= i).sum()/US_mask_count*100:.2f}% of US')
            print('\n')



    #
    # STEP 5: For each grid point, find the average value of all points within 25 miles.
    #
    if plot_type=='smooth':
        # Trim the data to the area you want to plot.
        x,y,category,US_mask,grids = trim_to_extent(x,y,category,US_mask,extent)

        # Smooth the categories.
        category = get_average_values(x,y,category,US_mask,grids)

        if where=='data':
            leg_loc = legend_location(category)

        # Mask the smoothed data by what's on land.
        #category = np.ma.masked_where(US_mask==False, category)


    #
    # STEP 5: Plot the maps.
    #
    start_timer = dt.now()
    print('--> Making maps')

    # Set Coordinate Reference System from the Shapefile Data
    data_crs = ccrs.PlateCarree()

    # Get time data
    start_time = cat_gdf['VALID'][0]
    end_time = cat_gdf['EXPIRE'][0]

    start_time_dt = dt.strptime(start_time, '%Y%m%d%H%M')

    from_zone = tz.gettz('UTC')
    to_zone = tz.gettz('America/New_York')

    polygon = cat_gdf.iloc[-1]['geometry']

    start_time,tweet_valid_time = convert_datetime_from_spc_to_local(polygon,start_time,'start',from_zone,to_zone)
    end_time,dummy = convert_datetime_from_spc_to_local(polygon,end_time,'end',from_zone,to_zone)

    # Generate legend patches
    legend_patches = []
    for i,risk in enumerate(reversed(list(cat_plot_colors.keys()))):
        if risk == 'General Thunderstorms Risk':
            patch = mpatches.Patch(color=cat_plot_colors[risk], label='Lightning Storms')
        elif plot_day<4:
            patch = mpatches.Patch(color=cat_plot_colors[risk], label=f'{5-i}: {risk}')
        else:
            patch = mpatches.Patch(color=cat_plot_colors[risk], label=f'{risk}')
        legend_patches.append(patch)

    # Setup matplotlib figure
    fig = plt.figure(1, figsize=(1024/48, 512/48))
    ax = plt.subplot(1, 1, 1, projection=map_crs)

    print('--> Plot SPC polygons')

    # Create colormap
    cmap_name='SPC'
    n_bin = 100
    cm = LinearSegmentedColormap.from_list(cmap_name, cat_fill_colors, N=n_bin)

    ########################
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

        for i,key in enumerate(cat_plot_colors.keys()):
            geometries = cat_gdf[cat_gdf['LABEL2'] == key]

            # Check to see if there an area outlooked at all. If so, add the polygons to the map.
            if len(geometries) > 0:
                ax.add_geometries(geometries['geometry'], crs=data_crs,
                              facecolor=cat_fill_colors[i+1],alpha=1)

    # Set plot extent
    ax.set_extent(extent, data_crs)

    # Add map features
    print('--> Adding cfeatures')
    if plot_type=='exact': ax.add_feature(cfeature.OCEAN.with_scale('50m'))
    elif plot_type=='smooth': ax.add_feature(cfeature.OCEAN.with_scale('50m'),zorder=2,edgecolor='k')
    ax.add_feature(cfeature.LAND.with_scale('50m'),facecolor='w')
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
    ax.add_feature(cfeature.STATES.with_scale('50m'))

    # Add Mexico and Canada to mask output in those areas.
    mexico_canada = Mexico_Canada()
    for country in mexico_canada:
        ax.add_geometries( [country], crs=data_crs, facecolor='w', edgecolor='k')

    # Add Great Lakes
    great_lakes = Great_Lakes()
    for lake in great_lakes:
        ax.add_geometries( [lake], crs=data_crs, facecolor=cfeature.COLORS['water'], edgecolor='k' )

    print('--> Legend, Title')
    # Plot the legend
    kwargs = {'loc':leg_loc,
                'fontsize':'medium',
                'framealpha':0.9
             }
    plt.legend(handles=legend_patches,**kwargs).set_zorder(7)

    # Plot the titles
    mid = (fig.subplotpars.right + fig.subplotpars.left)/2
    kwargs = {'fontsize':18,
                'fontweight':'bold',
                'x':mid,
                'y':0.93
            }
    plt.suptitle(f'SPC Day {plot_day} Severe Storm Outlook', **kwargs)
    if plot_day==1:
        plt.title(f'{start_time} through tomorrow morning.', fontsize=12, loc='center')
    else:
        plt.title(f'{start_time} through the next morning.', fontsize=12, loc='center')


    # Put attribution text in opposite corner as the legend.
    if leg_loc==4:
        att_x = 0.006; att_y = 0.01; att_ha='left'
    else: att_x = 1-0.006; att_y = 0.01; att_ha='right'

    # Add attribution
    text = '@SmoothedPC'
    kwargs = {'weight':'bold',
                    #'bbox':dict(boxstyle="round",ec='white',fc="white",alpha=0.75),
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


    # Save
    print('--> Save')
    start_timer = dt.now()

    if plot_type=='exact':
        # Save the figure.
        plt.savefig(f'spc/day{plot_day}_categorical.png', dpi=96, bbox_inches='tight')

        # Copy file to places.
        shutil.copy2(f'spc/day{plot_day}_categorical.png',f'spc/latest_day{plot_day}_categorical.png')
        shutil.copy2(f'spc/day{plot_day}_categorical.png',f'spc/latest_exact.png')

    elif plot_type=='smooth':
        # Save the figure.
        plt.savefig(f'spc/day{plot_day}_grid_categorical.png', dpi=96, bbox_inches='tight')

        # Copy file to places.
        shutil.copy2(f'spc/day{plot_day}_grid_categorical.png',f'spc/latest_day{plot_day}_categorical.png')
        shutil.copy2(f'spc/day{plot_day}_grid_categorical.png',f'spc/latest_smooth.png')

        # If high risk, keep it!
        if len(cat_gdf)==6:
            shutil.copy2(f'spc/day{plot_day}_grid_categorical.png',f'spc/latest_high.png')

        # Tweet the result.
        print('  --> Smooth. Tweet.')
        US_time_dt = start_time_dt.astimezone(tz.gettz('America/New_York'))

        if plot_day==1:
            if start_time_dt.strftime('%-H')=='1':
                print(f'    --> Tweet: {tweet_valid_time}: SPC forecast for TONIGHT, {US_time_dt:%A, %B %-d}.')
                tweet(f'{tweet_valid_time}: SPC forecast for TONIGHT, {US_time_dt:%A, %B %-d}.', f'spc/day{plot_day}_grid_categorical.png')
            else:
                print(f'    --> Tweet: {tweet_valid_time}: SPC forecast for TODAY, {US_time_dt:%A, %B %-d}.')
                tweet(f'{tweet_valid_time}: SPC forecast for TODAY, {US_time_dt:%A, %B %-d}.', f'spc/day{plot_day}_grid_categorical.png')
        elif plot_day==2:
            print(f'    --> Tweet: {tweet_valid_time}: SPC forecast for tomorrow, {start_time_dt:%A, %b %-d}.')
            tweet(f'{tweet_valid_time}: SPC forecast for tomorrow, {start_time_dt:%A, %b %-d}.', f'spc/day{plot_day}_grid_categorical.png')
        else:
            print(f'    --> Tweet: SPC forecast for {start_time_dt:%A, %b %-d}.')
            tweet(f'SPC forecast for {start_time_dt:%A, %b %-d}.', f'spc/day{plot_day}_grid_categorical.png')

        print('  --> Smooth. Tweeted.')

    # Clear figure.
    plt.clf()


    time_elapsed = dt.now() - start_timer
    tsec = round(time_elapsed.total_seconds(),2)
    print(f'--> Saved ({tsec:.2f} seconds)')


    big_time_elapsed = dt.now() - big_start_timer
    tsec = round(big_time_elapsed.total_seconds(),2)
    print(f'\n--> Done ({tsec:.2f} seconds)')






if __name__ == '__main__':
    time = dt.utcnow()
    if time.hour in [1,6,12,13,16,20]: h1=1; h2=2; st=1
    elif time.hour in [17]: h1=2; h2=3; st=1
    elif time.hour in [7]: h1=2; h2=4; st=1
    else: h1=8; h2=3; st=-1

    if override: h1=plot_day; h2=plot_day+1; st=1

    for plot_day in range(h1,h2,st):
        print(f'\n*** Day {plot_day} ***')
        grid_SPC_outlook(where,plot_type,plot_type_override,plot_day,setting,override)

