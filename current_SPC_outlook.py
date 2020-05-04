################################################################
# Plots SPC Outlooks for the area you are concerned about.     #
#   URLs exist for plotting WPC, NHC, drought data.            #
#                                                              #
# Code by: Stephen Mullens                                     #
# May 2020                                                     #
#                                                              #
# Original code from https://github.com/frontogenesis/metpy/   #
#   look at .ipynb_checkpoints/SPC Outlooks-checkpoint.ipynb   #
#   and at .ipynb_checkpoints/download_shapes-checkpoint.ipynb #
################################################################

import shutil
import urllib.request as request
from contextlib import closing
import zipfile
import os
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import geopandas
from cartopy import crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
from dateutil import tz


# What area do you want to plot?
where='data'
#where='CONUS'
#where = 'Southeast'
#where='Florida'

#
# STEP 1: Get the files
#

# Downloads a zip file and extracts it in a target directory of choice
def download_zip_file(file_url: str, root_folder: str):

    file = file_url.split('/')[-1]
    folder = file.split('.')[0]

    with closing(request.urlopen(file_url)) as r:
        with open(file, 'wb') as f:
            shutil.copyfileobj(r, f)

    with zipfile.ZipFile(file, "r") as zip_ref:
        zip_ref.extractall(f"{root_folder}/{folder}")
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


shapefiles = {
    "spc": {
        "categorical_day1": "https://www.spc.noaa.gov/products/outlook/day1otlk-shp.zip",
        "categorical_day2": "https://www.spc.noaa.gov/products/outlook/day2otlk-shp.zip",
        "categorical_day3": "https://www.spc.noaa.gov/products/outlook/day3otlk-shp.zip",
        "categorical_day4": "https://www.spc.noaa.gov/products/exper/day4-8/day4prob-shp.zip",
        "categorical_day5": "https://www.spc.noaa.gov/products/exper/day4-8/day5prob-shp.zip",
        "categorical_day6": "https://www.spc.noaa.gov/products/exper/day4-8/day6prob-shp.zip",
        "categorical_day7": "https://www.spc.noaa.gov/products/exper/day4-8/day7prob-shp.zip",
        "categorical_day8": "https://www.spc.noaa.gov/products/exper/day4-8/day8prob-shp.zip"
    }} #,
"""
    "wpc": {
        "excessive_rain_day1": "https://ftp.wpc.ncep.noaa.gov/shapefiles/qpf/excessive/EXCESSIVERAIN_Day1_latest.zip",
        "excessive_rain_day2" : "https://ftp.wpc.ncep.noaa.gov/shapefiles/qpf/excessive/EXCESSIVERAIN_Day2_latest.zip",
        "excessive_rain_day3" : "https://ftp.wpc.ncep.noaa.gov/shapefiles/qpf/excessive/EXCESSIVERAIN_Day3_latest.zip"
    },
    "nhc": {
        "tropical_wx_outlook": "https://www.nhc.noaa.gov/xgtwo/gtwo_shapefiles.zip"
    },
    "drought": {
        "latest": "https://droughtmonitor.unl.edu/data/shapefiles_m/USDM_current_M.zip",
    }
}
"""

# Clear directories before getting data
target_folders = ['nhc', 'spc', 'wpc', 'drought', 'zips']
for folder in target_folders:
    if not os.path.exists(folder):
        os.makedirs(folder)
    else:
        clear_folder_contents(folder)

# Get the shapefiles
print("--> Getting ZIPs")
for issuing_center in shapefiles:
    for product in shapefiles[issuing_center]:
        download_zip_file(file_url = shapefiles[issuing_center][product], root_folder = issuing_center)

for zip_file in glob.glob(f"{os.getcwd()}/*.zip"):
    head, tail = os.path.split(zip_file)
    shutil.move(zip_file, f"zips/{tail}")



#
# STEP 2: Plot the maps.
#

print("--> Making maps")

# Read in Shapefile
cat_gdf = geopandas.read_file('spc/day1otlk-shp/day1otlk_cat.shp')

# Set colors
# For SPC outlooks: "2" references a "General Thunderstorm" risk, "3" a Marginal Risk, "4" a Slight Risk, and so on

# cat_plot_colors = { #2: 'palegreen',
#                    3: 'green',
#                    4: 'yellow',
#                    5: 'brown',
#                    6: 'red',
#                    7: 'magenta'}

# Set colors
cat_plot_colors = {'Marginal Risk': 'green',
                   'Slight Risk': 'yellow',
                   'Enhanced Risk': 'orange',
                   'Moderate Risk': 'red',
                   'High Risk': 'magenta'}


# Set Coordinate Reference System for the map
if where=='CONUS':
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
        print(extent)
        clon = (extent[0]+extent[1])/2
        clat = (extent[2]+extent[3])/2
        print(clat,clon)
        map_crs = ccrs.Orthographic(central_latitude=clat, central_longitude=clon)
        leg_loc = 0
    elif len(extents)==1:
        # Plot extent of only category.
        extent = extents[0]
        clat = (extent[0]+extent[1])/2
        clon = (extent[2]+extent[3])/2
        map_crs = ccrs.Orthographic(central_latitude=clat, central_longitude=clon)
        leg_loc = 0
    else:
        # Plot CONUS with no categories.
        extent = [-125,-66,24,50]
        map_crs = ccrs.Orthographic(central_latitude=39.833333, central_longitude=-98.583333)
        leg_loc = 0

# Set Coordinate Reference System from the Shapefile Data
data_crs = ccrs.PlateCarree()

# Get time data
start_time = cat_gdf['VALID'][0]
end_time = cat_gdf['EXPIRE'][0]

from_zone = tz.gettz('UTC')
to_zone = tz.gettz('America/New_York')

def convert_datetime_from_spc_to_local(string,start_end):
    utc_time = datetime.strptime(string, '%Y%m%d%H%M').replace(tzinfo=from_zone)
    eastern = utc_time.astimezone(to_zone)
    if start_end=='start':
        date_time = datetime.strftime(eastern, '%a, %b %d, %Y %I:%M %p').lstrip('0').replace(' 0', ' ')
    elif start_end=='end':
        date_time = datetime.strftime(eastern, '%I:%M %p').lstrip('0').replace(' 0', ' ')
    return date_time

start_time = convert_datetime_from_spc_to_local(start_time,'start')
end_time = convert_datetime_from_spc_to_local(end_time,'end')

# Generate legend patches
legend_patches = []
for risk in cat_plot_colors.keys():
    patch = mpatches.Patch(color=cat_plot_colors[risk], label=risk)
    legend_patches.append(patch)

# Setup matplotlib figure
fig = plt.figure(1, figsize=(1024/96, 512/96))
ax = plt.subplot(1, 1, 1, projection=map_crs)

#Extents
ax.set_extent(extent, data_crs)

print("--> Adding cfeatures")
ax.add_feature(cfeature.OCEAN.with_scale('50m'))
ax.add_feature(cfeature.LAND.with_scale('50m'))
ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))
#ax.add_feature(cfeature.LAKES.with_scale('50m'))

print("--> Plot SPC polygons")
for key in cat_plot_colors.keys():
    geometries = cat_gdf[cat_gdf['LABEL2'] == key]
    # Check to see if there an area outlooked at all. If so, add the polygons to the map.
    if len(geometries) > 0:
        ax.add_geometries(geometries['geometry'], crs=data_crs,
                          facecolor=cat_plot_colors[key],
                          #edgecolor=cat_plot_colors[key],
                          alpha=0.5)

print("--> Legend, Title")
# Plot the legend
plt.legend(handles=legend_patches,loc=leg_loc)

plt.suptitle('Severe Storm Outlook', fontsize=18, fontweight='bold')
plt.title(f'{start_time} through {end_time}', fontsize=12, loc='center')

print("--> Save")
plt.savefig('spc/day1_categorical.png', dpi=96, bbox_inches='tight')

print("--> Done!")

#plt.show()

