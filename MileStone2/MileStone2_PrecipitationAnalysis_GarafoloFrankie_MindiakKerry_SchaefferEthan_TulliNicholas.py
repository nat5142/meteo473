# MileStone 2
#*********************************************************************
# Authors: Frankie Garafolo, Kerry Mindiak, Ethan Schaeffer and Nick Tulli
# Emails: fjg5053@psu.edu, kvm5478@psu.edu, eds5206@psu.edu and nat5142@psu.edu
#*********************************************************************
# This script will gather montly averaged precipitation data over 1979-
# 2016 from Earth System Research Laboratory (ESRL) for the whole world.
# Then it will isolate the precipitation data for Italy, and run a 
# regression model at each grid point corresponding to Italy to estimate
# trends in seasonal preciptation over the years 1979-2016. 

# Four sub maps on a first figure will be created depicting seasonal 
# spatial average precipitation over the 1979-2016 for Italy. 

# Four other sub maps, on the second figure, will show spatial trend of 
# the precipitation for each season over Italy. 

# The third figure will contain p-values at each gridpoint corresponding 
# to the trends in Figure 2. 

# Finally, four sub plots on the fourth figure will illustrate time 
# series of the averaged precipitation over all of the grid points in 
# Italy for each season.

# There will be a total of four figures with four subplots each.
#*********************************************************************
import netCDF4 as nc
import ncdump
import datetime as dt
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import shapefile
from scipy import stats
import numpy.ma as ma
from pylab import *

## Extract data from website 
## (http://www.esrl.noaa.gov/psd/data/gridded/data.cmap.html) and inspect
# We chose the enhanced data because all gaps are filled
url = 'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/cmap/enh/precip.mon.mean.nc'
url = r'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/cmap/enh/precip.mon.mean.nc'
f = nc.Dataset(url)
# ncdump.ncdump(f)   would display the variables and dimension information

## Define variables from data
longitude = f.variables['lon'][:]
latitude = f.variables['lat'][:]
time = f.variables['time'][:]
precip = f.variables['precip'][:]

## Focus on Italy longitudes and latitudes
# Italy is between 32N to 48N and 6E to 21E
Italy_lat = latitude[17:23]
Italy_lon = longitude[2:9]
[nt, nlat, nlon] = precip.shape

# Change the time from hours since 1800-1-1 to dates in YYYY-DD-MM
dt_time = [dt.date(1800, 1, 1) + dt.timedelta(hours=t)\
           for t in time]

#*********************************************************************

#*********************************************************************
## Read downloaded Italy Shapefile for use in plots
## Website: http://www.gadm.org
sf = shapefile.Reader("ITA_adm0.shp")
shapes = sf.shapes()

## Iterate through the elements of the file
#for name in dir(shapes[0]):
#    if not name.startswith('__'):
#        print(name)

## Get different attributes from shape file
bbox = shapes[0].bbox # boundary coordinates 
points = shapes[0].points # actual border coordinates
bbox = np.array(bbox)

# Make latitude
bbox_lat = []
for j in range(len(latitude)):
    if latitude[j] > bbox[1]-2.5 and latitude[j] < bbox[3]+2.5:
        # Add and subtract 2.5 in order to capture all of Italy
        bbox_lat.append(latitude[j])

# Make longitude
bbox_lon = []
for i in range(len(longitude)):
    if longitude[i] > bbox[0]-2.5 and longitude[i] < bbox[2]+2.5:
        # Add and subtract 2.5 in order to capture all of Italy
        bbox_lon.append(longitude[i])

## Plot shapefile
plt.figure(1)

## Setup basemap projection
m = Basemap(llcrnrlon=bbox[0],llcrnrlat=bbox[1],urcrnrlon=bbox[2], \
           urcrnrlat=bbox[3], resolution='c', projection='cyl', \
           lat_0 = 41.8719, lon_0 = 12.5674)
lon, lat = np.meshgrid(bbox_lon, bbox_lat)
xi, yi = m(lon, lat)

# Set contour levels from 0 to 4 and set colorbar ticks
clevs = np.linspace(0, 4, 400)
clev_ticks = np.linspace(0, 4, 5)

#*********************************************************************

#*********************************************************************
## Create arrays for Spring (March, April, May), Summer (June, July, August)
## Fall (September, October, November) and Winter (December, January, February)
precip_array = np.array(precip)

## Loop through precip array and store fall values
ii = 0

# Create individual arrays for September, October and November for 37 years,
# 7 latitudes and 6 longitudes
fall_sep=np.zeros((37,7,6))
fall_oct=np.zeros((37,7,6))
fall_nov=np.zeros((37,7,6))

# Begin loop
for i in range(8, len(precip)-6, 12):
    jj = 0
    for j in range(len(latitude)):
        if latitude[j] > bbox[1]-2.5 and latitude[j] < bbox[3]+2.5:            
            kk = 0
            for k in range(len(longitude)):
                if longitude[k] > bbox[0]-2.5 and longitude[k] < bbox[2]+2.5:
                    fall_sep[ii,jj,kk] = precip_array[i, j, k]
                    fall_oct[ii,jj,kk] = precip_array[i+1, j, k]
                    fall_nov[ii,jj,kk] = precip_array[i+2, j, k] 
                    kk += 1
            jj += 1
    ii += 1

# Put all fall months together into one season array
fall = np.concatenate((fall_sep, fall_oct, fall_nov), axis=0)

## Loop through precip array and store summer values
ii = 0

# Create individual arrays for June, July and August for 37 years,
# 7 latitudes and 6 longitudes
summer_jun=np.zeros((37,7,6))
summer_jul=np.zeros((37,7,6))
summer_aug=np.zeros((37,7,6))

# Begin loop
for i in range(5, len(precip)-12, 12):
    jj = 0
    for j in range(len(latitude)):
        if latitude[j] > bbox[1]-2.5 and latitude[j] < bbox[3]+2.5:            
            kk = 0
            for k in range(len(longitude)):
                if longitude[k] > bbox[0]-2.5 and longitude[k] < bbox[2]+2.5:
                    summer_jun[ii,jj,kk] = precip_array[i, j, k]
                    summer_jul[ii,jj,kk] = precip_array[i+1, j, k]
                    summer_aug[ii,jj,kk] = precip_array[i+2, j, k] 
                    kk += 1
            jj += 1
    ii += 1

# Put all summer months together into one season array
summer = np.concatenate((summer_jun, summer_jul, summer_aug), axis=0)

# Loop through precip array and store spring values
ii = 0

# Create individual arrays for March, April and May for 37 years,
# 7 latitudes and 6 longitudes
spring_mar=np.zeros((37,7,6))
spring_apr=np.zeros((37,7,6))
spring_may=np.zeros((37,7,6))

# Begin loop
for i in range(2, len(precip)-12, 12):
    jj = 0
    for j in range(len(latitude)):
        if latitude[j] > bbox[1]-2.5 and latitude[j] < bbox[3]+2.5:            
            kk = 0
            for k in range(len(longitude)):
                if longitude[k] > bbox[0]-2.5 and longitude[k] < bbox[2]+2.5:
                    spring_mar[ii,jj,kk] = precip_array[i, j, k]
                    spring_apr[ii,jj,kk] = precip_array[i+1, j, k]
                    spring_may[ii,jj,kk] = precip_array[i+2, j, k] 
                    kk += 1
            jj += 1
    ii += 1

# Put all spring months together into one season array
spring = np.concatenate((spring_mar, spring_apr, spring_may),axis=0)

# Loop through precip array and store winter values
ii = 0

# Create individual arrays for December, January and February for 37 years,
# 7 latitudes and 6 longitudes
winter_dec=np.zeros((37,7,6))
winter_jan=np.zeros((37,7,6))
winter_feb=np.zeros((37,7,6))

# Begin loop
for i in range(11, len(precip)-3, 12):
    jj = 0
    for j in range(len(latitude)):
        if latitude[j] > bbox[1]-2.5 and latitude[j] < bbox[3]+2.5:            
            kk = 0
            for k in range(len(longitude)):
                if longitude[k] > bbox[0]-2.5 and longitude[k] < bbox[2]+2.5:
                    winter_dec[ii,jj,kk] = precip_array[i, j, k]
                    winter_jan[ii,jj,kk] = precip_array[i+1, j, k]
                    winter_feb[ii,jj,kk] = precip_array[i+2, j, k] 
                    kk += 1
            jj += 1
    ii += 1        

# Put all winter months together into one season array
winter = np.concatenate((winter_dec, winter_jan, winter_feb), axis=0)

#*********************************************************************

#*********************************************************************
## Figure 1: Four sub maps depicting seasonal spatial average 
## precipitation over the 1979-2016 period for Italy
#*********************************************************************
## Winter
# Put in the 2x2 plot in the top left corner (1)
plt.subplot(2, 2, 1)

# Calculate average precipitation for each grid point
ii = 0
# Make array of 7 latitudes and 6 longitudes
winter_avg = np.zeros((7,6))

# Begin loop
for j in range(7):
    kk = 0
    for k in range(6):
        winter_avg[ii,kk] = np.mean(winter[:,j,k])
        kk += 1
    ii += 1

# Draw map with contours
m.drawmapboundary(fill_color='white')
m.fillcontinents(color='w',lake_color='white',zorder=0)
# Read shapefile
m.readshapefile('ITA_adm0','ITA_adm0')
# Make contours
cs = m.contourf(xi, yi, winter_avg, clevs)
# Make colorbar
cbar = m.colorbar(cs,location='bottom',pad="10%",ticks=clev_ticks)
cbar.set_label('Precipitation (mm/day)')
# Set plot title
plt.title('Winter')

## Spring 
# Put in the 2x2 plot in the top right corner (2)       
plt.subplot(2, 2, 2)

# Calculate average precipitation for each grid point
ii = 0
# Make array of 7 latitudes and 6 longitudes
spring_avg = np.zeros((7,6))

# Begin loop
for j in range(7):
    kk = 0
    for k in range(6):
        spring_avg[ii,kk] = np.mean(spring[:,j,k])
        kk += 1
    ii += 1

# Draw map with contours
m.drawmapboundary(fill_color='white')
m.fillcontinents(color='w',lake_color='white',zorder=0)
# Read shapefile
m.readshapefile('ITA_adm0','ITA_adm0')
# Make contours
cs = m.contourf(xi, yi, spring_avg, clevs)
# Make colorbar
cbar = m.colorbar(cs,location='bottom',pad="10%",ticks=clev_ticks)
cbar.set_label('Precipitation (mm/day)')
# Set plot title
plt.title('Spring')

## Summer
# Put in the 2x2 plot in the bottom left corner (3)
plt.subplot(2, 2, 3)

# Calculate average precipitation for each grid point
ii = 0
# Make array of 7 latitudes and 6 longitudes
summer_avg = np.zeros((7,6))

# Begin loop
for j in range(7):
    kk = 0
    for k in range(6):
        summer_avg[ii,kk] = np.mean(summer[:,j,k])
        kk += 1
    ii += 1

# Draw map with contours
m.drawmapboundary(fill_color='white')
m.fillcontinents(color='w',lake_color='white',zorder=0)
# Read shapefile
m.readshapefile('ITA_adm0','ITA_adm0')
# Make contours
cs = m.contourf(xi, yi, summer_avg, clevs)
# Make colorbar
cbar = m.colorbar(cs,location='bottom',pad="10%",ticks=clev_ticks)
cbar.set_label('Precipitation (mm/day)')
# Set plot title
plt.title('Summer')

## Fall
# Put in the 2x2 plot in the bottom right corner (4)
plt.subplot(2, 2, 4)

# Calculate average precipitation for each grid point
ii = 0
# Make array of 7 latitudes and 6 longitudes
fall_avg = np.zeros((7,6))

# Begin loop
for j in range(7):
    kk = 0
    for k in range(6):
        fall_avg[ii,kk] = np.mean(fall[:,j,k])
        kk += 1
    ii += 1

# Draw map with contours
m.drawmapboundary(fill_color='white')
m.fillcontinents(color='w',lake_color='white',zorder=0)
# Read shapefile
m.readshapefile('ITA_adm0','ITA_adm0')
# Make contours
cs = m.contourf(xi, yi, fall_avg, clevs)
# Make colorbar
cbar = m.colorbar(cs,location='bottom',pad="10%",ticks=clev_ticks)
cbar.set_label('Precipitation (mm/day)')
# Set plot title
plt.title('Fall')

# Make title for entire figure and show plot
# Add space between subplots
plt.subplots_adjust(hspace=.5, wspace=.5)
plt.suptitle('Seasonal Spatially Averaged Precipitation from 1979 to 2016', fontsize=15)
plt.savefig("MileStone2_figure1.png")
plt.show()
#*********************************************************************
######################################################################
#*********************************************************************
## Linear regression trend line for the seasonal data
#*********************************************************************
## Create seasonal average for each season to be used in the linear
## regression

## Spring
# Create array of 37 years, 7 latitudes and 6 longitudes
spring_seasonal_avg = np.zeros((37,7,6))

# Begin loop for seasonal average
ii = 0
for i in range(0, len(spring), 3):
    jj = 0
    for j in range(7):
        kk = 0
        for k in range(6):
            spring_seasonal_avg[ii,jj,kk] = (spring[i, j, k] + spring[i+1, j, k]\
            + spring[i+2, j, k])/3
            kk += 1
        jj += 1
    ii += 1

## Summer
# Create array of 37 years, 7 latitudes and 6 longitudes
summer_seasonal_avg = np.zeros((37,7,6))

# Begin loop for seasonal average
ii = 0
for i in range(0, len(summer), 3):
    jj = 0
    for j in range(7):
        kk = 0
        for k in range(6):
            summer_seasonal_avg[ii,jj,kk] = (summer[i, j, k] + summer[i+1, j, k]\
            + summer[i+2, j, k])/3
            kk += 1
        jj += 1
    ii += 1

## Fall 
# Create array of 37 years, 7 latitudes and 6 longitudes
fall_seasonal_avg = np.zeros((37,7,6))

# Begin loop for seasonal average
ii = 0
for i in range(0, len(fall), 3):
    jj = 0
    for j in range(7):
        kk = 0
        for k in range(6):
            fall_seasonal_avg[ii,jj,kk] = (fall[i, j, k] + fall[i+1, j, k]\
            + fall[i+2, j, k])/3
            kk += 1
        jj += 1
    ii += 1

## Winter
# Create array of 37 years, 7 latitudes and 6 longitudes
winter_seasonal_avg = np.zeros((37,7,6))

# Begin loop for seasonal average
ii = 0
for i in range(0, len(winter), 3):
    jj = 0
    for j in range(7):
        kk = 0
        for k in range(6):
            winter_seasonal_avg[ii,jj,kk] = (winter[i, j, k] + winter[i+1, j, k]\
            + winter[i+2, j, k])/3
            kk += 1
        jj += 1
    ii += 1

#*********************************************************************
## Create variables and values for linear regression for each season
#*********************************************************************
# Create time variable
time_y = []
w=1979
for i in range(37):
    time_y = np.append(time_y, w+i)

## Spring
# Create arrays for slope, intercept, r-value, p-value and standard deviation
# from stats.linregress of size (latitude, longitude)
slope_spring = np.zeros((7,6))
intercept_spring = np.zeros((7,6))
rv_spring = np.zeros((7,6))
pv_spring = np.zeros((7,6))
stderr_spring = np.zeros((7,6))

# Begin loop
for j in range(7):
    for k in range(6):
        slope_spring[j,k], intercept_spring[j,k], rv_spring[j,k], pv_spring[j,k], \
        stderr_spring[j,k] = stats.linregress(time_y, spring_seasonal_avg[:,j,k])

# Summer
# Create arrays for slope, intercept, r-value, p-value and standard deviation
# from stats.linregress of size (latitude, longitude)
slope_summer = np.zeros((7,6))
intercept_summer = np.zeros((7,6))
rv_summer = np.zeros((7,6))
pv_summer = np.zeros((7,6))
stderr_summer = np.zeros((7,6))

# Begin loop
for j in range(7):
    for k in range(6):
        slope_summer[j,k], intercept_summer[j,k], rv_summer[j,k], pv_summer[j,k], \
        stderr_summer[j,k] = stats.linregress(time_y, summer_seasonal_avg[:,j,k])

# Fall
# Create arrays for slope, intercept, r-value, p-value and standard deviation
# from stats.linregress of size (latitude, longitude)
slope_fall = np.zeros((7,6))
intercept_fall = np.zeros((7,6))
rv_fall = np.zeros((7,6))
pv_fall = np.zeros((7,6))
stderr_fall = np.zeros((7,6))

# Begin loop
for j in range(7):
    for k in range(6):
        slope_fall[j,k], intercept_fall[j,k], rv_fall[j,k], pv_fall[j,k], \
        stderr_fall[j,k] = stats.linregress(time_y, fall_seasonal_avg[:,j,k])

# Winter
# Create arrays for slope, intercept, r-value, p-value and standard deviation
# from stats.linregress of size (latitude, longitude)
slope_winter = np.zeros((7,6))
intercept_winter = np.zeros((7,6))
rv_winter = np.zeros((7,6))
pv_winter = np.zeros((7,6))
stderr_winter = np.zeros((7,6))

# Begin loop
for j in range(7):
    for k in range(6):
        slope_winter[j,k], intercept_winter[j,k], rv_winter[j,k], pv_winter[j,k], \
        stderr_winter[j,k] = stats.linregress(time_y, winter_seasonal_avg[:,j,k])

#*********************************************************************
## Plot Figure 2: Four maps with the spatial trend of precipitation using
## slope_SEASON
#*********************************************************************

plt.figure(2)

## Define contour levels from -.1 to .1 and set ticks
clevs = np.linspace(-.1, .1, 100)
clev_ticks = np.linspace(-.1, .1, 5)
 
## Winter linear regression slope
# Put in the 2x2 plot in the top left corner (1)
plt.subplot(2,2,1)

# Draw map with contours
m.drawmapboundary(fill_color='white')
m.fillcontinents(color='w',lake_color='white',zorder=0)
# Read shapefile
m.readshapefile('ITA_adm0','ITA_adm0')
# Make contours
cs = m.contourf(xi, yi, slope_winter, clevs)
# Make colorbar
cbar = m.colorbar(cs,location='bottom',pad="10%",ticks=clev_ticks)
cbar.set_label('Precipitation (mm/day)')
# Set plot title
plt.title('Winter')

## Spring linear regression slope
# Put in the 2x2 plot in the top right corner (2)       
plt.subplot(2,2,2)

# Draw map with contours
m.drawmapboundary(fill_color='white')
m.fillcontinents(color='w',lake_color='white',zorder=0)
# Read shapefile
m.readshapefile('ITA_adm0','ITA_adm0')
# Make contours
cs = m.contourf(xi, yi, slope_spring, clevs)
# Make colorbar
cbar = m.colorbar(cs,location='bottom',pad="10%",ticks=clev_ticks)
cbar.set_label('Precipitation (mm/day)')
# Set plot title
plt.title('Spring')

## Summer linear regression slope
# Put in the 2x2 plot in the bottom left corner (3)
plt.subplot(2,2,3)

# Draw map with contours
m.drawmapboundary(fill_color='white')
m.fillcontinents(color='w',lake_color='white',zorder=0)
# Read shapefile
m.readshapefile('ITA_adm0','ITA_adm0')
# Make contours
cs = m.contourf(xi, yi, slope_summer, clevs)
# Make contours
cbar = m.colorbar(cs,location='bottom',pad="10%",ticks=clev_ticks)
cbar.set_label('Precipitation (mm/day)')
# Set plot title
plt.title('Summer')

## Fall linear regression slope
# Put in the 2x2 plot in the bottom right corner (4)
plt.subplot(2,2,4)

# Draw map with contours
m.drawmapboundary(fill_color='white')
m.fillcontinents(color='w',lake_color='white',zorder=0)
# Read shapefile
m.readshapefile('ITA_adm0','ITA_adm0')
# Make contours
cs = m.contourf(xi, yi, slope_fall, clevs)
# Make colorbar
cbar = m.colorbar(cs,location='bottom',pad="10%",ticks=clev_ticks)
cbar.set_label('Precipitation (mm/day)')
# Set plot title
plt.title('Fall')

## Make title for entire figure and show plot
# Add space between subplots
plt.subplots_adjust(hspace=.5, wspace=.5)
plt.suptitle('Linear Trend for Yearly Seasonal Average Precipitation', fontsize=15)
plt.savefig("MileStone2_figure2.png")
plt.show()

#*********************************************************************
######################################################################
#*********************************************************************
## Figure 3: Four maps with p-values at each grid point corresponding
## to trends from Figure 2
#*********************************************************************
# p-values are the pv_SEASON from linear regression

plt.figure(3)

## Define contour levels from 0 to 1 and set colorbar ticks
clevs = np.linspace(0, 1, 100)
clev_ticks = np.linspace(0, 1, 5)

## Winter
# Put in the 2x2 plot in the top left corner (1)
plt.subplot(2,2,1)

# Draw map with contours
m.drawmapboundary(fill_color='white')
m.fillcontinents(color='w',lake_color='white',zorder=0)
# Read shapefile
m.readshapefile('ITA_adm0','ITA_adm0')
# Make contours
cs = m.contourf(xi, yi, pv_winter, clevs)
# Make colorbar
m.colorbar(cs,location='bottom',pad="10%",ticks=clev_ticks)
# Set plot title
plt.title('Winter')

## Spring
# Put in the 2x2 plot in the top right corner (2)       
plt.subplot(2,2,2)

# Draw map with contours
m.drawmapboundary(fill_color='white')
m.fillcontinents(color='w',lake_color='white',zorder=0)
# Read shapefile
m.readshapefile('ITA_adm0','ITA_adm0')
# Make contours
cs = m.contourf(xi, yi, pv_spring, clevs)
# Make colorbar
m.colorbar(cs,location='bottom',pad="10%",ticks=clev_ticks)
# Set plot title
plt.title('Spring')

## Summer
# Put in the 2x2 plot in the bottom left corner (3)
plt.subplot(2,2,3)

# Draw map with contours
m.drawmapboundary(fill_color='white')
m.fillcontinents(color='w',lake_color='white',zorder=0)
# Read shapefile
m.readshapefile('ITA_adm0','ITA_adm0')
# Make contours
cs = m.contourf(xi, yi, pv_summer, clevs)
# Make colorbar
m.colorbar(cs,location='bottom',pad="10%",ticks=clev_ticks)
# Set plot title
plt.title('Summer')

## Fall
# Put in the 2x2 plot in the bottom right corner (4)
plt.subplot(2,2,4)

# Draw map with contours
m.drawmapboundary(fill_color='white')
m.fillcontinents(color='w',lake_color='white',zorder=0)
# Read shapefile
m.readshapefile('ITA_adm0','ITA_adm0')
# Make contours
cs = m.contourf(xi, yi, pv_fall, clevs)
# Make colorbar
m.colorbar(cs,location='bottom',pad="10%",ticks=clev_ticks)
# Set plot title
plt.title('Fall')

# Add space between subplots
plt.subplots_adjust(hspace=.5, wspace=.5)
# Make title for entire figure and show plot
plt.suptitle('P-values Corresponding to Trends in Figure 2', fontsize=15)
plt.savefig("MileStone2_figure3.png")
plt.show()

#*********************************************************************
######################################################################
#*********************************************************************
## Figure 4: Four subplots for each season illustrating the time series of
## averaged precipitation over all grid points, with a corresponding linear 
## trend line and text box with: slope of the regression line, intercept of
## the regression line, correlation coefficient, two-sided p-value for a 
## hypothesis test whose null hypothesis is that the slope is zero, standard
## error of the estimate, and coefficient of determination.
#*********************************************************************

plt.figure(4, figsize=(20,10))

# Set year array from 1979 to 2016 with 37 steps
year = np.linspace(1979, 2016, 37)

## Winter
# Put in the 2x2 plot in the top left corner (1)
plt.subplot(2,2,1)

# Make year average of length 37 for 37 years of data
winter_year_avg = np.zeros((37))

# Begin loop
for i in range(37):
    winter_year_avg[i] = np.mean(winter_seasonal_avg[i, :, :])

# Set variables: slope, intercept, r-value, p-value and std-err equal to 
# stats.linregress
winter_slope, winter_intercept, winter_r_value, winter_p_value, winter_std_err = \
stats.linregress(winter_year_avg, year)

# Plot trendline
plt.plot(year, winter_year_avg)
(m,b) = np.polyfit(year,winter_year_avg,1)
yp = np.polyval([m,b],year)
plt.plot(year, yp)

# Create correlation coefficient 
winter_corr = np.corrcoef(year, winter_year_avg)

# Create R^2 by squaring r_value
winter_r = winter_r_value**2

# Create textbox with variables: slope, intercept, correlation coefficient, p-value
# r-value, SEE, and R^2
plt.text(1976, 1.0, ' Slope = %.2f\n Intercept  = %.2f\n Correlation Coefficient = %.2f\n P-Value = %.2f\n SEE = %.2f\n R^2 = %.2f\n' % (winter_slope, winter_intercept, \
winter_corr[0,1], winter_p_value, winter_std_err, winter_r), fontsize = 8)

plt.xlabel('Year', fontsize=15)
# Rotate xlabels
plt.xticks(rotation=70)
plt.ylabel('Average Precip (mm/day)', fontsize=15)
# Set plot title
plt.title('Winter', fontsize=15)

## Spring
# Put in the 2x2 plot in the top right corner (2)       
plt.subplot(2,2,2)

# Make year average of length 37 for 37 years of data
spring_year_avg = np.zeros((37))

# Begin loop
for i in range(37):
    spring_year_avg[i] = np.mean(spring_seasonal_avg[i, :, :])

# Set variables: slope, intercept, r-value, p-value and std-err equal to 
# stats.linregress
spring_slope, spring_intercept, spring_r_value, spring_p_value, spring_std_err = \
stats.linregress(spring_year_avg, year)

# Plot trendline
plt.plot(year, spring_year_avg)
(m,b) = np.polyfit(year,spring_year_avg,1)
yp = np.polyval([m,b],year)
plt.plot(year, yp)

# Create correlation coefficient 
spring_corr = np.corrcoef(year, spring_year_avg)

# Create R^2 by squaring r_value
spring_r = spring_r_value**2

# Create textbox with variables: slope, intercept, correlation coefficient, p-value
# r-value, SEE, and R^2
plt.text(1976, 1.0, ' Slope = %.2f\n Intercept =  %.2f\n Correlation Coefficient = %.2f\n \P-Value = %.2f\n SEE = %.2f\n R^2 = %.2f\n' % (spring_slope, spring_intercept, spring_corr[0,1],\
spring_p_value, spring_std_err, spring_r), fontsize = 8)

plt.xlabel('Year', fontsize=15)
# Rotate xlabels
plt.xticks(rotation=70)
plt.ylabel('Average Precip (mm/day)', fontsize=15)
# Set plot title
plt.title('Spring', fontsize=15)

## Summer
# Put in the 2x2 plot in the bottom left corner (3)
plt.subplot(2,2,3)

# Make year average of length 37 for 37 years of data
summer_year_avg = np.zeros((37))

# Begin loop
for i in range(37):
    summer_year_avg[i] = np.mean(summer_seasonal_avg[i, :, :])

# Set variables: slope, intercept, r-value, p-value and std-err equal to 
# stats.linregress
summer_slope, summer_intercept, summer_r_value, summer_p_value, summer_std_err = \
stats.linregress(summer_year_avg, year)

# Plot trendline
plt.plot(year, summer_year_avg)
(m,b) = np.polyfit(year,summer_year_avg,1)
yp = np.polyval([m,b],year)
plt.plot(year, yp)

# Create correlation coefficient 
summer_corr = np.corrcoef(year, summer_year_avg)

# Create R^2 by squaring r_value
summer_r = summer_r_value**2

# Create textbox with variables: slope, intercept, correlation coefficient, p-value
# r-value, SEE, and R^2
plt.text(1976, 0.8, ' Slope = %.2f\n Intercept =  %.2f\n Correlation Coefficient = %.2f\n \P-Value = %.2f\n SEE = %.2f\n R^2 = %.2f\n' % (summer_slope, summer_intercept, summer_corr[0,1],\
summer_p_value, summer_std_err, summer_r), fontsize = 8)

plt.xlabel('Year', fontsize=15)
# Rotate xlabels
plt.xticks(rotation=70)
plt.ylabel('Average Precip (mm/day)', fontsize=15)
# Set plot title
plt.title('Summer', fontsize=15)

## Fall
# Put in the 2x2 plot in the bottom right corner (4)
plt.subplot(2,2,4)

# Make year average of length 37 for 37 years of data
fall_year_avg = np.zeros((37))

# Begin loop
for i in range(37):
    fall_year_avg[i] = np.mean(fall_seasonal_avg[i, :, :])

# Set variables: slope, intercept, r-value, p-value and std-err equal to 
# stats.linregress
fall_slope, fall_intercept, fall_r_value, fall_p_value, fall_std_err = \
stats.linregress(fall_year_avg, year)

# Plot trendline
plt.plot(year, fall_year_avg)
(m,b) = np.polyfit(year,fall_year_avg,1)
yp = np.polyval([m,b],year)
plt.plot(year, yp)

# Create correlation coefficient 
fall_corr = np.corrcoef(year, fall_year_avg)

# Create R^2 by squaring r_value
fall_r = fall_r_value**2

# Create textbox with variables: slope, intercept, correlation coefficient, p-value
# r-value, SEE, and R^2
plt.text(2007, 1.0, ' Slope = %.2f\n Intercept =  %.2f\n Correlation Coefficient = %.2f\n \P-Value = %.2f\n SEE = %.2f\n R^2 = %.2f\n' % (fall_slope, fall_intercept, fall_corr[0,1],\
fall_p_value, fall_std_err, fall_r), fontsize = 8)

plt.xlabel('Year', fontsize=15)
# Rotate xlabels
plt.xticks(rotation=70)
plt.ylabel('Average Precip (mm/day)', fontsize=15)
# Set plot title
plt.title('Fall', fontsize=15)

## Make title for entire figure and show plot
# Add space between subplots
plt.subplots_adjust(hspace=.5, wspace=.5)

# Make title for entire figure and show plot
plt.suptitle('Time Series of Averaged Precipitation with Corresponding Linear Trendline', fontsize=15)
plt.savefig("MileStone2_figure4.png")
plt.show()

#*********************************************************************
######################################################################
#*********************************************************************
# DONE WITH MILESTONE 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! :D
#*********************************************************************
