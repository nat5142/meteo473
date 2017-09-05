##############################################################################################################
# Nicholas Tulli - nick.tulli95@gmail.com
# Script completed in April 2016 for METEO 473 - The Pennsylvania State University
# Purpose: Read HadCRUT global average temperature data for linear regression analysis of full dataset 
#          as well as abbreviated subsection of HadCRUT dataset
# File 'HadCRUT.4.4.0.0.annual_ns_avg.txt' provided for use by course administrators for METEO 473
#
# Script written for Python 3.4.0 on Department of Atmospheric Sciences and Meteorology servers
##############################################################################################################

# module library
import numpy as np
import random
import matplotlib.pyplot as plt
import random

# open and read HadCRUT file
f = open('HadCRUT.4.4.0.0.annual_ns_avg.txt','r')
year = []
temp = []
i = 1

# loop through lines in file to create arrays with appropriate data
for l in f:
    row = l.split()
    year.append(row[0])
    temp.append(row[1])
    years = [int(x) for x in year]
    temps = [float(y) for y in temp]

# run linear regression test to determine slope y-intercept, z-score and p-value of trendline
z = np.polyfit(years,temps,1)
m,b = np.polyfit(years,temps,1)
p = np.poly1d(z)

# plot figure of global average temperature departure from normal
plt.figure(1, figsize=(12,8))
plt.plot(years,temps)
plt.plot(years,p(years), "r--")
plt.grid()
plt.xlabel('Years')
plt.ylabel('Temperature Departure from 1966-1990 Average (Celsius)')
plt.xticks(years[0::10])
plt.hlines(0,1850,2016, "gray")
plt.title('Global Average Temperature Departure from Normal (Full Dataset')
plt.text(1860,1.02, 'Slope of Trendline = %f degrees C per year' % (m), fontsize = 12)
plt.savefig('MileStone3_figure1.png')
plt.show()
plt.close()

# run linear regression test on narrowed portion of extracted data
z_short = np.polyfit(years[117:167],temps[117:167],1)
m_short,b_short = np.polyfit(years[117:167],temps[117:167],1)
p_short = np.poly1d(z_short)

# plot figure of global average temperature departure from normal in narrowed subset of years
plt.figure(2, figsize=(12,8))
plt.plot(years[117:167],temps[117:167])
plt.plot(years[117:167],p_short(years[117:167]), "r--")
plt.grid()
plt.xlabel('Years')
plt.ylabel('Temperature Departure from 1966-1990 Average (Celsius)')
plt.xticks(years[0::5])
plt.hlines(0,1966,2016, "gray")
plt.title('Global Average Temperature Departure from Normal (1966-2016)')
plt.text(1970,1.02, 'Slope of Trendline = %f degrees C per year' % (m_short), fontsize = 12)
plt.savefig('MileStone3_figure2.png')
plt.show()
plt.close()

'''
# prompt user for input of statistical values for generation of random noise to be applied to data
mu = 0.0
sigma = float(input('Input a value for standard deviation of temperature noise: '))
alpha = float(input('Input a value for the noise parameter, alpha: '))

# create gaussian distribution of data given 'sigma' input
gauss = np.random.normal(mu,sigma, 50)

temp_short = temps[117:167]
slope_array = []
num_samples = 500

# begin loop to create gaussian distribution
for g in range(num_samples):
    gaussian = np.random.normal(mu,sigma, 50)
    gauss = gaussian
    i = 2
    for n in range(len(gauss) - 1):
        gauss = gaussian + (alpha * gaussian[i-1])
        i+=1
    temperatures = temp_short + gauss
    temp_slope, temp_intercept = np.polyfit(years[117:167],temperatures,1)
    temp_slope = float(temp_slope)
    slope_array.append(temp_slope)
#print(slope_array)


min_slope = min(slope_array)
max_slope = max(slope_array)

x_location = np.linspace(min_slope,max_slope, num = 100)
x_loc = x_location[1]

slopes = slope_array
slopes = np.split(slopes,x_location)
nbins = 75
y_loc = (num_samples/nbins) + 5

xticks = np.linspace(min_slope,max_slope, num = 5)
text_string = ('Alpha = %f \n' \
               'Sigma = %f \n' \
               'Trend for perfect case = %f \n' \
               'Trend for perfect case shown by red vertical dashed bar \n' \
               'Time Period: 1966-2016' % (alpha,sigma,m_short))

# display figure of distribution of randomly generated noise
plt.figure(3, figsize=(20,10))
plt.hist(slope_array, bins=nbins)
plt.axvline(m_short,0,30,color="red",linestyle="--")
plt.xticks(xticks)
plt.grid()
plt.xlabel('Trends of Global Mean temperature with 500 Generations of Random Noise')
plt.ylabel('Frequency of Trend Values')
plt.title('Trends in Global Mean Temperature from 1966-2016 in degrees Celsius', y=1.03)
plt.text(x_loc,y_loc,text_string,bbox=dict(facecolor='white',alpha=0.5))
plt.show()
plt.close()
'''



