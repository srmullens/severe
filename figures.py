#!/usr/bin/python

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import from_levels_and_colors
import numpy as np
import urllib3

#################
# Get some Data #
#################
year_start = 2003
year_end = 2019

##########################
# Data for First Subplot #
##########################
# List of max outlooks
outlooks_list = [[] for x in range(year_end-year_start+1)]
implied_outlooks_list = [[] for x in range(year_end-year_start+1)]
tornadoes_list = [[] for x in range(year_end-year_start+1)]
hail_list = [[] for x in range(year_end-year_start+1)]
wind_list = [[] for x in range(year_end-year_start+1)]

######################
# Defining Functions #
######################
def report_progress(year,month,day):
    if day==1: print(str(year)+str(month).zfill(2)+str(day).zfill(2))
    #print str(year)+str(month).zfill(2)+str(day).zfill(2)

###################
# For each day... #
###################
for year in range(year_start,year_end+1):
    for month in range(1,13):
        for day in range(1,32):
            if day==31 and month in (2,4,6,9,11): continue
            if day==30 and month==2: continue

            # Search these SPC Convective Outlook URLs...
            url1 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year)+'/KWNSPTSDY1_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'1200.txt'
            url2 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year)+'/KWNSPTSDY1_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'1300.txt'
            url3 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year)+'/KWNSPTSDY1_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'1630.txt'
            url4 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year)+'/KWNSPTSDY1_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'2000.txt'
            if day==31 and month==12:
              url5 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year+1)+'/KWNSPTSDY1_'+str(year+1)+'01010100.txt'
            elif day==28 and month==2 and year%4==0:
                url5 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year)+'/KWNSPTSDY1_'+str(year)+str(month).zfill(2)+str(day+1).zfill(2)+'0100.txt'
            elif day==28 and month==2 and year%4>0:
                url5 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year)+'/KWNSPTSDY1_'+str(year)+'03010100.txt'
            elif (day==30 and month in (4,6,9,11)) or (day==29 and month==2 and year%4==0) or (day==31): 
              url5 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year)+'/KWNSPTSDY1_'+str(year)+str(month+1).zfill(2)+'010100.txt'	
            else:
              url5 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year)+'/KWNSPTSDY1_'+str(year)+str(month).zfill(2)+str(day+1).zfill(2)+'0100.txt'

            urls = [url1,url2,url3,url4,url5]

            # These URLs must be used prior to 24 March 2005. No TOR/HAIL/WIND percentage breakdowns.
            url6 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year)+'/day1otlk_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'_1200.html'
            url7 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year)+'/day1otlk_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'_1300.html'
            url8 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year)+'/day1otlk_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'_1630.html'
            url9 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year)+'/day1otlk_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'_2000.html'
            if day==31 and month==12:
              url10 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year+1)+'/day1otlk_'+str(year+1)+'0101_0100.html'
            elif day==28 and month==2 and year%4==0:
              url10 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year)+'/day1otlk_'+str(year)+str(month).zfill(2)+str(day+1).zfill(2)+'_0100.html'
            elif day==28 and month==2 and year%4>0:
              url10 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year)+'/day1otlk_'+str(year)+'0301_0100.html'
            elif (day==30 and month in (4,6,9,11)) or (day==29 and month==2 and year%4==0) or (day==31): 
              url10 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year)+'/day1otlk_'+str(year)+str(month+1).zfill(2)+str(1).zfill(2)+'_0100.html'	
            else:
              url10 = 'http://www.spc.noaa.gov/products/outlook/archive/'+str(year)+'/day1otlk_'+str(year)+str(month).zfill(2)+str(day+1).zfill(2)+'_0100.html'

            urls_old = [url6,url7,url8,url9,url10]

            if year>2005: url_old_or_new = 1
            elif year==2005 and month>3: url_old_or_new = 1
            elif year==2005 and month==3 and day>=24: url_old_or_new = 1
            else: url_old_or_new = 0

            # Find the maximum outlook level for each day...
            # First, new URLs with TOR/WIND/HAIL breakdown
            if url_old_or_new==1:
                report_progress(year,month,day)
                outlook = 0; max_outlook = -1
                implied_outlook = 0; max_implied_outlook = -1
                tor = 0;  sig_tor = 0;  max_tor = 0
                hail = 0; sig_hail = 0; max_hail = 0
                wind = 0; sig_wind = 0; max_wind = 0
                read = 0; read_tor = 0; read_hail = 0; read_wind = 0

                #print year,month,day

                for url in urls:
                    try:
                        #if month==12 and day==31: print url
                        for line in urllib3.urlopen(url):
                            outlook = 0; implied_outlook = 0
                            tor = 0; sig_tor = 0
                            hail = 0; sig_hail = 0
                            wind = 0; sig_wind = 0

                            #print "LINE",line
                            if "... CATEGORICAL ..." in line: read = 1
                            if read==1 and "&&" in line: read = 0; break
                            if read==1 and "HIGH" in line: outlook = 6; #print url[-8:-4],"HIGH"
                            if read==1 and "MDT" in line: outlook = 5; #print url[-8:-4],"MDT"
                            if read==1 and "ENH" in line: outlook = 4; #print url[-8:-4],"ENH"
                            if read==1 and "SLGT" in line: outlook = 3; #print url[-8:-4],"SLGT"
                            if read==1 and "MRGL" in line: outlook = 2; #print url[-8:-4],"MRGL"
                            if read==1 and "TSTM" in line: outlook = 1; #print url[-8:-4],"TSTM"
                            #if read==1: print "CAT",outlook,line

                            if "... TORNADO ..." in line: read_tor = 1
                            if read_tor==1 and "&&" in line: read_tor = 0
                            if read_tor==1 and "SIGN" in line: sig_tor = 1; #print "SIG TOR"
                            if read_tor==1 and "0.60" in line: tor = 7
                            if read_tor==1 and "0.45" in line: tor = 6
                            if read_tor==1 and "0.30" in line: tor = 5
                            if read_tor==1 and "0.15" in line: tor = 4
                            if read_tor==1 and "0.10" in line: tor = 3
                            if read_tor==1 and "0.05" in line: tor = 2
                            if read_tor==1 and "0.02" in line: tor = 1
                            #if read_tor==1: print "TOR",tor,line

                            if "... HAIL ..." in line: read_hail = 1
                            if read_hail==1 and "&&" in line: read_hail = 0
                            if read_hail==1 and "SIGN" in line: sig_hail = 1; #print "SIG HAIL"
                            if read_hail==1 and "0.60" in line: hail = 5
                            if read_hail==1 and "0.45" in line: hail = 4
                            if read_hail==1 and "0.30" in line: hail = 3
                            if read_hail==1 and "0.15" in line: hail = 2
                            if read_hail==1 and "0.05" in line: hail = 1
                            #if read_hail==1: print "HAIL",hail,line

                            if "... WIND ..." in line: read_wind = 1
                            if read_wind==1 and "&&" in line: read_wind = 0
                            if read_wind==1 and "SIGN" in line: sig_wind = 1; #print "SIG WIND"
                            if read_wind==1 and "0.60" in line: wind = 5
                            if read_wind==1 and "0.45" in line: wind = 4
                            if read_wind==1 and "0.30" in line: wind = 3
                            if read_wind==1 and "0.15" in line: wind = 2
                            if read_wind==1 and "0.05" in line: wind = 1
                            #if read_wind==1: print "WIND",wind,line

                            #Implied 2018 Category based on TOR/HAIL/WIND percentages
                            if max(tor,hail,wind)==1: implied_outlook = 2  # MRGL
                            if max(tor,hail,wind)==2: implied_outlook = 3  # SLGT
                            if max(tor,wind,hail)==3 or (max(tor+sig_tor,wind+sig_wind,hail+sig_hail)==4): implied_outlook = 4  # ENH 
                            if max(tor+sig_tor,wind+sig_wind,hail+sig_hail)==5: implied_outlook = 5  # MDT
                            if max(tor+sig_tor,wind+sig_wind) >= 6: implied_outlook = 6  # HIGH

                            max_outlook = max(max_outlook,outlook)
                            max_implied_outlook = max(max_implied_outlook,implied_outlook,outlook)
                            max_tor  = max(max_tor,tor)
                            max_hail = max(max_hail,hail)
                            max_wind = max(max_wind,wind)

                        #print month,day,max_outlook

                        diff = max_implied_outlook-max_outlook
                        if diff>1: print(year,month,day,diff)

                        #if read==0: print max_tor,max_hail,max_wind,max_outlook,max_implied_outlook

                    # If the website is missing that day.
                    except:
                        outlook = -1
                        implied_outlook = -1
                        tor = -1
                        hail = -1
                        wind = -1

                        max_outlook = max(max_outlook,outlook)
                        max_implied_outlook = max(max_implied_outlook,implied_outlook,outlook)
                        max_tor  = max(max_tor,tor)
                        max_hail = max(max_hail,hail)
                        max_wind = max(max_wind,wind)

                #if month==12 and day>28: print month,day,max_outlook

                # Record the maximum outlook level in a list...
                outlooks_list[year-year_start].append(max_outlook)
                implied_outlooks_list[year-year_start].append(max_implied_outlook)
                tornadoes_list[year-year_start].append(max_tor)
                hail_list[year-year_start].append(max_hail)
                wind_list[year-year_start].append(max_wind)
                #if month==12 and day==31: print "Added",max_outlook,year,month,day,len(outlooks_list[0])

            # Old URLs prior to TOR/WIND/HAIL breakdown
            elif url_old_or_new==0:
                report_progress(year,month,day)
                outlook = 0
                max_outlook = -1
                read = 0

                for url in urls_old:
                    # print url
                    try:
                        for line in urllib3.urlopen(url):
                            ### These need to change. ###
                            if "VALID" in line: read = 1; #print "VALID"
                            #if read==1 and "&&" in line: read = 0
                            if read==1 and "THERE IS A HIGH RISK" in line: outlook = 6; #print month,day,"Found HIGH"
                            if read==1 and "THERE IS A MDT RISK" in line: outlook = 5; #print month,day,"Found MDT"
                            #if read==1 and "ENH" in line: outlook = 4; break
                            if read==1 and "THERE IS A SLGT RISK" in line: outlook = 3; #print "Fount SLGT"; break
                            #if read==1 and "MRGL" in line: outlook = 2; break
                            if read==1 and "GEN TSTMS ARE FCST" in line: outlook = 1; #print "Found TSTMS"; break
                            max_outlook = max(max_outlook,outlook)

                    # If the website is missing that day.
                    except:
                        print(month,day)
                        outlook = -1
                        max_outlook = max(max_outlook,outlook)

                #if month==12 and day==31: print month,day,max_outlook

                # Record the maximum outlook level in a list...
                outlooks_list[year-year_start].append(max_outlook)
                implied_outlooks_list[year-year_start].append(max_outlook)
                #if month==12 and day==31: print "Added",max_outlook,year,month,day,len(outlooks_list[0])

            else:
                #print "Added -1",year,month,day
                outlooks_list[year-year_start].append(-1)
                implied_outlooks_list[year-year_start].append(-1)

#####################################
### Report results of 2D array... ###
#####################################
for i,n in enumerate(outlooks_list):
    a = outlooks_list[i].count(-1)
    b = outlooks_list[i].count(0)
    c = outlooks_list[i].count(1)
    d = outlooks_list[i].count(2)
    e = outlooks_list[i].count(3)
    f = outlooks_list[i].count(4)
    g = outlooks_list[i].count(5)
    h = outlooks_list[i].count(6)
    print(str(i)+": "+str(len(outlooks_list[i]))+" - ("+str(a)+", "+str(b)+", "+str(c)+", "+str(d)+", "+str(e)+", "+str(f)+", "+str(g)+", "+str(h)+")")

print("-----------------------------------------")

for i,n in enumerate(implied_outlooks_list):
    a = implied_outlooks_list[i].count(-1)
    b = implied_outlooks_list[i].count(0)
    c = implied_outlooks_list[i].count(1)
    d = implied_outlooks_list[i].count(2)
    e = implied_outlooks_list[i].count(3)
    f = implied_outlooks_list[i].count(4)
    g = implied_outlooks_list[i].count(5)
    h = implied_outlooks_list[i].count(6)
    print(str(i)+": "+str(len(implied_outlooks_list[i]))+" - ("+str(a)+", "+str(b)+", "+str(c)+", "+str(d)+", "+str(e)+", "+str(f)+", "+str(g)+", "+str(h)+")")

###########################
# Data for Second Subplot #
###########################
max_outlook_by_day = np.array([])

for day in range(len(outlooks_list[0])):
    max_outlook_by_day = np.append(max_outlook_by_day,-2)
    for row in range(year_end-year_start+1):
        outlook = outlooks_list[row][day]
        if max_outlook_by_day[day] < outlook: max_outlook_by_day[day] = outlook

max_outlook_by_day = np.expand_dims(max_outlook_by_day, axis=0)

# Repeat, but for implied_outlooks
max_implied_outlook_by_day = np.array([])

for day in range(len(implied_outlooks_list[0])):
    max_implied_outlook_by_day = np.append(max_implied_outlook_by_day,-2)
    for row in range(year_end-year_start+1):
        implied_outlook = implied_outlooks_list[row][day]
        if max_implied_outlook_by_day[day] < implied_outlook: max_implied_outlook_by_day[day] = implied_outlook

max_implied_outlook_by_day = np.expand_dims(max_implied_outlook_by_day, axis=0)

##########################
# Data for Third Subplot #
##########################
min_outlook_by_day = np.array([])

for day in range(len(outlooks_list[0])):
    min_outlook_by_day = np.append(min_outlook_by_day,7)
    for row in range(year_end-year_start+1):
        outlook = outlooks_list[row][day]
        #if day==58 or day==59: print row,day,outlook
        #if (day==58 or day==59) and min_outlook_by_day[day] > outlook and outlook > -1: print min_outlook_by_day[day],"Will change min"
        if min_outlook_by_day[day] > outlook and outlook > -1: min_outlook_by_day[day] = outlook

min_outlook_by_day = np.expand_dims(min_outlook_by_day, axis=0)

#for thing in range(58,60):
#    print outlooks_list[0][thing],outlooks_list[1][thing],outlooks_list[2][thing],outlooks_list[3][thing],outlooks_list[4][thing],outlooks_list[5][thing],outlooks_list[6][thing],outlooks_list[7][thing],outlooks_list[8][thing],outlooks_list[9][thing],outlooks_list[10][thing],outlooks_list[11][thing],outlooks_list[12][thing],outlooks_list[13][thing],outlooks_list[14][thing],outlooks_list[15][thing]
#    print min(outlooks_list[0][thing],outlooks_list[1][thing],outlooks_list[2][thing],outlooks_list[3][thing],outlooks_list[4][thing],outlooks_list[5][thing],outlooks_list[6][thing],outlooks_list[7][thing],outlooks_list[8][thing],outlooks_list[9][thing],outlooks_list[10][thing],outlooks_list[11][thing],outlooks_list[12][thing],outlooks_list[13][thing],outlooks_list[14][thing],outlooks_list[15][thing])
#    print min_outlook_by_day[0][thing]

# Repeat, but for implied_outlooks
min_implied_outlook_by_day = np.array([])

for day in range(len(implied_outlooks_list[0])):
    min_implied_outlook_by_day = np.append(min_implied_outlook_by_day,7)
    for row in range(year_end-year_start+1):
        implied_outlook = implied_outlooks_list[row][day]
        if min_implied_outlook_by_day[day] > implied_outlook and implied_outlook > -1: min_implied_outlook_by_day[day] = implied_outlook

min_implied_outlook_by_day = np.expand_dims(min_implied_outlook_by_day, axis=0)

###########################
# Data for Fourth Subplot #
###########################
# Make lists
sum_tstmlove = []     # I want TSTM to mean something!
sum_counter = []

x = range(len(outlooks_list))
smooth_tstmlove = []

# Add cumulative risk from each day.
for i,year in enumerate(outlooks_list):
    for j,day in enumerate(year):
        if day == -1: tstmlove=0.0; counter=0;  #Missing
        if day == 0:  tstmlove=0.0; counter=1;  #None
        if day == 1:  tstmlove=1.0; counter=1;  #TSTM
        if day == 2:  tstmlove=2.0; counter=1;  #MRGL
        if day == 3:  tstmlove=3.0; counter=1;  #SLGT
        if day == 4:  tstmlove=4.0; counter=1;  #ENH
        if day == 5:  tstmlove=5.0; counter=1;  #MDT
        if day == 6:  tstmlove=6.0; counter=1;  #HIGH

        if i==0:
            sum_tstmlove.append(tstmlove)
            sum_counter.append(counter)
        else:
            sum_tstmlove[j] = sum_tstmlove[j] + tstmlove
            sum_counter[j] = sum_counter[j] + counter

# Normalize the plots by the amount of data gathered on that calendar day.
for i,day in enumerate(sum_counter):
    if sum_counter[i]>0:
        sum_tstmlove[i]=sum_tstmlove[i]/sum_counter[i]

# Smooth the data
for i,day in enumerate(sum_tstmlove):
    s2s_tstmlove = 0
    for j in range(-4,5):
        if i+j>365:
            s2s_tstmlove = s2s_tstmlove + sum_tstmlove[i+j-366]

        else:
            s2s_tstmlove = s2s_tstmlove + sum_tstmlove[i+j]
    smooth_tstmlove.append(s2s_tstmlove / 9)

y4_max = max(smooth_tstmlove)

### Repeat, but for implied_outlooks
# Make lists
sum_implied_tstmlove = []     # I want TSTM to mean something!
sum_implied_counter = []

x = range(len(implied_outlooks_list[0]))
smooth_implied_tstmlove = []

# Add cumulative risk from each day.
for i,year in enumerate(implied_outlooks_list):
    for j,day in enumerate(year):
        if day == -1: implied_tstmlove=0.0; implied_counter=0;  #Missing
        if day == 0:  implied_tstmlove=0.0; implied_counter=1;  #None
        if day == 1:  implied_tstmlove=1.0; implied_counter=1;  #TSTM
        if day == 2:  implied_tstmlove=2.0; implied_counter=1;  #MRGL
        if day == 3:  implied_tstmlove=3.0; implied_counter=1;  #SLGT
        if day == 4:  implied_tstmlove=4.0; implied_counter=1;  #ENH
        if day == 5:  implied_tstmlove=5.0; implied_counter=1;  #MDT
        if day == 6:  implied_tstmlove=6.0; implied_counter=1;  #HIGH

        if i==0:
            sum_implied_tstmlove.append(implied_tstmlove)
            sum_implied_counter.append(implied_counter)
        else:
            sum_implied_tstmlove[j] = sum_implied_tstmlove[j] + implied_tstmlove
            sum_implied_counter[j] = sum_implied_counter[j] + implied_counter

# Normalize the plots by the amount of data gathered on that calendar day.
for i,day in enumerate(sum_implied_counter):
    if sum_implied_counter[i]>0:
        sum_implied_tstmlove[i]=sum_implied_tstmlove[i]/sum_implied_counter[i]

# Smooth the data
for i,day in enumerate(sum_implied_tstmlove):
    s2s_implied_tstmlove = 0
    for j in range(-4,5):
        if i+j>365:
            s2s_implied_tstmlove = s2s_implied_tstmlove + sum_implied_tstmlove[i+j-366]

        else:
            s2s_implied_tstmlove = s2s_implied_tstmlove + sum_implied_tstmlove[i+j]
    smooth_implied_tstmlove.append(s2s_implied_tstmlove / 9)


####################
# Make the Figures #
####################
#Make a figure with subplots (4 figs in a column)
fig, axes = plt.subplots(4,1,sharex="True",figsize=(12,9))
ax1, ax2, ax3, ax4 = [axes[0],axes[1],axes[2],axes[3]]

#Make the figures different sizes
gs = gridspec.GridSpec(4,1,height_ratios=[14,1,1,6])
gs.update(hspace=0.04)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
ax3 = plt.subplot(gs[2])
ax4 = plt.subplot(gs[3])

#Set the y-axis ranges for each subplot
ax1.set_ylim(0,14)
ax2.set_ylim(0,1)
ax3.set_ylim(0,4)
ax4.set_ylim(0,4)

#Set the x-axis ranges for all subplots
ax1.set_xlim(-1,366)
ax2.set_xlim(-1,366)
ax3.set_xlim(-1,366)
ax4.set_xlim(-1,366)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)

#Set the x-axis labels
#months_days = [31,60,91,121,152,182,213,244,274,305,335]
months_days = [30,59,90.5,120.5,151.5,181.5,212.5,243.5,273.5,304.5,334.5]
x_ticks =     [15.0,44.5,74.5,105,135.5,166,196.5,227.5,258,288.5,319,349.5]
months_lables = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
plt.vlines(months_days,0,14,colors='k')
plt.xticks(x_ticks,months_lables)

#Remove axis ticks
ax1.tick_params(axis='both',length=0)
ax2.tick_params(axis='both',length=0)
ax3.tick_params(axis='both',length=0)
ax4.tick_params(axis='both',length=0)

#################
# First Subplot #
#################
#Set the y-axis labels for the first subplot
plt.sca(ax1)
y1_labels = ['2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018']
y1_ticks = range(len(y1_labels))
plt.yticks(y1_ticks,y1_labels)
ax1.set_ylim(y1_ticks[0]-0.5,y1_ticks[-1]+0.5)
plt.gca().invert_yaxis()

#Set the vertical lines for the first subplot
plt.vlines(months_days,y1_ticks[0]-0.5,y1_ticks[-1]+0.5,colors='k')

##################
# Second Subplot #
##################
#Set the y-axis label for the second subplot
plt.sca(ax2)
y2_labels = ['Max']
y2_ticks = [0]
plt.yticks(y2_ticks,y2_labels)
ax2.set_ylim(y2_ticks[0]-0.5,y2_ticks[-1]+0.5)

#Set the vertial lines for the second subplot
plt.vlines(months_days,y2_ticks[0]-0.5,y2_ticks[-1]+0.5,colors='k')

#################
# Third Subplot #
#################
#Set the y-axis label for the second subplot
plt.sca(ax3)
y3_labels = ['Min']
y3_ticks = [0]
plt.yticks(y3_ticks,y3_labels)
ax3.set_ylim(y3_ticks[0]-0.5,y3_ticks[-1]+0.5)

#Set the vertial lines for the second subplot
plt.vlines(months_days,y3_ticks[0]-0.5,y3_ticks[-1]+0.5,colors='k')

##################
# Fourth Subplot #
##################
#Set the y-axis label for the third subplot
plt.sca(ax4)
y4_labels = ['Lower Category','Higher Category']
y4_ticks = [0.3,y4_max]
plt.yticks(y4_ticks,y4_labels)
ax4.set_ylim(0,y4_max+0.5)

############################################
# Make colors for a colorbar and plot data #
############################################
# SPC Colors
# Gray (Missing Data), White (No TSTM), TSTM, MRGL, SLGT, ENH, MDT, HIGH
colors = [(0.784,0.784,0.784),(1,1,1),(0.753,0.910,0.753),(0.498,0.773,0.498),(0.965,0.965,0.498),(0.902,0.761,0.498),(0.902,0.498,0.498),(1,.498,1)]
levels = [-1, 0, 1, 2, 3, 4, 5, 6, 7]
cm, norm = from_levels_and_colors(levels, colors)

# Put data on the first subplot
cax = ax1.imshow(outlooks_list, aspect='auto', interpolation='nearest', cmap=cm, norm=norm)

# Put data on the second subplot
ax2.imshow(max_outlook_by_day, aspect='auto', interpolation='nearest', cmap=cm, norm=norm)

# Put data on the third subplot
ax3.imshow(min_outlook_by_day, aspect='auto', interpolation='nearest', cmap=cm, norm=norm)

# Put data on the fourth subplot
ax4.plot(x,smooth_tstmlove,'b-',label="TSTM=1, HIGH=6, Smoothed")

# Colorbar
cbar = plt.colorbar(cax,ticks=[-1,0,1,2,3,4,5,6],orientation='horizontal',spacing='uniform',pad=0.25,fraction=0.21,use_gridspec=False,anchor=(0.5,0.1))
cbar.ax.get_xaxis().set_ticks([])
tick_labels = ['Missing','None','TSTM','MRGL','SLGT','ENH','MDT','HIGH']
for j, label in enumerate(tick_labels):
        cbar.ax.text((2.0*j+1.0)/16, .5, label, ha='center', va='center')

#Make a title for the whole plot.
ax1.set_title('SPC Outlooks History')

#Add a signature
plt.text(366,-(y4_max+0.5)*0.375,'Stephen Mullens @srmullens',fontsize=10,horizontalalignment='right')

#Save
plt.savefig('figure.png')





### Attempt to make a second plot ###


####################
# Make the Figures #
####################
#Make a figure with subplots (3 figs in a column)
fig, axes = plt.subplots(4,1,sharex=True,figsize=(12,9))
ax1, ax2, ax3, ax4 = [axes[0],axes[1],axes[2],axes[3]]

#Make the figures different sizes
gs = gridspec.GridSpec(4,1,height_ratios=[14,1,1,6])
gs.update(hspace=0.04)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
ax3 = plt.subplot(gs[2])
ax4 = plt.subplot(gs[3])

#Set the y-axis ranges for each subplot
ax1.set_ylim(0,14)
ax2.set_ylim(0,1)
ax2.set_ylim(0,1)
ax4.set_ylim(0,4)

#Set the x-axis ranges for all subplots
ax1.set_xlim(-1,366)
ax2.set_xlim(-1,366)
ax3.set_xlim(-1,366)
ax4.set_xlim(-1,366)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)

#Set the x-axis labels
#months_days = [31,60,91,121,152,182,213,244,274,305,335]
months_days = [30,59,90.5,120.5,151.5,181.5,212.5,243.5,273.5,304.5,334.5]
x_ticks =     [15.0,44.5,74.5,105,135.5,166,196.5,227.5,258,288.5,319,349.5]
months_lables = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
plt.vlines(months_days,0,14,colors='k')
plt.xticks(x_ticks,months_lables)

#Remove axis ticks
ax1.tick_params(axis='both',length=0)
ax2.tick_params(axis='both',length=0)
ax3.tick_params(axis='both',length=0)
ax4.tick_params(axis='both',length=0)

#################
# First Subplot #
#################
#Set the y-axis labels for the first subplot
plt.sca(ax1)
y1_labels = ['2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018']
y1_ticks = range(len(y1_labels))
plt.yticks(y1_ticks,y1_labels)
ax1.set_ylim(y1_ticks[0]-0.5,y1_ticks[-1]+0.5)
plt.gca().invert_yaxis()

#Set the vertical lines for the first subplot
plt.vlines(months_days,y1_ticks[0]-0.5,y1_ticks[-1]+0.5,colors='k')

##################
# Second Subplot #
##################
#Set the y-axis label for the second subplot
plt.sca(ax2)
y2_labels = ['Max']
y2_ticks = [0]
plt.yticks(y2_ticks,y2_labels)
ax2.set_ylim(y2_ticks[0]-0.5,y2_ticks[-1]+0.5)

#Set the vertial lines for the second subplot
plt.vlines(months_days,y2_ticks[0]-0.5,y2_ticks[-1]+0.5,colors='k')

#################
# Third Subplot #
#################
#Set the y-axis label for the second subplot
plt.sca(ax3)
y3_labels = ['Min']
y3_ticks = [0]
plt.yticks(y3_ticks,y3_labels)
ax3.set_ylim(y3_ticks[0]-0.5,y3_ticks[-1]+0.5)

#Set the vertial lines for the second subplot
plt.vlines(months_days,y3_ticks[0]-0.5,y3_ticks[-1]+0.5,colors='k')

##################
# Fourth Subplot #
##################
#Set the y-axis label for the third subplot
plt.sca(ax4)
y4_labels = ['Lower Category','Higher Category']
y4_ticks = [0.3,y4_max]
plt.yticks(y4_ticks,y4_labels)
ax4.set_ylim(0,y4_max+0.5)

############################################
# Make colors for a colorbar and plot data #
############################################
# SPC Colors
# Gray (Missing Data), White (No TSTM), TSTM, MRGL, SLGT, ENH, MDT, HIGH
colors = [(0.784,0.784,0.784),(1,1,1),(0.753,0.910,0.753),(0.498,0.773,0.498),(0.965,0.965,0.498),(0.902,0.761,0.498),(0.902,0.498,0.498),(1,.498,1)]
levels = [-1, 0, 1, 2, 3, 4, 5, 6, 7]
cm, norm = from_levels_and_colors(levels, colors)

# Put data on the first subplot
cax = ax1.imshow(implied_outlooks_list, aspect='auto', interpolation='nearest', cmap=cm, norm=norm)

# Put data on the second subplot
ax2.imshow(max_implied_outlook_by_day, aspect='auto', interpolation='nearest', cmap=cm, norm=norm)

# Put data on the second subplot
ax3.imshow(min_implied_outlook_by_day, aspect='auto', interpolation='nearest', cmap=cm, norm=norm)

# Put data on the third subplot
ax4.plot(x,smooth_implied_tstmlove,'b-',label="TSTM=1, HIGH=6, Smoothed")

# Colorbar
cbar = plt.colorbar(cax,ticks=[-1,0,1,2,3,4,5,6],orientation='horizontal',spacing='uniform',pad=0.25,fraction=0.21,use_gridspec=False,anchor=(0.5,0.1))
cbar.ax.get_xaxis().set_ticks([])
tick_labels = ['Missing','None','TSTM','MRGL','SLGT','ENH','MDT','HIGH']
for j, label in enumerate(tick_labels):
        cbar.ax.text((2.0*j+1.0)/16, .5, label, ha='center', va='center')

#Make a title for the whole plot.
ax1.set_title('SPC Implied Outlooks History')

#Add a signature
plt.text(366,-(y4_max+0.5)*0.375,'Stephen Mullens @srmullens',fontsize=10,horizontalalignment='right')

#Save
plt.savefig('figure_implied.png')
