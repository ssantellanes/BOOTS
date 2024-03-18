#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 15:14:13 2020

@author: seansantellanes
"""

from pandas import read_csv
from numpy import array,arange
import numpy as np
from multitaper import MTSpec
from matplotlib import pyplot as plt
from scipy.stats import linregress
import os,shutil
import sys,datetime
import pandas as pd
import os.path
from scipy.interpolate import interp1d
import plotly.graph_objects as go
#from mudpy.forward import lowpass
region='Cascadia'
station='46407'
date='2010_2022'
depth=5812
regionDir=region+'/csv/'
path1=regionDir+station+'_'+date+'_'
path2=region+'/'+station+'/'+date+'/biweekly/figures/'
filename=path1+'biweekly_slope_intercept.csv'
x_array=[]
y_array=[]
z_array=[]
if os.path.exists(filename):
    os.remove(filename)
else:
    print('cannot delete file as it does not exist')
if os.path.exists(path2):
    shutil.rmtree(path2)
    
else:
    print('This directory did not exist!')
os.makedirs(path2)
if os.path.exists(regionDir):
    pass
else:
    print('Making directory:',regionDir)
    os.makedirs(regionDir)
startTime = datetime.datetime.now()
#weekStart=datetime.datetime.strptime(str(sys.argv[1]),'%Y-%m-%d').date()
weekStart=datetime.datetime.strptime('2010-01-01','%Y-%m-%d').date()
df=read_csv(station+'/pressure_res_long.csv',parse_dates=['Time'],index_col=['Time'])
for i in range(1, 365*13):
    weekItr=datetime.timedelta(days=1)
    nextWeek=weekStart+weekItr
    weekStartStr=str(weekStart)
    nextWeekStr=str(nextWeek)
    print(weekStart)
    print((weekStartStr), (nextWeekStr))
    try:
        
        df2=df[weekStartStr:nextWeekStr]
    # except AssertionError:
    #     weekStart=weekStart+datetime.timedelta(days=14)
    #     print('There was an AssertionError!')
    #     mls=np.nan
    #     bls=np.nan
    #     slope_inter={'Time':[weekStartStr],'MLS':[mls],'BLS':[bls]}
    #     df3=pd.DataFrame(slope_inter,columns=['Time','MLS','BLS'])
    #     if not os.path.isfile(filename):
    #         df3.to_csv(filename, header='column_names')
    #     else: # else it exists so append without writing the header
    #         df3.to_csv(filename, mode='a', header=False)
    #     continue
    # print('Made it through slice.')
    except KeyError:
        weekStart=weekStart+datetime.timedelta(days=1)
        print('There was an AssertionError!')
        mls=np.nan
        bls=np.nan
        slope_inter={'Time':[weekStartStr],'MLS':[mls],'BLS':[bls]}
        df3=pd.DataFrame(slope_inter,columns=['Time','MLS','BLS'])
        if not os.path.isfile(filename):
            df3.to_csv(filename, header='column_names')
        else: # else it exists so append without writing the header
            df3.to_csv(filename, mode='a', header=False)
        continue
    time_epi = np.datetime64(weekStart)
    #Check for empty dataframe
    if df2.empty == True:
        print('This beh emteh!!!')
        weekStart=weekStart+datetime.timedelta(days=1)
        mls=np.nan
        bls=np.nan
        slope_inter={'Time':[weekStartStr],'MLS':[mls],'BLS':[bls]}
        df3=pd.DataFrame(slope_inter,columns=['Time','MLS','BLS'])
        if not os.path.isfile(filename):
            df3.to_csv(filename, header='column_names')
        else: # else it exists so append without writing the header
            df3.to_csv(filename, mode='a', header=False)
        continue
    #Check for null values in the dataframe
    if df2.isnull().values.any() == False:
        pass
    #Move on to next dataframe if this one has nulls
    else:
        print('Oh, crap!')
        dt=0 
        weekStart=weekStart+datetime.timedelta(days=1)
        df2=df2[pd.notnull(df2['Pressure [dbar]'])]
        mls=np.nan
        bls=np.nan
        slope_inter={'Time':[weekStartStr],'MLS':[mls],'BLS':[bls]}
        df3=pd.DataFrame(slope_inter,columns=['Time','MLS','BLS'])
        if not os.path.isfile(filename):
            df3.to_csv(filename, header='column_names')
        else: # else it exists so append without writing the header
            df3.to_csv(filename, mode='a', header=False)
        continue
    #Check for any sets of mismatching dts in the two week period
    try:
        dtls = []
        for j in range(1,df2.index.size):
            dt = df2.index.values[j]-df2.index.values[j-1]
            dtls.append(dt)
        dtset = set(dtls)
        if len(dtset) == 1:
            seconds = dt / np.timedelta64(1, 's')
            print(seconds)
        elif len(dtset) > 1:
            dt = min(dtset)
            seconds = dt / np.timedelta64(1, 's')
            weekStart=weekStart+datetime.timedelta(days=14)
            print(weekStart,i,' Welp looks like it was too big.')
            continue
    except IndexError as e:
        print(e)
        print('Ope!')
        dt=0 
        weekStart=weekStart+datetime.timedelta(days=1)
        print(weekStart,i)
        continue
    dt=int(seconds) #seconds
    pressure=array(df2['Pressure [dbar]'])
    #Check to see if Hi-Res before proceeding
    if dt > 15:
        print('Not Hi-Res')
        weekStart=weekStart+datetime.timedelta(days=1)
        print(weekStart,i)
        mls=np.nan
        bls=np.nan
        slope_inter={'Time':[weekStartStr],'MLS':[mls],'BLS':[bls]}
        df3=pd.DataFrame(slope_inter,columns=['Time','MLS','BLS'])
        if not os.path.isfile(filename):
            df3.to_csv(filename, header='column_names')
        else: # else it exists so append without writing the header
            df3.to_csv(filename, mode='a', header=False)
        continue
    else:
        print('Hi-Res. Proceeding...')
    #slice out fall season
    
    time=arange(0,len(pressure)*dt,dt)
    t=np.zeros(len(time))
    if len(dtset) > 1:
        print('Hope this interpolates')
        for k in range(len(time)):
            t[k]=((df2.index.values[k])-time_epi)/np.timedelta64(1, 's')
        print(t[0], t[-1])
        t_in=arange(t[0],t[-1],15)
        f=interp1d(t,pressure)
        data_in=f(t_in)
        pressure=data_in-np.mean(data_in)
        dt=t_in[2]-t_in[1]
    else:
        pass
    #Apply corection factors, first convert from dBar to Pa
    pressure*=1e4
    
    #nowd ivide by rho*g with rho=1037kg/m^3 and g=9.8m/s/s 
    sea_surface_height = pressure / (1037*9.8)
    
    #demean before spectra
    sea_surface_height -= sea_surface_height.mean()
   
    #get unfiltered spectra
    # psd1= mtspec(
    #     data=sea_surface_height, delta=dt, time_bandwidth=3.5,
    #     number_of_tapers=5, nfft=len(pressure), statistics=True)
    psd= MTSpec(
        sea_surface_height,nw=3.5,
        kspec=5, dt=dt)
    f=psd.freq
    psd1=psd.spec
    psd1=psd1*(10000./60)
    period=(1/f)/60
    #least squares inversion
    period=period[1:]
    #new line of code
    #Converting to wavenumber to check something
    wavenumber=1./(period*np.sqrt(9.8*3772))
    # print((wavenumber))
    # xl=(1./(12*np.sqrt(9.8*3772)))
    # xr=(1./(100*np.sqrt(9.8*3772)))
    y=psd1[1:]
    x=np.asarray(period)
    y=np.asarray(y)
    x_array.append(x)
    z_array.append(y)
    y_array.append(str(weekStart))
    if dt > 15:
        
        continue
    else:
        mask=(x>12) & (x<100)
        # mask=(x>xr) & (x<xl)
    xls=x[mask]
    yls=y[mask]
    print(len(xls),len(yls))
    # if len(xls) == 0 and len(yls) ==0:
    #     weekStart=weekStart+datetime.timedelta(days=14)
    #     continue
    # else:
    #     pass
    stats=linregress(np.log10(xls),np.log10(yls))
    #sys.exit()
    mls=stats.slope
    bls=stats.intercept
    ys=np.power(10,bls)*xls**mls
    print('MLS: ' + str(mls))
    print('BLS: ' + str(bls))
    mls_legend = round(mls, 5)
    bls_legend = round(bls, 5)
    #convert fromf requency to period in minutes
    
    
    #get rid of the first elemnt which is == inf
    period=period[1:]
    psd1=psd1[1:]
    
    #Set up plot path
    plotOut=path2+'week_'+str(weekStart).replace('-','_')+'.png'

    #make plots
#    plt.subplot(121)
#    time = time/3600 #to hours for convenience in plotting
#    plt.plot(time,sea_surface_height,label='unfitlered')
#    
#    plt.legend()
#    plt.xlabel('Time (hours)')
#    plt.ylabel('Sea surface height (m)')
    
#    plt.subplot(122)
    if i%20==0:
        plt.figure()
        plt.axvspan(5,120,color='green',alpha=0.25)
        plt.loglog(x,y)
        plt.plot(xls,ys,linewidth=3,label='Least-squares model:\nslope = ' + str(mls_legend) + ',\n intercept = ' + str(bls_legend))
        # # y1=x**2
        # # y2=x**1.8
        # # y3=x**2.2
        #plt.figure()
        # plt.loglog(x,y1,label='n=2')
        # plt.loglog(x,y2,linestyle='dashed',label='n=1.8')
        # plt.loglog(x,y3,linestyle='dashed',label='n=2.2')
    
        # plt.xlabel('Period (min)')
        plt.xlabel('Period (min)')
        plt.ylabel(r'$cm^{2} cpm^{-1}$')
        plt.xlim([1,1000])
        plt.ylim([10e-5,10e2])
        
        
        plt.grid()
        plt.legend(loc = 4)
        plt.savefig(plotOut)
    weekStart=weekStart+datetime.timedelta(days=1)
    print(weekStart,plotOut,i)
    slope_inter={'Time':[weekStartStr],'MLS':[mls],'BLS':[bls]}
    df3=pd.DataFrame(slope_inter,columns=['Time','MLS','BLS'])
    if not os.path.isfile(filename):
        df3.to_csv(filename, header='column_names')
    else: # else it exists so append without writing the header
        df3.to_csv(filename, mode='a', header=False)
#set legend, size,...
print(datetime.datetime.now() - startTime)