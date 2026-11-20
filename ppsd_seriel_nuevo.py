#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 15:14:13 2020

@author: seansantellanes
"""

from pandas import read_hdf
from numpy import array,arange
import numpy as np

from matplotlib import pyplot as plt
import scipy.integrate as integrate

import datetime
import pandas as pd
import warnings
warnings.simplefilter('ignore',UserWarning)
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from multitaper import MTSpec
from obspy.signal import PPSD
from obspy import Stream,Trace
from mudpy.forward import highpass
from mudpy.forward import lowpass
import logging
from random import gauss
plt.close("all")
def daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days)):
        yield start_date + datetime.timedelta(n)
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

startTime = datetime.datetime.now()
meteotsunami_count=pd.read_csv("/Users/seansantellanes/Documents/infragravity_waves/water_level_codes/meteotsunami_count.csv")
count=meteotsunami_count.iloc[:,1].tolist()
sta_names=meteotsunami_count.iloc[:,0].tolist()
for station in sta_names:
    logging.basicConfig(filename=f'LOG/ppsd_{station}.log',encoding='utf-8',level=logging.INFO,filemode='w',force=True)
    fcorner=[1./(3600*2),1./(32)]
    sta_info=pd.read_csv("/Users/seansantellanes/Documents/infragravity_waves/water_level_codes/regions_sep_num_2026.csv")
    start=[]
    end=[]
    for i in range(len(sta_info.iloc[:,8])):
        start.append(sta_info.iloc[i,8][:10])
        end.append(sta_info.iloc[i,8][11:])
    sta_info["Start"]=pd.Series(start)
    sta_info["End"]=pd.Series(end)
    sta_info["Start"]=pd.to_datetime(sta_info["Start"])
    sta_info["End"]=pd.to_datetime(sta_info["End"])
    idx=list(np.where(sta_info["Station ID"]==int(station))[0])
    start_date=sta_info["Start"].iloc[idx].dt.date.item()
    end_date=sta_info["End"].iloc[idx].dt.date.item()
    print(start_date)
    #start_date = datetime.date(2016, 1, 1)
    infragravity_wave_list=[]
    day_list=[]
    #end_date = datetime.date(2024, 5, 31)

    counter=0
    tsu_counter=0


    weekStart=start_date
    df=pd.read_hdf(f"/Users/seansantellanes/Documents/infragravity_waves/PTHA_codes/{station}/pressure_res_long.h5",parse_dates=['1'],index_col=['1'])
    df["datetime"] = pd.to_datetime(df.iloc[:,1],unit="s")
    df=df.set_index("datetime")
    flag=0
    ppsd_init=0
    while start_date < end_date:
        dayItr=datetime.timedelta(days=1)
        nextday=start_date+dayItr
        dayStartStr=str(start_date)
        nextDayStr=str(nextday)

        
        try:
            
            #df2=df[dayStartStr:nextDayStr]
            df2=df.loc[dayStartStr]
        except KeyError:
            start_date=start_date+dayItr
            continue
        #print('Made it through slice.')
        time_epi = np.datetime64(weekStart)
        #Check for empty dataframe
        if df2.empty == True:
            #print('This beh emteh!!!')
            start_date=start_date+dayItr
            continue
        #Check for null values in the dataframe
        if df2.isnull().values.any() == False:
            pass
        #Move on to next dataframe if this one has nulls
        else:
            #print('Oh, crap!')
            dt=0 
            start_date=start_date+dayItr
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
                #print(seconds)
            elif len(dtset) > 1:
                dt = min(dtset)
                seconds = dt / np.timedelta64(1, 's')
                start_date=start_date+dayItr
                continue
        except IndexError as e:
            print(e)
            print('Ope!')
            dt=0 
            start_date=start_date+dayItr
            continue
        dt=int(seconds) #seconds
        pressure=array(df2.iloc[:,0])
        #Check to see if Hi-Res before proceeding
        # if start_date in tsunami:
        #     #print(f'Tsunami at {start_date}, skipping 14 days')
        #     dayItr=datetime.timedelta(days=13)
        #     start_date=start_date+dayItr
        #     continue
        if dt > 15:
            print('Not Hi-Res')
            start_date=start_date+dayItr
            continue
        else:
            flag+=1
            #print('Hi-Res. Proceeding...')
        #slice out fall season
        
        time=arange(0,len(pressure)*dt,dt)
        t=np.zeros(len(time))
        if len(dtset) > 1:
            #print('Hope this interpolates')
            for k in range(len(time)):
                t[k]=((df2.index.values[k])-time_epi)/np.timedelta64(1, 's')
            #print(t[0], t[-1])
            t_in=arange(t[0],t[-1],15)
            f=interp1d(t,pressure)
            data_in=f(t_in)
            #data_smooth=smooth(data_in,30)
            pressure=data_in-np.mean(data_in)
            dt=t_in[2]-t_in[1]
        else:
            pass
        '''
        try:
            new_data_try=pd.read_csv(f"/Users/seansantellanes/Documents/infragravity_waves/water_level_codes/outputs/csv/{station}_meteotsunami.csv")
            ms_list=new_data_try.iloc[:,0].to_list()
            
            if (start_date.month >=1 and start_date.month < 5) or (start_date.month >=10 and start_date.month <= 12):
                counter+=1
                print(f"Could be a meteotsunami here!\n")
                print(f"There have been {counter} meteotsunamis!")
                dayItr=datetime.timedelta(days=4)
                nextday=start_date+dayItr
                start_date=nextday
                continue
            else:
                
                pass
        except pd.errors.EmptyDataError:
            start_date=nextday
            pass
        except FileNotFoundError:
            start_date=nextday
            pass
        '''
        new_data_tsu=pd.read_csv(f"/Users/seansantellanes/Documents/infragravity_waves/water_level_codes/tsunami_database_pacific.csv")
        mt_list=new_data_tsu.iloc[:,0].to_list()
        
        if dayStartStr in mt_list:
            tsu_counter+=1
            print(f"Tsunami here at {dayStartStr}!\n")
            print(f"There have been {tsu_counter} tsunamis!")
            dayItr=datetime.timedelta(days=14)
            nextday=start_date+dayItr
            start_date=nextday
            continue
        else:
            
            pass
        if dayStartStr in mt_list:
            tsu_counter+=1
            print(f"Tsunami here at {dayStartStr}!\n")
            print(f"There have been {tsu_counter} tsunamis!")
            dayItr=datetime.timedelta(days=14)
            nextday=start_date+dayItr
            start_date=nextday
            continue
        else:
            
            pass
        #Apply corection factors, first convert from dBar to Pa
        pressure*=1e4
        
        #nowd ivide by rho*g with rho=1037kg/m^3 and g=9.8m/s/s 
        sea_surface_height = pressure / (1037*9.8)
        
        #demean before spectra
        sea_surface_height -= sea_surface_height.mean()
        P1= MTSpec(
            sea_surface_height,nw=3.5,
            kspec=5, dt=dt)
        f,spec=P1.rspec()
        
        file="regions_sep_num.csv"
        dr=pd.read_csv(file)
        dr=dr.set_index("Station ID")
        depth=dr.loc[int(station),"Depth [m]"]
        
        k=(2*np.pi*f)/np.sqrt(9.81*depth)
        M=(1037*9.81)/np.cosh(k*depth)
        #M=(1037*9.81)
        E=M**2*spec
        '''
        try:
            x_1=f[144]
            x_2=f[1977]
        except IndexError:
            print("Did I get IndexError?")
            start_date=start_date+dayItr
            continue
        
        dx=(f[1977]-f[144])/len(f[144:1977])
        #H_integrate=integrate.quad(E,8.3e-4,1.1e-2)
        H_integrate=np.trapezoid(E[144:1977],f[144],f[1977])
        H_ig=4*np.sqrt(H_integrate.sum())
        print(E[144:1977])
        if H_ig > 0.:
            infragravity_wave_list.append(H_ig)
            day_list.append(start_date)
        '''
        st=Stream(Trace())
        #print(E)
        st[0].data=sea_surface_height
        #st[0].data=E[:,0]
        st[0].stats.delta=15
        st[0].stats.starttime=start_date
        #print(st)
        '''
        if flag==10:
            import matplotlib.pyplot as plt
            import sys
            plt.plot(spec[3:])
            plt.savefig("debug.png")
            sys.exit()
        #st[0].stats.endtime=weekStart+datetime.timedelta(days=14)
        '''
        paz = {'gain': 1,'poles': [1],'sensitivity': 1,'zeros': [0j, 0j]}
        #print(ppsd_init)
        if ppsd_init==0:
            ppsd=PPSD(st[0].stats,paz,db_bins=(-80, 0, 1.0),period_limits=[30,3600*2],ppsd_length=3600*24,overlap=0.0,special_handling='ringlaser')
            ppsd_init+=1
        elif ppsd_init >= 1:
            if ppsd_init%10==0:
                #print(start_date)
                #series=np.array([gauss(0,0.0001) for i in range(len(sea_surface_height))])
                #st[0].data+=series
                ppsd.add(st)
                ppsd_init+=1
            else:
                ppsd.add(st)
                ppsd_init+=1
        if len(ppsd.psd_values) > 0:
            #print("Accessing ppsd values")
            
            min_psd=np.min(ppsd.psd_values[0])
            if min_psd < 40.:
                logging.info(f'{start_date} idx {len(ppsd.psd_values)-1} at flag={flag} in bottom branch.')

                if str(start_date) == '2021-12-24':
                    '''
                    ssh_fil=lowpass(sea_surface_height,fcorner,1/15,2,zerophase=True)
                    print(np.max(ssh_fil))
                    plt.plot(ssh_fil,label=f'{start_date} low branch')
                    #plt.xlim([2000,2500])
                    plt.legend()
                    plt.savefig(f'time_series_png/{start_date}_{station}_bottom_branch.pdf')
                    '''
                else:
                    pass
            elif min_psd >= 40.:
                logging.info(f'{start_date} idx {len(ppsd.psd_values)-1} at flag={flag} in top branch.')
            #    pass
                if str(start_date) == '2021-12-24':
                    print(len(ppsd.psd_values)-1)
                    '''
                    #ssh_fil=lowpass(sea_surface_height,fcorner,1/15,2,zerophase=True)
                    ssh_fil=sea_surface_height
                    print(np.max(ssh_fil))
                    plt.plot(ssh_fil, label=f'{start_date} top branch')
                    #plt.xlim([2000,2500])
                    plt.legend()
                    plt.savefig(f'time_series_png/{start_date}_{station}_top_branch.pdf')
                    '''
                else:
        
                    pass
        
        else:
            pass
    
        start_date=start_date+dayItr
        

    # plt.legend() 
    # plt.savefig(f'time_series_png/{start_date}.png')   
    if tsu_counter > counter:
        print("Tsunamis have won!")
    elif counter > tsu_counter:
        print("Meteotsunamis have won!")
    else:
        print("It's a draw!")
    infragravity_wave_arr=np.array(infragravity_wave_list)
    days_ig=pd.to_datetime(day_list)
    df_dict={'Time':days_ig,'infragravity heights':infragravity_wave_arr}
    df=pd.DataFrame(data=df_dict)
    df.to_csv(f'/Users/seansantellanes/Documents/infragravity_waves/water_level_codes/ig_uncorrected/{station}_ig_uncorrected.csv')
    ppsd.save_npz(f'/Users/seansantellanes/Documents/infragravity_waves/water_level_codes/NPZ_w_meteo/{station}_psd.npz')
    #set legend, size,...
print(datetime.datetime.now() - startTime)
