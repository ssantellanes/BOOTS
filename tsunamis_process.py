#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 11:34:25 2020

@author: seansantellanes
"""
import pandas as pd

def main():
    df=pd.read_csv('/Users/seansantellanes/Documents/Background_Spectra/tsunamis-2024-03-14_13-43-06_-0500.csv')
    df=df.dropna()
    df=df.applymap(str)

    df=df.drop(['Hr','Mn','Sec'],axis=1)
    df=df.rename(columns={'Year':'year','Mo':'month','Dy':'day'})
    df=pd.to_datetime(df,errors='ignore')
    df.to_csv('tsunami_database_pacific_2024.csv',index=False,header=['Tsunamis'])
    print(df)
    
main()