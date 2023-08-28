#Optimize scan pattern using the LiDAR Barnes Statistical Objective Analysis tool [Letizia et al., Atmos Meas Tech, 2021a,b; Letizia et al., WAKE conf., 2023]
#The copyright of this software belongs to Stefano Letizia from NREL. 
#The author's contribution must be acknoledged in all the pubblications that include results form this software.
#The author is not responsible for any erroneous and/or financially or scientific harmful results produced by the present software.
#08/28/2023: created, finalized

import os
cd=os.path.dirname(__file__)
import sys
sys.path.append('functions')
import LiSBOA_functions as LiS
import pandas as pd
import numpy as np
import warnings
import xarray as xr
from datetime import datetime
warnings.filterwarnings('ignore')

#%% Inputs
source='data/Scan_info.xlsx'

#%% Initialization
scan_info=pd.read_excel(source)
now=datetime.strftime(datetime.now(),'%Y%m%d_%H%M')
LiS.mkdir('data/'+now+'_'+os.path.basename(source)[:-5])

#%% Main
for id_scan in range(len(scan_info.index)):

    if scan_info['Run'].iloc[id_scan]:
        scan_name=scan_info['Scan name'].iloc[id_scan]
        vec_dazi=[np.float64(d) for d in str(scan_info['vec_dazi [deg]'].iloc[id_scan]).split(';')]
        vec_azi_max=[np.float64(d) for d in str(scan_info['vec_azi_max [deg]'].iloc[id_scan]).split(';')]
        vec_dele=[np.float64(d) for d in str(scan_info['vec_dele [deg]'].iloc[id_scan]).split(';')]
        vec_ele_max=[np.float64(d) for d in str(scan_info['vec_ele_max [deg]'].iloc[id_scan]).split(';')]
        d_cartesian=scan_info['d_cartesian'].iloc[id_scan]
        max_cartesian=scan_info['max_cartesian'].iloc[id_scan]
        Dn0=[np.float64(d) for d in scan_info['Dn0 [m]'].iloc[id_scan].split(';')]
        dr=scan_info['dr [m]'].iloc[id_scan]
        rmin=scan_info['rmin [m]'].iloc[id_scan]
        rmax=scan_info['rmax [m]'].iloc[id_scan]
        mins=[np.float64(d) for d in scan_info['mins [m]'].iloc[id_scan].split(';')]
        maxs=[np.float64(d) for d in scan_info['maxs [m]'].iloc[id_scan].split(';')]
        sigma=scan_info['sigma'].iloc[id_scan]
        T_tot=scan_info['T_tot [s]'].iloc[id_scan]
        tau=scan_info['tau [s]'].iloc[id_scan]
        sampl_time=scan_info['sampl_time [s]'].iloc[id_scan]
        U=scan_info['U [m/s]'].iloc[id_scan]
        
        r=np.arange(rmin,rmax,dr)
        
        epsilon1=np.zeros((len(vec_dazi),len(vec_azi_max),len(vec_dele),len(vec_ele_max)))+np.nan
        epsilon2=epsilon1.copy()
        epsilon3=epsilon1.copy()
        
        i_dazi=0
        for dazi in vec_dazi:
            i_azi_max=0
            for azi_max in vec_azi_max:
                theta=np.arange(-azi_max,azi_max+dazi/2+10**-10,dazi+10**-10)
                i_dele=0
                for dele in vec_dele:
                    i_ele_max=0
                    if (d_cartesian==False and i_dele==i_dazi) or d_cartesian==True:         
                        for ele_max in vec_ele_max:
                            if (max_cartesian==False and i_ele_max==i_azi_max) or  max_cartesian==True:     
                                    beta=np.arange(-ele_max,ele_max+dele/2+10**-10,dele+10**-10)
                                
                                    e1,e2,e3,X2,Dd,excl,T,x1=LiS.Pareto_v3(r,theta,beta,mins,maxs,Dn0,sigma,T_tot,tau,sampl_time,U)
                                    epsilon1[i_dazi,i_azi_max,i_dele,i_ele_max]=np.round(e1,1)
                                    epsilon2[i_dazi,i_azi_max,i_dele,i_ele_max]=np.round(e2,1)
                                    epsilon3[i_dazi,i_azi_max,i_dele,i_ele_max]=np.round(e3,1)
                            i_ele_max+=1
                    i_dele+=1
                i_azi_max+=1
            i_dazi+=1
            
        #Output
        Output = xr.Dataset({
            'epsilon1': xr.DataArray(
                        data   = epsilon1,   # enter data here
                        dims   = ['dazi','azi_max','dele','ele_max'],
                        coords = {'dazi': vec_dazi,'azi_max': vec_azi_max,'dele': vec_dele,'ele_max': vec_ele_max},
                        attrs  = {
                            '_FillValue': 'Nan',
                            'units'     : '%',
                            'description': 'Undersampled domain ratio'
                            }
                        ),
            'epsilon2': xr.DataArray(
                        data   = epsilon2,   # enter data here
                        dims   = ['dazi','azi_max','dele','ele_max'],
                        coords = {'dazi': vec_dazi,'azi_max': vec_azi_max,'dele': vec_dele,'ele_max': vec_ele_max},
                        attrs  = {
                            '_FillValue': 'Nan',
                            'units'     : '%',
                            'description': 'Error on the mean'
                            }
                        ),
            'epsilon3': xr.DataArray(
                        data   = epsilon3,   # enter data here
                        dims   = ['dazi','azi_max','dele','ele_max'],
                        coords = {'dazi': vec_dazi,'azi_max': vec_azi_max,'dele': vec_dele,'ele_max': vec_ele_max},
                        attrs  = {
                            '_FillValue': 'Nan',
                            'units'     : '%',
                            'description': 'Unresolved energy ratio'
                            }
                        ),
            'd_cartesian': xr.DataArray(
                        data   = d_cartesian,   # enter data here
                        attrs  = {'units':'boolean',
                                  'description': 'If TRUE: azimuth and elevation resolution are tested independently. If FALSE, azimuth and elevation resolution are tested together'}
                        ),
            'max_cartesian': xr.DataArray(
                        data   = max_cartesian,   # enter data here
                        attrs  = {'units':'boolean',
                                  'description': 'If TRUE: azimuth and elevation limits are tested independently. If FALSE, azimuth and elevation limits are tested together'}
                        ),
            'Dn0': xr.DataArray(
                        data   = Dn0,   # enter data here
                        attrs  = {'units':'m',
                                  'description': 'Fundamental half-wavelength'}
                        ),
            'rmin': xr.DataArray(
                        data   = rmin,   # enter data here
                        attrs  = {'units':'m',
                                  'description': 'Minimum lidar range'}
                        ),
            'rmax': xr.DataArray(
                        data   = rmax,   # enter data here
                        attrs  = {'units':'m',
                                  'description': 'Maximum lidar range'}
                        ),
            'dr': xr.DataArray(
                        data   = dr,   # enter data here
                        attrs  = {'units':'m',
                                  'description': 'Lidar range date length'}
                        ),
            'mins': xr.DataArray(
                        data   = mins,   # enter data here
                        attrs  = {'units':'m',
                                  'description': 'Lower boundary of the domain (x,y,z)'}
                        ),
            'maxs': xr.DataArray(
                        data   = maxs,   # enter data here
                        attrs  = {'units':'m',
                                  'description': 'Upper boundary of the domain (x,y,z)'}
                        ),
            'sigma': xr.DataArray(
                        data   = sigma,   # enter data here
                        attrs  = {'units':'none',
                                  'description': 'Barnes scheme smoothing parameter (1/4 recommended)'}
                        ),
            'T_tot': xr.DataArray(
                        data   = T_tot,   # enter data here
                        attrs  = {'units':'s',
                                  'description': 'Total duration of the experiment'}
                        ),
            'tau': xr.DataArray(
                        data   = tau,   # enter data here
                        attrs  = {'units':'s',
                                  'description': 'Mean integral timescale of the flow'}
                        ),
            'sampl_time': xr.DataArray(
                        data   = sampl_time,   # enter data here
                        attrs  = {'units':'s',
                                  'description': 'Sampling time of the lidar for each beam'}
                        ),
            'U': xr.DataArray(
                        data   = U,   # enter data here
                        attrs  = {'units':'s',
                                  'description': 'Mean wind speed'}
                        )},
            attrs = {'Contact': 'stefano.letizia@nrel.gov',
                      'Description':'Pareto front result from LiSBOA algorithm',
                      'Source':source,
                      'Scan name':scan_name})
        
        Output.to_netcdf('data/'+now+'_'+os.path.basename(source)[:-5]+'/'+scan_name+'.nc')
        sys.stdout.flush()
        sys.stdout.write('\r')
        sys.stdout.write(scan_name+' completed')
        sys.stdout.write('\n')
        sys.stdout.flush()
