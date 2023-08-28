#Plot results of scan optimization from the LiDAR Barnes Statistical Objective Analysis tool [Letizia et al., Atmos Meas Tech, 2021a,b; Letizia et al., WAKE conf., 2023]
#The copyright of this software belongs to Stefano Letizia from NREL. 
#The author's contribution must be acknoledged in all the pubblications that include results form this software.
#The author is not responsible for any erroneous and/or financially or scientific harmful results produced by the present software.
#08/28/2023: created, finalized

import os
import sys
sys.path.append('functions')
cd=os.path.dirname(__file__)
import LiSBOA_functions as LiS
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import warnings
import xarray as xr
import matplotlib as mpl
mpl.rcParams.update({'font.size': 14})

warnings.filterwarnings('ignore')
plt.close('all')

#%% Inputs
root='data/20230828_0958_Scan_info/'#location of Pareto results

#PPI case
# source='PPI.nc'#specific file name
# sel=[0.5,20,0,0]#sequence of dazi,azi_max,dele,ele_max that is plotted

# #first volumetric case
# source='Volumetric rect1.nc'#specific file name
# sel=[2,10,1,10]#sequence of dazi,azi_max,dele,ele_max that is plotted

#second volumetric case
# source='Volumetric rect2.nc'#specific file name
# sel=[2,15,2,5]#sequence of dazi,azi_max,dele,ele_max that is plotted

# #third volumetric case
source='Volumetric symm.nc'#specific file name
sel=[2,15,2,15]#sequence of dazi,azi_max,dele,ele_max that is plotted

save_fig=True

#graphics
cmap = cm.get_cmap('viridis')
markers=['o','v','s','p','*']

#%% Initialization
Data=xr.open_dataset(root+source)

epsilon1=Data['epsilon1'].values
epsilon2=Data['epsilon2'].values
epsilon3=Data['epsilon3'].values
vec_dazi=Data['dazi'].values
vec_azi_max=Data['azi_max'].values
vec_dele=Data['dele'].values
vec_ele_max=Data['ele_max'].values
Dn0=Data['Dn0'].values
dr=np.float64(Data['dr'].values)
rmin=np.float64(Data['rmin'].values)
rmax=np.float64(Data['rmax'].values)
mins=Data['mins'].values
maxs=Data['maxs'].values
sigma=np.float64(Data['sigma'].values)
T_tot=np.float64(Data['T_tot'].values)
tau=np.float64(Data['tau'].values)
sampl_time=np.float64(Data['sampl_time'].values)
U=np.float64(Data['U'].values)
d_cartesian=Data['d_cartesian'].values
max_cartesian=Data['max_cartesian'].values

r=np.arange(rmin,rmax,dr)

#%% Main
if len(sel)==4:
    assert np.sum(vec_dazi==sel[0])==1, 'azimuth resolution not available'
    assert np.sum(vec_azi_max==sel[1])==1, 'maximum azimuth not available'
    assert np.sum(vec_dele==sel[2])==1, 'elevation resolution not available'
    assert np.sum(vec_ele_max==sel[3])==1, 'maximum elevation not available'
    
    i_azi_max_sel=np.where(vec_azi_max==sel[1])[0][0]
    azi_max=sel[1]
    ele_max=sel[3]
    
    theta=np.arange(-azi_max,azi_max+sel[0]/2+10**-10,sel[0]+10**-10)
    beta= np.arange(-ele_max,ele_max+sel[2]/2+10**-10,sel[2]+10**-10)
    
    e1,e2,e3,X2,Dd,excl,T,x1=LiS.Pareto_v3(r,theta,beta,mins,maxs,Dn0,sigma,T_tot,tau,sampl_time,U)
    
    e1=np.round(e1,1)
    e2=np.round(e2,1)
    e3=np.round(e3,1)
    T=np.round(T,1)
        
    
    X_plot=np.transpose(X2[0],(1,2,0))
    Y_plot=np.transpose(X2[1],(1,2,0))
    Z_plot=np.transpose(X2[2],(1,2,0))
    excl_plot=np.transpose(excl,(1,2,0))
    Dd_plot=np.transpose(Dd,(1,2,0))
    Dd_plot[excl_plot]=np.nan

else:
    X_plot=np.nan
    Y_plot=np.nan
    Z_plot=np.nan
    Dd_plot=np.nan
    excl_plot=np.nan
    e1=np.nan
    e2=np.nan
    e3=np.nan
    i_azi_max_sel=0

#%% Plots
    
#Pareto 1
if d_cartesian==True and max_cartesian==True:
    plt.figure(figsize=(18,10))
    i_dele=0
    for dele in vec_dele:
        i_ele_max=0
        for ele_max in vec_ele_max:
            ax=plt.subplot(len(vec_ele_max),len(vec_dele),i_ele_max*(len(vec_dele))+i_dele+1)
            plt.plot(epsilon1.ravel(),epsilon2.ravel(),'.k',markersize=15,alpha=0.2)
            for i_dazi in range(len(vec_dazi)):
                for i_azi_max in range(len(vec_azi_max)):
                    p=plt.plot(epsilon1[i_dazi,i_azi_max,i_dele,i_ele_max],epsilon2[i_dazi,i_azi_max,i_dele,i_ele_max],'.',markersize=10,alpha=1,markeredgecolor='k',marker=markers[i_azi_max],color=cmap(i_dazi/(len(vec_dazi)-1)))
            
            if len(sel)==4:
                if dele==sel[2] and ele_max==sel[3]:
                    plt.plot(e1,e2,'.',markersize=12,alpha=1,markeredgecolor='r',marker=markers[i_azi_max_sel],markeredgewidth=2,fillstyle='none')
                
            for i_azi_max in range(len(vec_azi_max)):
                plt.plot(1000,1000,'.',markeredgecolor='k',marker=markers[i_azi_max],label=r'$\theta_{max} ='+LiS.str_dec(vec_azi_max[i_azi_max])+'^\circ$',color='k')
            
            for i_dazi in range(len(vec_dazi)):
                plt.plot(1000,1000,'.',marker='o',markeredgecolor='k',color=cmap(i_dazi/(len(vec_dazi)-1)),label=r'$\Delta \theta ='+LiS.str_dec(vec_dazi[i_dazi])+'^\circ$')
            i_ele_max+=1
        
            plt.xlim([0,101])
            plt.ylim([0,101])
            
            ax.grid(visible=True) 
            plt.title(r'$\Delta \beta ='+LiS.str_dec(dele)+'^\circ$, '+r'$\beta_{max}='+LiS.str_dec(ele_max)+'^\circ$')
            plt.xlabel(r'$\epsilon_I$ [%]')
            plt.ylabel(r'$\epsilon_{II}$ [%]')
        i_dele+=1

if d_cartesian==False and max_cartesian==True:
    plt.figure(figsize=(18,6))
    i_ele_max=0
    for ele_max in vec_ele_max:
        ax=plt.subplot(1,len(vec_ele_max),i_ele_max+1)
        plt.plot(epsilon1.ravel(),epsilon2.ravel(),'.k',markersize=15,alpha=0.2)
        
        for i_dazi in range(len(vec_dazi)):
            for i_azi_max in range(len(vec_azi_max)):
                p=plt.plot(epsilon1[i_dazi,i_azi_max,i_dazi,i_ele_max],epsilon2[i_dazi,i_azi_max,i_dazi,i_ele_max],'.',markersize=10,alpha=1,markeredgecolor='k',marker=markers[i_azi_max],color=cmap(i_dazi/(len(vec_dazi)-1)))
        
        if len(sel)==4:         
            if ele_max==sel[3]:
                plt.plot(e1,e2,'.',markersize=12,alpha=1,markeredgecolor='r',marker=markers[i_azi_max_sel],markeredgewidth=2,fillstyle='none')
            
        for i_azi_max in range(len(vec_azi_max)):
            plt.plot(1000,1000,'.',markeredgecolor='k',marker=markers[i_azi_max],label=r'$\theta_{max} ='+LiS.str_dec(vec_azi_max[i_azi_max])+'^\circ$',color='k')
        
        for i_dazi in range(len(vec_dazi)):
            plt.plot(1000,1000,'.',marker='o',markeredgecolor='k',color=cmap(i_dazi/(len(vec_dazi)-1)),label=r'$\Delta \theta ='+LiS.str_dec(vec_dazi[i_dazi])+r'^\circ$, $\Delta \beta ='+LiS.str_dec(vec_dele[i_dazi])+'^\circ$ ')
        i_ele_max+=1
        
        plt.xlim([0,101])
        plt.ylim([0,101])
        
        ax.grid(visible=True) 
        plt.title(r'$\beta_{max}='+LiS.str_dec(ele_max)+'^\circ$')
        plt.xlabel(r'$\epsilon_I$ [%]')
        plt.ylabel(r'$\epsilon_{II}$ [%]')

if d_cartesian==True and max_cartesian==False:
    plt.figure(figsize=(18,6))
    i_dele=0
    for dele in vec_dele:       
        ax=plt.subplot(1,len(vec_dele),i_dele+1)
        plt.plot(epsilon1.ravel(),epsilon2.ravel(),'.k',markersize=15,alpha=0.2)
        
        for i_dazi in range(len(vec_dazi)):
            for i_azi_max in range(len(vec_azi_max)):
                p=plt.plot(epsilon1[i_dazi,i_azi_max,i_dele,i_azi_max],epsilon2[i_dazi,i_azi_max,i_dele,i_azi_max],'.',markersize=10,alpha=1,markeredgecolor='k',marker=markers[i_azi_max],color=cmap(i_dazi/(len(vec_dazi)-1)))
        
        if len(sel)==4:        
            if dele==sel[2]:
                plt.plot(e1,e2,'.',markersize=12,alpha=1,markeredgecolor='r',marker=markers[i_azi_max_sel],markeredgewidth=2,fillstyle='none')
            
        for i_azi_max in range(len(vec_azi_max)):
            plt.plot(1000,1000,'.',markeredgecolor='k',marker=markers[i_azi_max],label=r'$\theta_{max} ='+LiS.str_dec(vec_azi_max[i_azi_max])+r'^\circ$, $\beta_{max} ='+LiS.str_dec(vec_ele_max[i_azi_max])+'^\circ$',color='k')
        
        for i_dazi in range(len(vec_dazi)):
            plt.plot(1000,1000,'.',marker='o',markeredgecolor='k',color=cmap(i_dazi/(len(vec_dazi)-1)),label=r'$\Delta \theta ='+LiS.str_dec(vec_dazi[i_dazi])+r'^\circ$')
        i_dele+=1
        
        plt.xlim([0,101])
        plt.ylim([0,101])
        
        ax.grid(visible=True) 
        plt.title(r'$\Delta \beta ='+LiS.str_dec(dele)+'^\circ$')
        plt.xlabel(r'$\epsilon_I$ [%]')
        plt.ylabel(r'$\epsilon_{II}$ [%]')

if d_cartesian==False and max_cartesian==False:
    plt.figure(figsize=(12,10))
    for i_dazi in range(len(vec_dazi)):
        for i_azi_max in range(len(vec_azi_max)):
            p=plt.plot(epsilon1[i_dazi,i_azi_max,i_dazi,i_azi_max],epsilon2[i_dazi,i_azi_max,i_dazi,i_azi_max],'.',markersize=10,alpha=1,markeredgecolor='k',marker=markers[i_azi_max],color=cmap(i_dazi/(len(vec_dazi)-1)))
    
    plt.plot(e1,e2,'.',markersize=12,alpha=1,markeredgecolor='r',marker=markers[i_azi_max_sel],markeredgewidth=2,fillstyle='none')
    
    for i_azi_max in range(len(vec_azi_max)):
        plt.plot(1000,1000,'.',markeredgecolor='k',marker=markers[i_azi_max],label=r'$\theta_{max} ='+LiS.str_dec(vec_azi_max[i_azi_max])+r'^\circ$, $\beta_{max} ='+LiS.str_dec(vec_ele_max[i_azi_max])+'^\circ$',color='k')
    
    for i_dazi in range(len(vec_dazi)):
        plt.plot(1000,1000,'.',marker='o',markeredgecolor='k',color=cmap(i_dazi/(len(vec_dazi)-1)),label=r'$\Delta \theta = '+LiS.str_dec(vec_dazi[i_dazi])+r'^\circ$, $\Delta \beta='+LiS.str_dec(vec_dele[i_dazi])+'^\circ$')
    
    plt.xlim([0,101])
    plt.ylim([0,101])
    
    plt.gca().grid(visible=True) 
    plt.xlabel(r'$\epsilon_I$ [%]')
    plt.ylabel(r'$\epsilon_{II}$ [%]')
    
plt.legend().set_draggable(state=True)
plt.tight_layout()

#Pareto 2
if d_cartesian==True and max_cartesian==True:
    plt.figure(figsize=(18,10))
    i_dele=0
    for dele in vec_dele:
        i_ele_max=0
        for ele_max in vec_ele_max:
            
            ax=plt.subplot(len(vec_ele_max),len(vec_dele),i_ele_max*(len(vec_dele))+i_dele+1)
            plt.plot(epsilon1.ravel(),epsilon3.ravel(),'.k',markersize=15,alpha=0.2)
            
            for i_dazi in range(len(vec_dazi)):
                for i_azi_max in range(len(vec_azi_max)):
                    p=plt.plot(epsilon1[i_dazi,i_azi_max,i_dele,i_ele_max],epsilon3[i_dazi,i_azi_max,i_dele,i_ele_max],'.',markersize=10,alpha=1,markeredgecolor='k',marker=markers[i_azi_max],color=cmap(i_dazi/(len(vec_dazi)-1)))
            
            if len(sel)==4:
                if dele==sel[2] and ele_max==sel[3]:
                    plt.plot(e1,e3,'.',markersize=12,alpha=1,markeredgecolor='r',marker=markers[i_azi_max_sel],markeredgewidth=2,fillstyle='none')
                
            for i_azi_max in range(len(vec_azi_max)):
                plt.plot(1000,1000,'.',markeredgecolor='k',marker=markers[i_azi_max],label=r'$\theta_{max} ='+LiS.str_dec(vec_azi_max[i_azi_max])+'^\circ$',color='k')
            
            for i_dazi in range(len(vec_dazi)):
                plt.plot(1000,1000,'.',marker='o',markeredgecolor='k',color=cmap(i_dazi/(len(vec_dazi)-1)),label=r'$\Delta \theta ='+LiS.str_dec(vec_dazi[i_dazi])+'^\circ$')
            i_ele_max+=1
        
            plt.xlim([0,101])
            plt.ylim([0,101])
            
            ax.grid(visible=True) 
            plt.title(r'$\Delta \beta ='+LiS.str_dec(dele)+'^\circ$, '+r'$\beta_{max}='+LiS.str_dec(ele_max)+'^\circ$')
            plt.xlabel(r'$\epsilon_I$ [%]')
            plt.ylabel(r'$\epsilon_{III}$ [%]')
        i_dele+=1

if d_cartesian==False and max_cartesian==True:
    plt.figure(figsize=(18,6))
    i_ele_max=0
    for ele_max in vec_ele_max:
        ax=plt.subplot(1,len(vec_ele_max),i_ele_max+1)
        plt.plot(epsilon1.ravel(),epsilon3.ravel(),'.k',markersize=15,alpha=0.2)
        
        for i_dazi in range(len(vec_dazi)):
            for i_azi_max in range(len(vec_azi_max)):
                p=plt.plot(epsilon1[i_dazi,i_azi_max,i_dazi,i_ele_max],epsilon3[i_dazi,i_azi_max,i_dazi,i_ele_max],'.',markersize=10,alpha=1,markeredgecolor='k',marker=markers[i_azi_max],color=cmap(i_dazi/(len(vec_dazi)-1)))
        
        if len(sel)==4:        
            if ele_max==sel[3]:
                plt.plot(e1,e3,'.',markersize=12,alpha=1,markeredgecolor='r',marker=markers[i_azi_max_sel],markeredgewidth=2,fillstyle='none')
            
        for i_azi_max in range(len(vec_azi_max)):
            plt.plot(1000,1000,'.',markeredgecolor='k',marker=markers[i_azi_max],label=r'$\theta_{max} ='+LiS.str_dec(vec_azi_max[i_azi_max])+'^\circ$',color='k')
        
        for i_dazi in range(len(vec_dazi)):
            plt.plot(1000,1000,'.',marker='o',markeredgecolor='k',color=cmap(i_dazi/(len(vec_dazi)-1)),label=r'$\Delta \theta ='+LiS.str_dec(vec_dazi[i_dazi])+r'^\circ$, $\Delta \beta ='+LiS.str_dec(vec_dele[i_dazi])+'^\circ$ ')
        i_ele_max+=1
        
        plt.xlim([0,101])
        plt.ylim([0,101])
        
        ax.grid(visible=True) 
        plt.title(r'$\beta_{max}='+LiS.str_dec(ele_max)+'^\circ$')
        plt.xlabel(r'$\epsilon_I$ [%]')
        plt.ylabel(r'$\epsilon_{III}$ [%]')

if d_cartesian==True and max_cartesian==False:
    plt.figure(figsize=(18,6))
    i_dele=0
    for dele in vec_dele:       
        ax=plt.subplot(1,len(vec_dele),i_dele+1)
        plt.plot(epsilon1.ravel(),epsilon3.ravel(),'.k',markersize=15,alpha=0.2)
        
        for i_dazi in range(len(vec_dazi)):
            for i_azi_max in range(len(vec_azi_max)):
                p=plt.plot(epsilon1[i_dazi,i_azi_max,i_dele,i_azi_max],epsilon3[i_dazi,i_azi_max,i_dele,i_azi_max],'.',markersize=10,alpha=1,markeredgecolor='k',marker=markers[i_azi_max],color=cmap(i_dazi/(len(vec_dazi)-1)))
        
        if len(sel)==4:        
            if dele==sel[2]:
                plt.plot(e1,e3,'.',markersize=12,alpha=1,markeredgecolor='r',marker=markers[i_azi_max_sel],markeredgewidth=2,fillstyle='none')
            
        for i_azi_max in range(len(vec_azi_max)):
            plt.plot(1000,1000,'.',markeredgecolor='k',marker=markers[i_azi_max],label=r'$\theta_{max} ='+LiS.str_dec(vec_azi_max[i_azi_max])+r'^\circ$, $\beta_{max} ='+LiS.str_dec(vec_ele_max[i_azi_max])+'^\circ$',color='k')
        
        for i_dazi in range(len(vec_dazi)):
            plt.plot(1000,1000,'.',marker='o',markeredgecolor='k',color=cmap(i_dazi/(len(vec_dazi)-1)),label=r'$\Delta \theta ='+LiS.str_dec(vec_dazi[i_dazi])+r'^\circ$')
        i_dele+=1
        
        plt.xlim([0,101])
        plt.ylim([0,101])
        
        ax.grid(visible=True) 
        plt.title(r'$\Delta \beta ='+LiS.str_dec(dele)+'^\circ$')
        plt.xlabel(r'$\epsilon_I$ [%]')
        plt.ylabel(r'$\epsilon_{III}$ [%]')

if d_cartesian==False and max_cartesian==False:
    plt.figure(figsize=(12,10))
    for i_dazi in range(len(vec_dazi)):
        for i_azi_max in range(len(vec_azi_max)):
            p=plt.plot(epsilon1[i_dazi,i_azi_max,i_dazi,i_azi_max],epsilon3[i_dazi,i_azi_max,i_dazi,i_azi_max],'.',markersize=10,alpha=1,markeredgecolor='k',marker=markers[i_azi_max],color=cmap(i_dazi/(len(vec_dazi)-1)))
    
    plt.plot(e1,e3,'.',markersize=12,alpha=1,markeredgecolor='r',marker=markers[i_azi_max_sel],markeredgewidth=2,fillstyle='none')
    
    for i_azi_max in range(len(vec_azi_max)):
        plt.plot(1000,1000,'.',markeredgecolor='k',marker=markers[i_azi_max],label=r'$\theta_{max} = '+LiS.str_dec(vec_azi_max[i_azi_max])+r'^\circ$, $\beta_{max}='+LiS.str_dec(vec_ele_max[i_azi_max])+'^\circ$',color='k')
    
    for i_dazi in range(len(vec_dazi)):
        plt.plot(1000,1000,'.',marker='o',markeredgecolor='k',color=cmap(i_dazi/(len(vec_dazi)-1)),label=r'$\Delta \theta = '+LiS.str_dec(vec_dazi[i_dazi])+r'^\circ$, $\Delta \beta='+LiS.str_dec(vec_dele[i_dazi])+'^\circ$')
    
    plt.xlim([0,101])
    plt.ylim([0,101])
    
    plt.gca().grid(visible=True) 
    plt.xlabel(r'$\epsilon_I$ [%]')
    plt.ylabel(r'$\epsilon_{III}$ [%]')
    
plt.legend().set_draggable(state=True)
plt.tight_layout()

#scan geometry
fig = plt.figure(figsize=(18,9))
ds=np.max(np.array([maxs[0]-mins[0],maxs[1]-mins[1],maxs[2]-mins[2]]))/100
ax = fig.add_subplot(121,projection='3d')
sel_x=(x1[0]>mins[0]-ds)*(x1[0]<maxs[0]+ds)
sel_y=(x1[1]>mins[1]-ds)*(x1[1]<maxs[1]+ds)
sel_z=(x1[2]>mins[2]-ds)*(x1[2]<maxs[2]+ds)
sel_plot=sel_x*sel_y*sel_z

sc=ax.scatter(x1[0][sel_plot],x1[1][sel_plot],x1[2][sel_plot],s=3, c='r', alpha=1,vmin=0,vmax=1,marker='x')

ax.set_xlim([mins[0]-ds,maxs[0]+ds])
ax.set_ylim([mins[1]-ds,maxs[1]+ds])
ax.set_zlim([mins[2]-ds,maxs[2]+ds])
ax.set_box_aspect([maxs[0]-mins[0]+2*ds,maxs[1]-mins[1]+2*ds,maxs[2]-mins[2]+2*ds])

ax.set_xlabel(r'$x/D$')
ax.set_ylabel(r'$y/D$')
ax.set_zlabel(r'$z/D$')
if (mins[2]-maxs[2])==0:
    ax.set_zticks([mins[2]])

ax.xaxis.labelpad=30

ax = fig.add_subplot(122,projection='3d')
sc=ax.scatter(X_plot[~excl_plot],Y_plot[~excl_plot],Z_plot[~excl_plot],s=2, c=Dd_plot[~excl_plot], alpha=0.5,vmin=0,vmax=1)
ax.set_xlim([mins[0]-ds,maxs[0]+ds])
ax.set_ylim([mins[1]-ds,maxs[1]+ds])
ax.set_zlim([mins[2]-ds,maxs[2]+ds])
ax.set_box_aspect([maxs[0]-mins[0]+2*ds,maxs[1]-mins[1]+2*ds,maxs[2]-mins[2]+2*ds])
plt.colorbar(sc,label=r'$\Delta \tilde{d}$',location='top')

ax.set_xlabel(r'$x/D$')
ax.set_ylabel(r'$y/D$')
ax.set_zlabel(r'$z/D$')
if (mins[2]-maxs[2])==0:
    ax.set_zticks([mins[2]])

ax.xaxis.labelpad=30

ax=fig.add_axes([0.05,0.4,0.1,0.1])
plt.text(0,0,r'$\Delta\theta = '+LiS.str_dec(sel[0])+' ^\circ$\n'+
             r'$\theta_{max}='+LiS.str_dec(sel[1])+' ^\circ$ \n'+
             r'$\Delta\beta = '+LiS.str_dec(sel[-2])+' ^\circ$\n'+
             r'$\beta_{max}='+LiS.str_dec(sel[-1])+' ^\circ$ \n'+
             r'$\Delta r='+LiS.str_dec(dr)+'$ m \n'+
             r'$\Delta n_0 = ['+LiS.str_dec(Dn0)+']$ m \n'+
             r'$T='+LiS.str_dec(T)+'$ s \n'+
             r'$T_{tot}='+LiS.str_dec(T_tot)+'$ s \n'+
             r'$\tau='+LiS.str_dec(tau)+'$ s \n'+
             r'$\Delta t='+LiS.str_dec(sampl_time)+'$ s \n'+
             r'$U_\infty='+LiS.str_dec(U)+'$ m/s \n'+
             r'$\epsilon_I='+LiS.str_dec(e1)+'\%$ \n'+
             r'$\epsilon_{II} ='+LiS.str_dec(e2)+'\%$ \n'+
             r'$\epsilon_{III} ='+LiS.str_dec(e3)+'\%$',fontsize='small') 
plt.axis('off')

if save_fig:
    LiS.mkdir('figures/'+root[5:]+Data.attrs['Scan name'])
    LiS.save_all_fig(root[5:]+'/'+Data.attrs['Scan name']+'/'+Data.attrs['Scan name']+'_LiSBOA')