import numpy as np

def Kaimal(f,U,Lu=340.2):
    return (4*Lu/U)/(1+6*f*Lu/U)**(5/3)  

def cosd(x):
    return np.cos(x/180*np.pi)

def sind(x):
    return np.sin(x/180*np.pi)

def mid(x):
    return (x[:-1]+x[1:])/2

def mkdir(path):
    #makes recursively folder from path, no existing error
    import os
    path.replace('\\','/')
    folders=path.split('/')
    upper=''
    for f in folders:
        try:
            os.mkdir(upper+f)           
        except:
            pass
                
        upper+=f+'/'
        

def len2(x):
    if 'int' in str(type(x)) or 'float' in str(type(x)):
        return 1
    elif 'list' in str(type(x)) or 'array' in str(type(x))  or 'str' in str(type(x)) or 'series' in str(type(x)):
        return len(x)
    else:
        raise ValueError
        
def str_dec(num):
    if len2(num)>1:
        s=''
        for n in num:
            if n!=np.round(n):
                s=s+ str(n)+', '
            else:
                s=s+ str(int(n))+', '
        return s[:-2]
    else:
        if num!=np.round(num):
            return str(num)
        else:
            return str(int(num))

def save_all_fig(name,newfolder=False,resolution=300):
    from matplotlib import pyplot as plt
    mkdir('figures')
    if newfolder:
        mkdir('figures/'+name)
    figs = [plt.figure(n) for n in plt.get_fignums()]
    inc=0
    for fig in figs:
        if newfolder:
            fig.savefig('figures/'+name+'/'+'{i:02d}'.format(i=inc)+'.png',dpi=resolution, bbox_inches='tight')
        else:
            fig.savefig('figures/'+name+'_'+'{i:02d}'.format(i=inc)+'.png',dpi=resolution, bbox_inches='tight')
        inc+=1

def Pareto_v3(r,theta,beta,mins,maxs,Dn0,sigma,T_tot,tau,sampl_time,U):
    #03/31/2023 (v 2): new halo simulator
    #04/03/2023: LiSBOA with real scanning points, finalized
    #08/28/2023 (v 3): adapted for release
    
    f=np.float64(10)**np.arange(-5,3,0.01)#[Hz] range of frequencies
    
    theta1,beta1=np.meshgrid(theta,beta)
    theta1=theta1.ravel() 
    beta1=beta1.ravel() 
    
    #virtual lidar (with ideal geometry)
    [R1,TH1]=np.meshgrid(r,theta1)
    [R1,B1]=np.meshgrid(r,beta1)
    X=R1*cosd(B1)*cosd(TH1)
    Y=R1*cosd(B1)*sind(TH1)
    Z=R1*sind(B1)
    x1=[np.ravel(X),np.ravel(Y),np.ravel(Z)]
    
    T=sampl_time*len(np.ravel(theta1))
    L=np.floor(T_tot/T)
    
    if L>0:
        #LiSBOA
        X2,Dd,excl,avg,HOM=LiSBOA_v7(x1,mins,maxs,Dn0,sigma)

        #Pareto
        epsilon1=np.sum(excl)/np.size(excl)*100
            
        p=np.arange(1,L)
        f_nqs=1/(2*T)
        
        epsilon2=(1/L+2/L**2*np.sum((L-p)*np.exp(-T/tau*p)))**0.5*100
    
        ER=np.trapz(Kaimal(f[f<f_nqs],U,U*tau),x=f[f<f_nqs])
        epsilon3=(1-ER)*100
    else:
        epsilon1=epsilon2=epsilon3=X2=Dd=excl=np.nan
        
    return epsilon1,epsilon2,epsilon3,X2,Dd,excl,T,x1

def LiSBOA_v7(x_exp,mins,maxs,Dn0,sigma,max_iter=None,calculate_stats=False,f=None,order=None,R_max=3,grid_factor=0.25,tol_dist=0.1,max_Dd=1):
    #03/01/2021 (v 5): undersampling checkes in hypercubes instead of hyperspheres (faster)
    #03/02/2022: finalized
    #08/28/2023 (v 7): added HOM,handled 0 dimensions, finalized
    
    from scipy.special import gamma
    from scipy.interpolate import interpn
    import itertools
    import sys
    import time
         
    #outliers rejection
    n=len(Dn0)    
    if calculate_stats:
        real=~np.isnan(np.sum(np.array(x_exp),axis=0)+f)
        for j in range(n):
            x_exp[j]=x_exp[j][real]
        f=f[real]
    else:
        real=~np.isnan(np.sum(np.array(x_exp),axis=0))
        for j in range(n):
            x_exp[j]=x_exp[j][real]
    
    #Initialization
    t0=time.time()
    Dn0=np.array(Dn0)   
    Dn0[Dn0<0.91]=1
    N=len(x_exp[0])
    x=np.zeros((n,N))
    xc=np.zeros(n)
    X_bin=[];
    X_vec=[];
    X2=[]
    avg=None
    HOM=None
    
    #LiSBOA setup
    dx=grid_factor*Dn0
    R_max=R_max*sigma
    V=np.pi**(n/2)/gamma(n/2+1)*R_max**n
    NoD=np.ceil(np.log10(np.max(np.abs(np.array(x_exp)/tol_dist))))+1
    
    for j in range(n):
        xc[j]=np.min(x_exp)
        x[j]=(x_exp[j]-xc[j])/Dn0[j]       
        X_bin.append((np.arange(mins[j]-dx[j]/2,maxs[j]+dx[j]/2+10**-10,dx[j])-xc[j])/Dn0[j])
        X_vec.append(mid(X_bin[j]))

    X=np.meshgrid(*[X_vec[j] for j in range(n)], indexing='ij')
    
    for j in range(n):
        X2.append(X[j]*Dn0[j]+xc[j])
      
    w=np.zeros(np.shape(X[0]),dtype=object)
    sel=np.zeros(np.shape(X[0]),dtype=object)
    val=np.zeros(np.shape(X[0]),dtype=object)
    Dd=np.zeros(np.shape(X[0]))
    N_grid=X[0].size
    dist_inf=np.zeros(n)
    
    #weights
    for j in range(n):
        dist_inf[j]=np.ceil(R_max/(dx[j]/Dn0[j]))
        
    nodes=np.where(X[0])
    counter=0
    for i in zip(*[xx for xx in nodes]):
        distSq=0
        for j in range(n):
            distSq+=(x[j]-X[j][i])**2
        s=np.where(distSq<R_max**2)   
        if len(s)>0:
            w[i]=np.exp(-distSq[s]/(2*sigma**2))
       
        #local spacing
        if Dd[i]!=10^99:   
            if len(s[0])>1:                
                pos_uni=np.around(x[0][s]/tol_dist)*tol_dist
                for j in range(1,n):                
                    pos_uni+=np.around(x[j][s]/tol_dist)*tol_dist*(10**NoD)**j          
                N_uni= len(np.unique(np.array(pos_uni)))
                
                if N_uni>1:
                    Dd[i]=V**(1/n)/(N_uni**(1/n)-1)
                else:
                    Dd[i]=np.inf
            else:
                Dd[i]=np.inf
                
            ind_inf=[]
            if Dd[i]>max_Dd:
                for j in range(n):
                    i1=max(i[j]-dist_inf[j],0)
                    i2=min(i[j]+dist_inf[j],np.shape(X[0])[j])
                    ind_inf.append(np.arange(i1,i2).astype(int))                
                for i_inf in itertools.product(*[ii for ii in ind_inf]):
                    Dd[i_inf]=10^99
        #store
        sel[i]=s
        
        counter+=1
        if np.floor(counter/N_grid*100)>np.floor((counter-1)/N_grid*100):
            est_time=(time.time()-t0)/counter*(N_grid-counter)
            sys.stdout.write('\r LiSBOA:'+str(np.floor(counter/N_grid*100).astype(int))+'% done, '+str(round(est_time))+' s left.') 
    sys.stdout.write('\r                                                                         ')
    sys.stdout.flush()
    excl=Dd>max_Dd
                
    #stats
    if calculate_stats:
        avg=[]
        HOM=[]
        df=f
        for m in range(max_iter+1):
            WM=np.zeros(np.shape(X[0]))+np.nan
            WM_HOM=np.zeros(np.shape(X[0]))+np.nan
            sys.stdout.write('\r Iteration #'+str(m))
            sys.stdout.flush()
            for i in zip(*[xx for xx in nodes]):
                val[i]=f[s]
                if not excl[i]:
                    fs=np.array(df[sel[i]])
                    ws=np.array(w[i])
                    reals=~np.isnan(fs+ws)
                    if sum(reals)>0:
                        fs=fs[reals]
                        ws=ws[reals]                     
                        WM[i]=sum(np.multiply(fs,ws))/sum(ws)              
                        if m>0:
                            WM_HOM[i]=sum(np.multiply(fs**order,ws))/sum(ws)  
            if m==0:
                avg.append(WM+0)
            else:
                avg.append(avg[m-1]+WM)
                HOM.append(WM_HOM)
            
            df=f-interpn(tuple(X_vec),avg[m],np.transpose(x),bounds_error=False,fill_value=np.nan)

    return X2,Dd,excl,avg,HOM