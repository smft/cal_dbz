"""
calculate reflectivity
@author:qzhang
"""
import numpy as np
from netCDF4 import Dataset
from matplotlib import pyplot as plt

def cal_dbz(t,t00,p,pb,qvp,qra,qsn,qgr):
    cp=1004
    rd=287.05
    prs=(p+pb)/100
    tmk=(t+t00)/((1000/prs)**(rd/cp))
    PI=np.pi
    r1=1.e-15
    ron=8.e6
    ron2=1.e10
    son=2.e7
    gon=5.e7
    ron_min=8.e6
    ron_qr0=0.00010
    ron_delqr0=0.25*ron_qr0
    ron_const1r=(ron2-ron_min)*0.5
    ron_const2r=(ron2+ron_min)*0.5
    gamma_seven=720.
    RHOWAT=1000
    rho_r=RHOWAT
    """1000. kg m^-3"""
    rho_s=100.
    """kg m^-3"""
    rho_g=400.
    """kg m^-3"""
    rhoair=(p+pb)/(rd*(1+0.608*qvp)*tmk)
    alpha=0.224
    celkel=273.15
    """"absolute 0 degree"""

    factor_r=gamma_seven*1.e18*(1./(PI*rho_r))**1.75
    dbz=np.zeros(np.shape(qra))
    tmc=0

    i=0
    while i<np.shape(qra)[0]:
        j=0
        while j<np.shape(qra)[1]:
            qra[i,j]=max(qra[i,j],0)
            qvp[i,j]=max(qvp[i,j],0)
            qsn[i,j]=max(qsn[i,j],0)
            qgr[i,j]=max(qgr[i,j],0)
            if tmk[i,j]<celkel:
                qra[i,j]=0
                qsn[i,j]=qra[i,j]
                factor_s=gamma_seven*1.e18*(1./(PI*rho_s))**1.75*(rho_s/RHOWAT)**2*alpha
                factor_g=gamma_seven*1.e18*(1./(PI*rho_g))**1.75*(rho_g/RHOWAT)**2*alpha
            else:
                factor_s=gamma_seven*1.e18*(1./(PI*rho_s))**1.75*(rho_s/RHOWAT)**2
                factor_g=gamma_seven*1.e18*(1./(PI*rho_g))**1.75*(rho_g/RHOWAT)**2
            tmc=min(-0.001,tmk[i,j]-celkel)
            sonv=min(2e08,2e06*np.exp(-0.12*tmc))
            gonv=gon
            if qgr[i,j]>=r1:
                gonv=2.38*(PI*rho_g/(rhoair*qgr[i,j]))**0.92
                gonv=max(1.e4,min(gonv,gon))
            ronv=ron
            if qra[i,j]>=r1:
                ronv=ron_const1r*np.tanh((ron_qr0-qra[i,j])/ron_delqr0)+ron_const2r
            dbz[i,j]=max(factor_r*(rhoair[i,j]*qra[i,j])**1.75/ronv**.75+\
                factor_s*(rhoair[i,j]*qsn[i,j])**1.75/sonv**.75+\
                factor_g*(rhoair[i,j]*qgr[i,j])**1.75/gonv**.75,0.001)
            dbz[i,j]=10*np.log10(dbz[i,j])
            j+=1
        i+=1
    return dbz

"""TEST!!!TEST"""
t=Dataset('./wrfout_d04_2011-07-18_06:00:00',mode='r').variables['T'][0][5][:][:]
t00=Dataset('./wrfout_d04_2011-07-18_06:00:00',mode='r').variables['T00'][:]
p=Dataset('./wrfout_d04_2011-07-18_06:00:00',mode='r').variables['P'][0][5][:][:]
pb=Dataset('./wrfout_d04_2011-07-18_06:00:00',mode='r').variables['PB'][0][5][:][:]
qvp=Dataset('./wrfout_d04_2011-07-18_06:00:00',mode='r').variables['QVAPOR'][0][5][:][:]
qra=Dataset('./wrfout_d04_2011-07-18_06:00:00',mode='r').variables['QRAIN'][0][5][:][:]
qsn=Dataset('./wrfout_d04_2011-07-18_06:00:00',mode='r').variables['QSNOW'][0][5][:][:]
size=np.shape(t)
qgr=np.zeros([size[0],size[1]])
plt.contourf(cal_dbz(t,t00,p,pb,qvp,qra,qsn,qgr))
plt.colorbar()
plt.show()
