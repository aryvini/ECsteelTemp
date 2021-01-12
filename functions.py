import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.ticker as tickers
import seaborn as sns
from math import *
from scipy.interpolate import interp1d



def param_reduction(fp,fy,Ea,temp):

    #Reduction factors deppending on temperature 
    # See Table 3.1
    X = [0,20,100,200,300,400,500,600,700,800,900,1000,1100,1200]
    ky = [1,1,1,1,1,1,0.78,0.47,0.23,0.11,0.06,0.04,0.02,0]
    kp = [1,1,1,0.807,0.613,0.42,0.36,0.18,0.075,0.050,0.0375,0.025,0.0125,0]
    kE = [1,1,1,0.9,0.8,0.7,0.6,0.31,0.13,0.09,0.0675,0.0450,0.0225,0]

    #Interpolation for intermediate values
    fky = interp1d(X,ky,kind='linear')
    fkp = interp1d(X,kp,kind='linear')
    fkE = interp1d(X,kE,kind='linear')

    #Applying the reduction factor to fyi,fpi,Eai to get theta values e.g for each temperature

    fyr = fy * fky(temp) * 1
    fpr = fp * fkp(temp) * 1
    Ear = Ea * fkE(temp) * 1

    return {"fp":fpr,"fy":fyr,"Ea":Ear}


def stress(e,temp,fp,fy,Ea):
    '''
    inputs:
    e: current strain value;
    temp: current temperature value in Celsius;
    fp: strength on the proportional limit at 20 C in Pa;
    fy: yield strength at 20 C in Pa
    Ea: Young modulus at 20 C in Pa
    '''

    #Definition of parameters
    
    reduced = param_reduction(fp,fy,Ea,temp)

    fp=reduced['fp']
    fy=reduced['fy']
    Ea= reduced['Ea']


    ep=fp/Ea
    ey=0.02
    et=0.15
    eu=0.2

    ## Auxiliary functions 
    def Fc(fy,fp,ey,ep,Ea):
        return ((fy-fp)**2)/((ey-ep)*Ea-(2*(fy-fp)))

    def Fa2(ey,ep,Ea):
        c=Fc(fy,fp,ey,ep,Ea)
        return (ey-ep)*(ey-ep+(c/Ea))

    def Fb2(ey,ep,Ea):
        c=Fc(fy,fp,ey,ep,Ea)
        return c*(ey-ep)*Ea+c**2

    #Strain ranges 
    if e<=ep:
        return e*Ea

    elif ep<e<ey:
        
        c=Fc(fy,fp,ey,ep,Ea)
        a2=Fa2(ey,ep,Ea)
        b2=Fb2(ey,ep,Ea)
        a=sqrt(a2)
        b=sqrt(b2)

        return fp-c+(b/a)*(a2-(ey-e)**2)**0.5
    
    elif ey<=e<=et:
        return fy
    
    elif et<e<eu:
        return fy*(1-(e-et)/(eu-et))
    
    elif e==eu:
        return 0
    else:
        pass


def stress_hardening(e,temp,fp,fy,Ea):
    '''
    inputs:
    e: current strain value;
    temp: current temperature value in Celsius;
    fp: strength on the proportional limit at 20 C in Pa;
    fy: yield strength at 20 C in Pa
    Ea: Young modulus at 20 C in Pa
    '''

    #item A(2) states that hardening valid for temperatures bellow 400 C
    if (temp >=400) or (e<=0.02):
        return stress(e,temp,fp,fy,Ea)
        
    else:

        #Definition of parameters
        
        reduced = param_reduction(fp,fy,Ea,temp)

        fp=reduced['fp']
        fy=reduced['fy']
        Ea=reduced['Ea']


        #fu definition

        if temp < 300:
            fu = 1.25*fy

        elif 300 <= temp < 400:
            fu = fy*(2-0.0025 * temp)

        elif temp >= 400:
            fu = fy
        else:
            print('fu definition failed')

        
        # #strain ranges
        # if e<=0.02:
        #     return stress(e,temp,fp,fy,Ea)

        if 0.02 < e < 0.04:
            return 50*(fu-fy)*e+2*fy-fu

        elif 0.04<=e<=0.15:
            return fu
        
        elif 0.15<=e<=0.2:
            return fu*(1-20*(e-0.15))
        
        elif e>=0.2:
            return 0
        