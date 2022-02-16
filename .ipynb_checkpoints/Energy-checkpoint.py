#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 12:56:37 2021

@author: matias
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# Programa que calcula las energias del atomo de hidrogeno en funcion
# del campo magnetico. El valor de p indica si se quiere aplicar
# correccion por masa finita del nucleo (p=1) o no (p=0).

# Igual que Energy.py excepto que tiene una ligera modificacion bajo la
# condicion xa>xb (cruce de niveles), de manera que no se anule la diferencia
# eb-ea (xs indeterminado para nu mayor o igual que 4).

def Energy(beta,m,ny,p):

    import numpy as np
    import math
    from Corresp_01 import Corresp_01
    
    n,l=Corresp_01(m,ny) 
        
    x=math.log10(beta)

    en=-2*math.log10(n)

    if (ny!=0):    
        var=(ny+1)/2
        eny=-2*math.log10(int(var)) 

   
    if ny==0:

        b0=np.array([-8.51584e-1,7.86224e-1,9.50091e-2,5.73409e-1,1.26974,1.70505e-1])
        b1=np.array([-2.90213,-2.28335,-1.97412,-1.54066,-3.78015e-1,5.16550e-2])
        b2=np.array([1.01555,9.37692e-1,1.00523,9.77581e-1,9.10852e-1,6.92991e-1])
        xa=b0[0]+b1[0]*(math.log10(1+abs(m)))**b2[0]
        xb=b0[1]+b1[1]*(math.log10(1+abs(m)))**b2[1]
        ea=b0[2]+b1[2]*(math.log10(1+abs(m)))**b2[2]
        eb=b0[3]+b1[3]*(math.log10(1+abs(m)))**b2[3]
        ec=b0[4]+b1[4]*(math.log10(1+abs(m)))**b2[4]
        ep=b0[5]+b1[5]*(math.log10(1+abs(m)))**b2[5]
        a1=(en-ea)/(ea-1)
        a2=-ep*((1+a1)**2/((en-1)*a1))

        if(x<xa):
            e=1+(en-1)/(1+a1*math.exp(a2*(x-xa-0.1*(x-xa)**2)))
        elif (xa<x<xb):
            e=ea+(eb-ea)*((x-xa)/(xb-xa))**1.22
        else:
            e=eb+(ec-eb)*((x-xb)/(3-xb))**0.92 
    
    elif ny==1:

        b0=np.array([-6.66302e-1,-3.61037e-1,2.13743e-1]) 
        b1=np.array([-1.50237,-9.80935e-1,2.23000e-1]) 
        b2=np.array([1.17845,1.22078,8.82388e-1])
        xb=b0[0]+b1[0]*(math.log10(1+abs(m)))**b2[0]
        eb=b0[1]+b1[1]*(math.log10(1+abs(m)))**b2[1]
        ep=b0[2]+b1[2]*(math.log10(1+abs(m)))**b2[2]    
        n=abs(m)+2          
        en=-2*math.log10(n)      

        d=0.26  

        if abs(m)==1:  
            d=0.20  
        elif abs(m)==2:     
            d=0.22  
        elif abs(m)==3:   
            d=0.24
        
        if (x>xb):
            d=0
            
        a1=(en-eb)/(eb-eny)
        a2=-ep*(1+a1)**2/((en-eny)*a1) 
                
        e=eny+(en-eny)/(1+a1*math.exp(a2*(x-xb)*(1-d*(x-xb))))  

    elif ny==2:           
        b0=np.array([5.28777e-2,-4.52254e-1,1.19340e-1])   
        b1=np.array([-2.38204,-1.00281,2.96234e-1])
        b2=np.array([9.60364e-1,1.23880,1.03199])
        xb=b0[0]+b1[0]*[math.log10(1+abs(m))]**b2[0]
        eb=b0[1]+b1[1]*[math.log10(1+abs(m))]**b2[1]
        ep=b0[2]+b1[2]*[math.log10(1+abs(m))]**b2[2]
        
        if(x<=xb):
            n=abs(m)+2
            d=0.2
            a1=(en-eb)/(eb-eny)
            a2=-ep*(1+a1)**2/((en-eny)*a1) 
            e=eny+(en-eny)/(1+a1*math.exp(a2*(x-xb)*(1-d*(x-xb))))            
        else:    
            mp=min(abs(m),4)   
            epsi=-0.0125+0.030456*(math.log10(1+mp))**1.134 
            en1=-2*math.log10(1)   
            e=eb+(2/math.pi)*(eb-en1)*math.atan(((math.pi*ep*(x-xb))/(2*(eb-eny)))*(1+epsi*(x-xb)**2))
    
        e=float(e);
    
    elif ny==3:
        b0=np.array([-7.10984e-1,-7.90709e-1,1.17903e-1])
        b1=np.array([-1.78597,-7.84790e-1,2.60062e-1]) 
        b2=np.array([1.16795,1.36181,1.11196])
        xb=b0[0]+b1[0]*[math.log10(1+abs(m))]**b2[0]
        eb=b0[1]+b1[1]*[math.log10(1+abs(m))]**b2[1]
        ep=b0[2]+b1[2]*[math.log10(1+abs(m))]**b2[2]
        n=abs(m)+3   
        en=-2*math.log10(n)

        if abs(m)==0:
            d=0
        elif abs(m)==1:   
            d=0.08
        elif abs(m)==2:
            d=0.13
        elif abs(m)==3:
            d=0.15
        elif abs(m)==4:
            d=0.165
        else:   
            d=0.17

        if(x<=xb):
            a1=(en-eb)/(eb-eny)
            a2=-ep*(1+a1)**2/((en-eny)*a1) 
            e=eny+(en-eny)/(1+a1*math.exp(a2*(x-xb)*(1-d*(x-xb)))) 
        else:                      
            q=2.5
            e=eb+((eny-eb)*(x-xb))/(((eny-eb)/ep)**q+(x-xb)**q)**(1/q)
       
        e=float(e);
        
    else:

        aux=ny/2-int(ny/2)                           
        
        if(aux==0):
            
            t=2*math.sqrt((1/2)*(ny+1/2))-2
        
            tau=t-int(t)
            
            Delta=math.log10(ny)-math.log10(4)
            
            b0=np.array([1.1*tau-1.154902-2.087178*Delta**1.082710, \
                         -0.522879,-0.68-1.176143*Delta**0.8685913, \
                             -0.8867395-1.744739*Delta**1.095173, \
                                 0.02-0.034*tau+0.2/(ny**1.1+3), \
                                     0.3780437*ny**(-0.9572978)])
 
            b1=np.array([-0.01890508,-0.95-1.1*ny**(-0.4),-0.2-1.1* \
                         ny**(-0.4),-2.487767*ny**(-0.9652760), \
                         0.1438085,0.8265754*ny**(-0.9347425)])

            b2=np.array([0.5904491,0.85+1.1*ny**(-0.4),0.6+0.8* \
                         Delta**0.4,1.209001,1.596943,1.114659])
                        
        else:

            t=2*math.sqrt((1/2)*(ny-1/2))-2
            
            tau=t-int(t)

            Delta=math.log10(ny)-math.log10(5)

            b0=np.array([-tau-1.18-2.312886*Delta**0.7737455, \
                         -1.154902,-0.9558838-1.069160*Delta**(0.8065575), \
                         -1.12-1.707775*Delta**1.119483,0.013- \
                          0.034*tau+0.2/(ny**0.74+6),0.3480917* \
                          ny**(-0.9508739)])

            b1=np.array([0.1044253,-0.1-2*ny**(-0.4),-0.2-0.9*ny**(-0.4), \
                         -2.7113701*ny**(-1.000845),0.1457385, \
                          1.0286911*ny**(-0.9818393)])
            
            b2=np.array([0.7094884,0.90+1.5*ny**(-0.4),0.7+0.73*Delta**(0.4), \
                        1.251972,1.603306,1.181565])  
        
        xa=b0[0]+b1[0]*(math.log10(1+abs(m)))**b2[0]
        xb=b0[1]+b1[1]*(math.log10(1+abs(m)))**b2[1]
        ea=b0[2]+b1[2]*(math.log10(1+abs(m)))**b2[2]
        eb=b0[3]+b1[3]*(math.log10(1+abs(m)))**b2[3]
        ya=b0[4]+b1[4]*(math.log10(1+abs(m)))**b2[4]
        yb=b0[5]+b1[5]*(math.log10(1+abs(m)))**b2[5]
        
        if(xa>xb):
            xa=(xa+xb)/2
            xb=xa
            ea=(ea+eb)/2
            eb=ea
            xs=xa
     
        if(ny==4 and m==0):
            ea=-0.823333
        if(ny==5 and m==0):
            ea=-1.1
 
        
        psi=xb-xa
        es=(eb+ea)/2
        ca=(ea-en)/ya
        cs=(es-ea)*(eb-es)*psi**2
        cb=(eny-eb)/yb

        if(xb>xa):    
            xs=(xa*(es-ea)+xb*(eb-es))/(eb-ea)
#            if(math.isnan(xs)==True):            
#                print([xa,xb])
        if(x<=xa):
            e=ea+(ea-en)*(x-xa)/(ca**2+(x-xa)**2)**(1/2)
        elif(x>=xb):
            e=eb+(eny-eb)*(x-xb)/(cb**2+(x-xb)**2)**(1/2)
        else:
            e=ea+(es-ea)*((eb-es)*psi+ \
             (eb-ea)*(x-xs))**2/(cs+(eb-ea)**2*(x-xs)**2)
    
    Ef=-10**e;
    
    if(m>0):            
        Ef=Ef+4*beta*float(m)

    if(p==1):
        Ef=Ef-2.1784e-3*beta*float(m)
        
    return Ef                      
    

# Prueba de graficas:

    
#index=1

# import pandas as pd

#if (index==0):

#    import matplotlib.pyplot as plt


#    N=500;

#    lbeta=[None]*(N+1)
#    E1=[None]*(N+1)
#    E2=[None]*(N+1)
#    E3=[None]*(N+1)
#    E4=[None]*(N+1)

#    lbeta[1]=-4

#    for I in range(2,N+1,1):
    
#        lbeta[I]=lbeta[I-1]+0.016
 
#        beta=10**lbeta[I]
  
#        E1[I-1]=(Energy(beta,0,0)) 
#        E2[I-1]=(Energy(beta,-1,0)) 
#        E3[I-1]=(Energy(beta,-2,0)) 
#        E4[I-1]=(Energy(beta,-3,0)) 

#    plt.plot(lbeta,E1)
#    plt.plot(lbeta,E2)
#    plt.plot(lbeta,E3)
#    plt.plot(lbeta,E4)
#    plt.show()    


    

