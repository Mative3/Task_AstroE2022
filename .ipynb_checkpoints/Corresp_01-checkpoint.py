#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 12:59:27 2021

@author: matias
"""

# Programa que calcula los numeros cuanticos (n,l) a partir de (m,ny).

def Corresp_01(m,ny):
    
   import math    
    
   aux=ny/2-int(ny/2)
        
   if(aux==0):    # ny par
    
        n=int(1+(2*ny)/(1+math.sqrt(2*ny+1)))+abs(m)
     
   else:
        
        n=int(2+(2*ny-2)/(1+math.sqrt(2*ny-1)))+abs(m)
   
   aux=n+m
   


   if((-1)**ny==-1):
       if(aux/2-int(aux/2)!=0):
           l=0.5*((n-abs(m))**2-1)+abs(m)-ny
       else:
           l=0.5*(n-abs(m))**2+abs(m)-ny
   elif((-1)**ny==1):
       if(aux/2-int(aux/2)!=0):
           l=0.5*((n-abs(m)+1)**2-4)+abs(m)-ny
       else:
           l=0.5*((n-abs(m)+1)**2-5)+abs(m)-ny
    
   l=int(l)
   
   return n,l

#m=0
#ny=4

#n,l=Corresp_01(m,ny)

#print(n,l)