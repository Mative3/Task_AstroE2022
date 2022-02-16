#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 18:53:54 2021

@author: matias
"""

# Programa que calcula los numeros cuanticos (m,ny) a partir de (n,l,m)

def Corresp_02(n,l,m):
    
    piz=(-1)**(l+m)
    
    aux=n+m
    
    if(piz>0):
        
        if(aux/2-int(aux/2)==0):
            
            ny=0.5*((n-abs(m)+1)**2-5)-l+abs(m)
            
        else:

            ny=0.5*((n-abs(m)+1)**2-4)-l+abs(m)
        
    else:
        
        if(aux/2-int(aux/2)==0):

            ny=0.5*(n-abs(m))**2-l+abs(m)
            
        else:
            
            ny=0.5*((n-abs(m))**2-1)-l+abs(m)

    ny=int(ny)            
    
    return ny            