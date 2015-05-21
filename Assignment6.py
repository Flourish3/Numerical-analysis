# -*- coding: utf-8 -*-
"""
Created on Thu May 21 13:41:11 2015
@author: Johan
"""
from  scipy import *
from  pylab import *

def ydot(y):
    m2=3.
    r=(y[0]**2+y[2]**2)**(3./2.)
    z=array([0.,0.,0.,0.])
    z[0]=y[1]
    z[1]=(-m2*y[0])/r
    z[2]=y[3]
    z[3]=(-m2*y[2])/r
    return z

def ynextdot(y):
    z=y    
    for i in range(10):     
        z=z+h*z    
    return z
    
    
y = array([0,1,2,0])
xplot = [0]
yplot = [2]
h=0.001

#loop
for i in range(100000):
    y=y+h*ydot(ynextdot(y))
    xplot += [y[0]]
    yplot += [y[2]]

plot(xplot,yplot)

