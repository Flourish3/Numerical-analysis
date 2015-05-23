# -*- coding: utf-8 -*-
"""
Created on Fri May 22 09:24:30 2015
@author: Johan
"""
from  scipy import *
from  pylab import *
from  Assignment4 import cubespline, cubsplineval
from  numman import *

def rev(x,y):
    
    #Lagrange polynom    
    xplot = linspace(-1,1,100)
    yplot = interpolation(xplot,x,y)
    #print(yplot)
    coeff = cubespline(x,y)
    yplot2 = cubsplineval(coeff,x,xplot)
    
#    plot(x,y,'o')

    
    yder = derive(xplot,derive(xplot,yplot))
    yder2 = derive(xplot,derive(xplot,yplot2))
    
#    plot(xplot,yder,'r')
#    plot(xplot,yder2,'g')

    yint = integral(xplot,yder)
    yint2 = integral(xplot,yder2)
    
    print(yint)
    print(yint2)
    
def integral(x,y):    
    h=x[1]-x[0]
    integral = 0
    
    for i in range(len(x)):
        if i == 0 or i==len(x)-1:
            integral+=1/2 * y[i]
        else:
            integral +=y[i] 
    return h*integral
    
def derive(x,y):

    yder = []
    
    for i in range(len(y)-1):
        yder += [((y[i+1]-y[i])/(x[i+1]-x[i]))]
    
    yder += [yder[len(yder)-1]+(yder[len(yder)-1]-yder[len(yder)-2])]
    return yder
    

rev([-1,0,1],[1,0,1])
#-------------------------------------------------------------
def ydot(y):
    g=9.82
    z=array([0.,0.])
    z[0]=y[1]
    z[1]=-g*sin(y[0])

    return z

def ynextdot(y):
    z=y    
    for i in range(10):
        z=y+h*z    
    return z
    
def nextf(y,lasty):
    z=y    
    for i in range(10):     
        z= 4/3*z-1/3 * lasty + h*z  
    return z
    
y = array([pi/4,0])
lasty=y
xplot = [pi/4]
yplot = [0]
h=0.01

y=y+h*ydot(ynextdot(y))
 
#loop

for i in range(50):
    y= 4/3*y-1/3 * lasty + 2*h/3*nextf(y,lasty)
    xplot += [y[0]]
#    yplot += [y[2]]
for i in range(len(xplot)):
    xplot[i] -=pi/2
print(xplot)
plot(cos(xplot),sin(xplot))