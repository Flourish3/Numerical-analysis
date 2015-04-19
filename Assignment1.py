# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 12:24:31 2015

@author: Patrik
"""
from  scipy import *
from  pylab import *
from scipy.optimize import fsolve

def plotIntervall(f, a, b):
    x = linspace(a,b)
    plot(x, f(x))

def bisec(f, a, b, tol):
    tolerance = 1;
    count = 0;
    xlist =[]
    if a < b and f(a) * f(b) < 0:
    
        while tolerance > tol:
            tolerance  = abs(a-b)        
            x = (a + b) / 2
            
            if f(a) * f(x) < 0:
                b = x
            else:
                a = x
            count = count + 1;
            xlist.append(x);
        return x, count
    
    else:
        print("illegal interval")
        return 0,0;
            
def f(x,*args):
    return x**2 - 2

def newton(f, x, tol):
    tolerance = 1
    h = 10**-7;
    count = 0;
    
    while(tolerance > tol):
        xp1 = x - (f(x) * h) /(f(x + h) - f(x))
        tolerance = abs(xp1 - x);
        x = xp1
        count = count + 1
        
    return x, count
    
def secant(f, x0, x1, tol):
    tolerance = 1
    count = 0;
    while  tolerance > tol :
        xt = x1 - (f(x1)*(x1 - x0))/(f(x1) - f(x0))
        x0 = x1
        x1 = xt
        
        tolerance = abs(x1 - x0)
        count = count + 1
    return x1,count
                
plotIntervall(f, -4,4)

print(bisec(f, 0.0, 2.0,10**-7))
    
print(fsolve(f,1.0,10**-7))
            
print(newton(f,1.0, 10**-7));
print(secant(f, 0.0, 1.0, 10**-7));
