# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 13:54:06 2015
@author: Johan
"""
from  scipy import *
from  pylab import *
from  numman import *
from scipy.linalg import solve


def w(x,n,xn):
    y=1
    for i in range(n):
            y*=(x-xn[i])
    return y
    
def f(x):
    return 1/(1+25*x**2) 
    
def T(x,n):
    return cos(n*arccos(x))
def W(x,n):
    return T(x,n)/2**(n-1)
    
days = array([1,2,3,4,5])
energy = array([27.93, 46.98, 31.95, 31.68, 21])


choice = raw_input("What task? ")

if choice == "1":
#---------------------------------------------------------------
#Task 1

    #energy = energy*1000
    #Using Polyfit/Polyval
    p=polyfit(days,energy,4)
    figure(1)
    #plot(days,energy,'*')
    x=linspace(1,5,100)
    #plot(x,polyval(p,x))
    print(p)
    print(polyval(p,9))
    
    #Using Vandermonde
    A=vander(days)
    b=energy
    a=solve(A,b)
    print(a)
    print(polyval(a,9))
    #plot(x,polyval(a,x))
    
    #Using Lagrange
    yplot=interpolation(x,days,energy)
    #plot(x,yplot)
    lag=polyfit(x,yplot,4)
    print(polyval(lag,9))
    
    #Sensitivity
    energy2 = array([27.93, 46.98, 31.95, 31.68, 21])
    
    p2=polyfit(days,energy2,4)
    
    b=energy2
    a2=solve(A,b)
    
    yplot=interpolation(x,days,energy2)
    lag2=polyfit(x,yplot,4)
    
    dayChanged = 4
    relDiffen=abs((energy2[dayChanged]-energy[dayChanged]))/abs(max(array([energy[dayChanged],energy2[dayChanged]])))
    
    relRes=abs(polyval(a2,9)-polyval(a,9))/abs(max(array([polyval(a,9),polyval(a2,9)])))
    
    print(relRes/relDiffen)
#-------------------------------------------------------------
#Task 2   
if choice == "2":
   
    n=5
    xn=linspace(-1,1,n)
    #xn=[-0.9,-0.7,-0.3,-0.2,0.5]
    x=linspace(-1,1,100)
    plot(x,w(x,n,xn))
#-------------------------------------------------------------    
#Task 3
if choice == "3":
  
    n=5
   
        
    xn=array([cos(((2.*k-1.)/(2.*n))*pi) for k in range(1,n+1)])
    print(xn)
    plot(x,w(x,n,xn))
#------------------------------------------------------------
#Task 4
if choice == "4":
        
    n=17
  
    x=linspace(-1,1,100)
    
    plot(x,f(x), label = 'f(x)')
    
    xn=array([cos(((2.*k-1.)/(2.*n))*pi) for k in range(1,n+1)])
    y=[]
    
    for i in range(n):
        y+= [f(xn[i])]
    
    p = polyfit(xn,y,n-1)
    
    plot(x,polyval(p,x), label='Chebychev')
    
    xn = linspace(-1,1,n)
    y=[]
    
    for i in range(n):
        y+= [f(xn[i])]
    
    p = polyfit(xn,y,n-1)
    
    plot(x,polyval(p,x), label='Lagrange')
    legend()
    
#-----------------------------------------------------------------
 





