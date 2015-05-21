# -*- coding: utf-8 -*-
"""
Created on Wed May 13 09:21:49 2015
@author: Johan
"""
from  scipy import *
from  pylab import *
from  scipy.interpolate import splrep, splev
from  scipy.interpolate import PPoly
 
def compSimpson(f,a,b,n):
    h=(b-a)/n
    integral = f(a) + 4*f(a+h/2) + f(b)
    sum = 0
    xj=a
    for i in range(1,n):
        xj+=h
        sum+=2*f(xj)+4*f(xj+h/2)
    integral+=sum
    integral*=(h/6)
#    (f(x+h)-f(x))/h
    
#    g = derivative(f,z,n=4)    
    return integral


def f(x):
    return x
def g(x):
    return 4./(1.+x**2)
    
print(compSimpson(g,0.,1.,5))

tolerance = 10**-6
tol = 1
n=1

while tol>tolerance:
  
    tol = abs(pi-compSimpson(g,0.,1.,n))
    
    n+=1
#print(n)

def read_acc_data(file='garden.txt'):
    r"""
    reads coordinates from a file (default: garden.txt)
    and returns three lists:
    
    t time  (sec)
    
    x acceleration in x-direction $(m/sec^2)$
    
    y acceleration in y-direction $(m/sec^2)$
    
    
    Call
    
    t,xacc,yacc=read_acc_data()
    """
    fi=open(file,'r')
    x,y,t=[],[],[0]
    for line in fi.readlines():
        if line.startswith('#') or line=='\n': continue
        xi,yi,zi,ti=line.split()
        x.append(float(xi))
        y.append(float(yi))
        t.append(t[-1]+float(ti)/1000.)  # in seconds
    t=t[:-1]
    return t,x,y

t,x,y = read_acc_data()

tck = splrep(t,x)
t = tck[0]
c = tck[1]

print(t)
print(c)

a = PPoly
a.from_spline(tck)

print(a.__call__(array([1])))
t = tck[0]
c = tck[1]
plot(t, c)


