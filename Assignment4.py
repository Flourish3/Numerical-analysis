# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 19:53:52 2015
@author: Johan
"""
from __future__ import division
from  scipy import *
from  pylab import *
import numpy as np
from matplotlib.widgets import Slider, Button, RadioButtons
from scipy.linalg import solve

from matplotlib.pyplot import *

from s1002 import *

close('all')
#---------------------------------------------------------------------------------------------
#Cubic splines

def cubespline(xint,yint):# should return mx4 matrix with coefficiants
    h=xint[1]-xint[0]
    a = len(xint)-1
    
    coeff = zeros((a,4))
    
    #Adding d to the matrix    
    for i in range(a):
        coeff[i][3] = yint[i]
    
    #To calculate sigma we need to solve diagonal 141 matrix
    A=diag((a+1)*[4]) + diag((a)*[1],-1) + diag((a)*[1],1)
    
    A[0][0]=1
    A[0][1]=0
    
    A[a][a]=1
    A[a][a-1]=0 
    
    b = [0]
    
    for i in range(1,a):
        b+= [(6/h**2)*(yint[i+1] - 2*yint[i] + yint[i-1])]
    
    b+=[0]
    
    #Solve fo the sigmas
    sigma = solve(A,b)
   
   
    for i in range(a):
         #Adding ai to the matrix
        coeff[i][0] = (1./(6.*h))*(sigma[i+1]-sigma[i])
        
        #Adding b to out matrix
        coeff[i][1] = sigma[i]/2.
        
        #Adding c to the matrix
        coeff[i][2]=(yint[i+1]-yint[i])/h - ((h/6.)*(2.*sigma[i]+sigma[i+1]))
   
    return coeff

def cubsplineval(coeff,xint,xval):
    yplot =[]
     
    for k in range(len(xval)):
        #Argmax gives false on all when k is the last element(6 in this case)        
        if k==len(xval)-1:
            i = argmax(xint>=xval[k])-1
        else:
            i = argmax(xint>xval[k]) -1
        tmp= (coeff[i][0]*(xval[k]-xint[i])**3 + coeff[i][1]*(xval[k]-xint[i])**2 + coeff[i][2]*(xval[k]-xint[i])**1 + coeff[i][3])
        
        
        yplot += [tmp]
    return yplot

x = array([0,1,2,3,4,5,6])
y = array([1,3,-2,0,1,0,1])

coeff = cubespline(x,y)
xplot = linspace(0,6,100)
yplot=cubsplineval(coeff,x,xplot)

#Plots the data
#plot(x,y,'o')

#plot(xplot,yplot)

#--------------------------------------------------------------------------------------------------

eps = .0000000000000000000001
def Bsplbasis(t,v,dt):
    #Evaluate the piecewise cubic B-spline curve at intervals of dt    
    #knots are t(i) and de-boor points v(i)

    #%  IMPORTANT: 
    #(1) there should be 4 more knots t(i) than points v(i).
    #(2) This function can be used to evaluate the B-spline curve over [t_3, t_{m-3}],
    #       where m is the number of knots. (Note the function is not 
    #       well-defined outside this interval.)
    m = len(t) #number of knots
    i = 4      #index of first knot
    q=[]
    for u in arange(t[3],t[m-3]+dt,dt):
        
        # check if u value has moved to the next knot interval
        # include small tolerance on knots to avoid round-off error in comparisons.
        
        while (u>(t[i]+eps)):
            i+=1
        # Now evaluate the spline at u using the deBoor algorithm.
        # Start with the relevant control points.
        # w used here to simplify indices.
        w = i-4
        qq = zeros(len(v))
        for j in arange(1,5,1):
            
            qq[j-1]=v[w+j-1]
        for j in arange(1,4,1):
            for k in arange(1,4-j+1,1):
                qq[k-1] = ((t[w + k + 4-1] - u)/(t[w + k + 4-1] - t[w + k + j-1])) * qq[k-1] + ((u - t[w + k + j-1])/(t[w + k + 4-1] - t[w + k + j-1])) * qq[k+1-1]
        #Create vector of points on the B-spline curve.
        q.append(qq[0])
    return q

#------------------------------------------------------------------------------
#Slider plot

#fig, ax = plt.subplots()
#plt.subplots_adjust(left=0.25, bottom=0.52)
#
#plt.plot(x,y,'o')
#
#init = 2
#min = -15.
#max = 15.
#
#x =      array([0,0,0,0,1,2,3,4,5,6,6,6,6]) #Same x values but with extra data on the side
#deBoor = array([0,0,0,0,0,0,0,0,0])         #Initial deBoor points
#
#y=Bsplbasis(x,deBoor,(x[len(x)-1]-x[0])/99.)
#plt.plot(xplot,yplot)
#l, = plt.plot(xplot,y, lw=2, color='red')
#plt.axis([0, 6, min, max])
#
##Adds axis and sliders
#axcolor = 'lightgoldenrodyellow'
#v0ax = plt.axes([0.25, 0.05, 0.65, 0.03], axisbg=axcolor)
#v1ax = plt.axes([0.25, 0.1,  0.65, 0.03], axisbg=axcolor)
#v2ax = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
#v3ax = plt.axes([0.25, 0.2,  0.65, 0.03], axisbg=axcolor)
#v4ax = plt.axes([0.25, 0.25, 0.65, 0.03], axisbg=axcolor)
#v5ax = plt.axes([0.25, 0.3,  0.65, 0.03], axisbg=axcolor)
#v6ax = plt.axes([0.25, 0.35, 0.65, 0.03], axisbg=axcolor)
#v7ax = plt.axes([0.25, 0.4,  0.65, 0.03], axisbg=axcolor)
#v8ax = plt.axes([0.25, 0.45, 0.65, 0.03], axisbg=axcolor)
#
#sv0 = Slider(v0ax, 'v0', min, max, valinit=init)
#sv1 = Slider(v1ax, 'v1', min, max, valinit=init)
#sv2 = Slider(v2ax, 'v2', min, max, valinit=init)
#sv3 = Slider(v3ax, 'v3', min, max, valinit=init)
#sv4 = Slider(v4ax, 'v4', min, max, valinit=init)
#sv5 = Slider(v5ax, 'v5', min, max, valinit=init)
#sv6 = Slider(v6ax, 'v6', min, max, valinit=init)
#sv7 = Slider(v7ax, 'v7', min, max, valinit=init)
#sv8 = Slider(v8ax, 'v8', min, max, valinit=init)
#
##Adds function that updates the values v0-v8 when you change the slider
#def update(val):
#    v0=sv0.val
#    v1=sv1.val
#    v2=sv2.val
#    v3=sv3.val
#    v4=sv4.val
#    v5=sv5.val
#    v6=sv6.val
#    v7=sv7.val
#    v8=sv8.val
#    deBoor=[v0,v1,v2,v3,v4,v5,v6,v7,v8]
#    
#    l.set_ydata(Bsplbasis(x,deBoor,(x[len(x)-1]-x[0])/99.))
#    fig.canvas.draw_idle()
#
#sv0.on_changed(update)
#sv1.on_changed(update)
#sv2.on_changed(update)
#sv3.on_changed(update)
#sv4.on_changed(update)
#sv5.on_changed(update)
#sv6.on_changed(update)
#sv7.on_changed(update)
#sv8.on_changed(update)
#
##Adds a reset button
#resetax = plt.axes([0.8, 0.005, 0.1, 0.04])
#button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
#def reset(event):
#    sv0.reset()
#    sv1.reset()
#    sv2.reset()
#    sv3.reset()
#    sv4.reset()
#    sv5.reset()
#    sv6.reset()
#    sv7.reset()
#    sv8.reset()
#button.on_clicked(reset)
#
##Adds a save values button
#saveax = plt.axes([0.25, 0.005, 0.1, 0.04])
#savebutton = Button(saveax,'Save', color=axcolor, hovercolor='0.975')
#def save(event):
#    v0=sv0.val
#    v1=sv1.val
#    v2=sv2.val
#    v3=sv3.val
#    v4=sv4.val
#    v5=sv5.val
#    v6=sv6.val
#    v7=sv7.val
#    v8=sv8.val
#    deBoor=[v0,v1,v2,v3,v4,v5,v6,v7,v8]
#    print(deBoor)
#savebutton.on_clicked(save)
#
##Changes the coulor of the line
#rax = plt.axes([0.025, 0.5, 0.15, 0.15], axisbg=axcolor)
#radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)
#def colorfunc(label):
#    l.set_color(label)
#    fig.canvas.draw_idle()
#radio.on_clicked(colorfunc)
#
#plt.show()

#-------------------------------------------------------------------------------------
#Task 3 - train data
#wheel = []
#mm = [i for i in range(-70,60+1)]
#
#for i in range(len(mm)):
#    wheel += [(-1.)*s1002(mm[i])] #Here I should not have to use the -1 to det it looking right?
#
#coeff = cubespline(mm,wheel)
#xplot = linspace(mm[0],mm[len(mm)-1],100)
#yplot=cubsplineval(coeff,mm,xplot)
#
##Plots the data
##plot(mm,wheel,'o')
#
#plot(xplot,yplot)

















