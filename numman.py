# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 12:42:47 2015
@author: Johan
"""
from  scipy import *
from  pylab import *

def LUSolve(A, b):
    L = eye(len(A))
    U = zeros((len(A), len(A)))
    U[0][0] = A[0][0] / L[0][0]
    
    if L[0][0] * U[0][0] == 0:
        print("factorization not possible!")
        return empty(len(b))
    
    for j in range(1, len(A)):
        U[0][j] = A[0][j]/L[0][0] 
        L[j][0] = A[j][0]/U[0][0]
    for i in range(1,len(A)-1):
        U[i][i] = A[i][i] / L[i][i]
        tmp = 0;
        for k in range(0,i):
            tmp += L[i,k]*U[k,i]
        U[i][i] -= tmp/L[i][i]
        
        if L[i][i]*U[i][i] == 0:
            print("LU factorization is not possible")
            return empty(len(A))
    
        for j in range(i+1, len(A)):            
            tmp = 0;
            for k in range(0,i):
                tmp += L[i,k]*U[k,j]
                                
            U[i][j] = (A[i][j] - tmp)/L[i][i]
                        
            tmp = 0;
            for k in range(0,i):
                tmp += L[j,k]*U[k,i]
                
            L[j][i] = (A[j][i] - tmp)/U[i][i]
    
    U[len(U)-1][len(U)-1] = A[len(A) - 1][len(A) - 1]
         
    tmp = 0;
    for k in range(0,len(A) - 1):
        tmp += L[len(A) -1][k]*U[k][len(A)-1]

    U[len(U)-1][len(U)-1] -= tmp        

    y = empty(len(b))    
    
    for i in range(len(A)):
        tmp = 0;
        for j in range(i):
            tmp += L[i][j]*y[j];
        y[i] = (b[i] - tmp) / L[i][i]
    
    x = empty(len(b))
        
    for i in range(len(A)):
        ii = len(A) - i - 1
        
        tmp = 0;
        for j in range(len(A)):
            ij = len(A) - j - 1            
            tmp += U[ii][ij]*x[ij]
            
        x[ii] = (y[ii] - tmp) / U[ii][ii];
        
    return x
    
def jacobi(A, b, x0, tolerance):
    Q = diag(diag(A))
    Qinv = zeros((len(Q), len(Q)))
    
    for i in range(len(Q)):
        Qinv[i][i] = 1.0/Q[i][i]
                
    G = dot(Qinv, (Q - A))
    C = dot(Qinv, b)

    tol = 1;
    
    while tol > tolerance:
        
        xn = dot(G,x0) + C 
        
        for i in range(len(x0)):
            tol += abs(xn[i] - x0[i])
        print(xn)
        x0=xn   
    return xn;

def lagrange (x ,i , xm ):
#"""
#Evaluates the i-th Lagrange polynomial at x
#based on grid data xm
#"""
    n=len( xm )-1
    y=1.
    for j in range ( n+1 ):
        if i!=j:
            y*=( x-xm[j])/( xm[i]-xm[j])
    return y

def interpolation (x , xm , ym ):
    n=len( xm )-1
    
    lagrpoly = array ([ lagrange (x ,i , xm ) for i in range ( n+1 )])
    #print(lagrpoly)    
    y = dot( ym , lagrpoly )
    return y