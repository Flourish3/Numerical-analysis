# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:20:56 2015

@author: Patrik
"""
from  scipy import *
from  pylab import *
from scipy.optimize import fsolve
from scipy.linalg import solve


def jacobian (f , x , eps = 1.e-8 , fx = []):
    """
    Computes the Jacobian of f
    eps is the increment and fx a reference value f(x)
    """
    if fx==[]: fx=f(x)
    jac = empty (( len( fx ) ,len( x )))
    for i in range (len( x )):
        x[i]+=eps
        jac [:,i]=( f(x) - fx )/eps
        x[i]-=eps
    return jac
    
    
def F(x,*args):
    y1 = 5*cos(x[0]) + 6*cos(x[0] + x[1]) - 10
    y2 = 5*sin(x[0]) + 6*sin(x[0] + x[1]) -4
    return array([y1,y2])
    
y = fsolve(F, array([.7,.7]),10**-7)

#print(F([y[0],y[1]]))

def newton(f, x0, tolerance):    
    tol = 1;
    xn = x0
    while tol > tolerance:
        xn1 = xn + xDelta(f, xn)
        tol = secondNormDiff(xn1,xn)
        xn = xn1;
    return xn;
    
def secondNormDiff(x, y):
    tmp = 0;
    for i in range(len(x)):
        tmp = (x[i] - y[i])**2        
    return sqrt(tmp)
    
def xDelta(f, xn):
   #return solve(jacobian(f,xn),-f(xn))
   
    return LUSolve(jacobian(f,xn),-f(xn))
    
#y = newton(F, array([.7,.7]), 10**-7)

#print(y)
#print(F([y[0],y[1]]))

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
    
A = array([ [1,2,3],
            [3,4,1],
            [1,0,1]])
b = array([0,6,1])

#print(LUSolve(A,b))

a = 1;
A = diag(30*[-(2+a)]) +  diag(29*[1],-1) +  diag(29*[1],1) 

b = zeros(30)
b[5] = 2;

xexact = LUSolve(A,b)

#print(xexact)

#print(dot(A,x))


def jacobi(A, b, x0, tolerance, xexact):
    Q = diag(diag(A))
    Qinv = zeros((len(Q), len(Q)))
    
    for i in range(len(Q)):
        Qinv[i][i] = 1.0/Q[i][i]
                
    G = dot(Qinv, (Q - A))
    C = dot(Qinv, b)

    tol = 1;
    xn = x0;
    
    while tol > tolerance:
        xn = dot(G,xn) + C 
        tol = 0;
        for i in range(len(xexact)):
            tol += abs(xn[i] - xexact[i])
            
    return xn;
    
def gausseidel(A, b, x0, tolerance, xexact):
    Q = tril(A)
    Qinv = zeros((len(Q), len(Q)))
                
    G = dot(Qinv, (Q - A))
    C = dot(Qinv, b)

    tol = 1;
    xn = x0;
    
    while tol > tolerance:
        xn = dot(G,xn) + C 
        tol = 0;
        for i in range(len(xexact)):
            tol += abs(xn[i] - xexact[i])
            
    return xn;

x0 = zeros(30)
print("jacobi")
#print(jacobi(A,b,x0,10**-7, xexact))    

    