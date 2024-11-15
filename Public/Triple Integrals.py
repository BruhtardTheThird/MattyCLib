import sympy as sp
import os
import sys

# go up one directory level from this file's directory:
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# prepend parent directory to the system path:
sys.path.insert(0, path)

from FunctionLib.MTHFunctions import Jacobian
x, y, z, u, v, w = sp.symbols('x,y,z,u,v,w')
r, theta = sp.symbols('r, theta')
p, phi = sp.symbols('p, phi')
pi = sp.pi

sin, cos, tan, sqrt, e= sp.sin, sp.cos, sp.tan, sp.sqrt, sp.exp


xBounds = 0,1-y
yBounds = 0,sp.sqrt(1-z)
zBounds = 0,25-r**2

uBounds = 0,1
vBounds = 0,pi/3
wBounds = 0,2*pi

rBounds = 0,1
pBounds = 0,4*cos(phi)

thetaBounds = 0,2*pi
phiBounds = 0,pi/2

f1 = x
f2 = 9*r
f3 = 1/(p**2)
f4 = 1

TX = u*sin(v)*cos(w) #u = rho, v = phi, w = theta
TY = u*sin(v)*sin(w)
TZ = u*cos(v)

#integral1 = sp.integrate(f1,(x,xBounds),(y,yBounds),(z,zBounds))
#print(integral1)
#print(integral1.evalf(5))

#integral2 = sp.integrate(r*f2,(z,zBounds),(r,rBounds),(theta, thetaBounds))
#print(integral2)
#print(integral2.evalf(5))

integral3 = sp.integrate(((p**2)*sp.sin(phi))*f3,(p,pBounds),(phi,phiBounds),(theta,thetaBounds))
print(integral3)
print(integral3.evalf(5))

#integral4 = sp.integrate((Jacobian(TX,TY,TZ))*f4,(u,uBounds),(w,wBounds),(v,vBounds))
#print(integral4)
