import sympy as sp
from MTHFunctions import Jacobian, Div, Curl, findPotntialFnction
import ParametricGraphing as pmg
import matplotlib.pyplot as plt
import numpy as np
x, y, z, u, v, w = sp.symbols('x,y,z,u,v,w')
r, theta = sp.symbols('r, theta')
p, phi = sp.symbols('p, phi')
pi = sp.pi
sin, cos, tan, sqrt, e= sp.sin, sp.cos, sp.tan, sp.sqrt, sp.exp

Field1 = sp.Array([
    p*sp.sin(phi)*sp.cos(theta),
    p*sp.sin(phi)*sp.sin(theta),
    p*sp.cos(phi)
    ])



Point1 = [7,0,0]

varList1 = [p,phi,theta]
varList2 = [u,v]

varTupleList1 = [(phi,theta,p)]

parBoundsTupleList1 = [(0,np.pi),(0,2*np.pi)]
#the bounds for a 3 parameter equation must exclude a parameter, otherwise the function will try to graph a 4-dimensional shape.
#also, make sure, if one parameter is constant, that it is at the end of the TupleLists, and the varLists don't necessarily reflect that.
#make sure for a 3-parameter equation, you substitute in the value of the constant parameter before executing the function. 
#This can be done by simply making a statement such that:

#InField1 = []
#Field2 = sp.Array(InField1.append(i.subs(sp.Symbol,'constant') for i in Field1))
#   and then using Field2 for the argument in NumpyParaGraph


div1 = Div(Field1)
ediv1 = Div(Field1,Point1)
print('Divergence expression is:\n',div1,'\nDivergence evaluation is:\n',ediv1)
curl1 = Curl(Field1,VarList=varList1)
ecurl1 = Curl(Field1,Point1,varList1)
print('Curl expression is:\n',curl1,'\nCurl evaluation is:\n',ecurl1)
#try:
#    potFun = findPotntialFnction(Field1,VarList=varList1)
#    print('A potential function is: \n',potFun)
#except:
#    print('Field:\n',Field1,'\nIs not conservative. Therefore a potential function does not exist for the field.')
pmg.NumpyParaGraph(Field1,varTupleList1,parBoundsTupleList1)
plt.show()