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
    3*x**3,
    5*y**5,
    -2*x*y+y
    ])



Point1 = [7,0,0]

varList1 = [x,y,z]
varList2 = [x,y]

varTupleList1 = [(x,y)]

parBoundsTupleList1 = [(0,1),(0,1)]



div1 = Div(Field1)
ediv1 = Div(Field1,Point1)
print('Divergence expression is:\n',div1,'\nDivergence evaluation is:\n',ediv1)
curl1 = Curl(Field1,VarList=varList1)
ecurl1 = Curl(Field1,Point1,varList1)
print('Curl expression is:\n',curl1,'\nCurl evaluation is:\n',ecurl1)
try:
    potFun = findPotntialFnction(Field1,VarList=varList1)
    print('A potential function is: \n',potFun)
except:
    print('Field:\n',Field1,'\nIs not conservative. Therefore a potential function does not exist for the field.')
pmg.NumpyParaGraph(Field1,varTupleList1,parBoundsTupleList1)
plt.show()
