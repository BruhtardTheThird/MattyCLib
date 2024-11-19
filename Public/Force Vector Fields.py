import sympy as sp
import os.path
import sys

# go up one directory level from this file's directory:
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# prepend parent directory to the system path:
sys.path.insert(0, path)

from FunctionLib.MTHFunctions import Jacobian, Div, Curl, findPotntialFnction, CountourInt
import FunctionLib.ParametricGraphing as pmg
import matplotlib.pyplot as plt
import numpy as np
x, y, z, u, v, w = sp.symbols('x,y,z,u,v,w')
r, theta = sp.symbols('r, theta')
p, phi = sp.symbols('p, phi')
pi = sp.pi
sin, cos, tan, sqrt, e= sp.sin, sp.cos, sp.tan, sp.sqrt, sp.exp

Field1 = sp.Array([
    x-4*y,
    4*x+3*y,
    0
])


Point1 = [7,0,0]

t = sp.Symbol('t')

parameterField = [
    r*sp.cos(theta),
    r*sp.sin(theta)
]

varList1 = [x,y,z]
varList2 = [x,y]


varTupleList1 = [(x,y)]

parVarList1 = [r,theta]

parBoundsTupleList1 = [(0,1),(0,2*sp.pi)]



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
Oint = CountourInt(Field1,VarList=varList2,Hard=False,ParVars=parVarList1,Eval=True,ParaField=parameterField,ParVarBounds=parBoundsTupleList1)
print(Oint)

#pmg.NumpyParaGraph(Field1,varTupleList1,parBoundsTupleList1)
#plt.show()

