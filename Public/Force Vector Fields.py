import sympy as sp
import os.path
import sys

# go up one directory level from this file's directory:
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# prepend parent directory to the system path:
sys.path.insert(0, path)


import FunctionLib.MTHFunctions as mtf
import FunctionLib.ParametricGraphing as pmg
import matplotlib.pyplot as plt
import numpy as np
x, y, z, u, v, w = sp.symbols('x,y,z,u,v,w')
a,b = sp.symbols('a,b',)
r, theta = sp.symbols('r, theta',)
p, phi = sp.symbols('p, phi')
pi = sp.pi
sin, cos, tan, sqrt, e= sp.sin, sp.cos, sp.tan, sp.sqrt, sp.exp

Field1 = sp.Array([
    x/2,
    y/2,
    0
])
Surface1 = sp.Array([
    x,
    y,
    25+x+y+x**2
])


Point1 = [7,0,0]

t = sp.Symbol('t')

parameterField1 = [
    r*sp.cos(theta),
    r*sp.sin(theta)
]


varList1 = [x,y,z]
varList2 = [x,y]


varTupleList1 = [(x,y)]

parVarList1 = [r,theta]

parBoundsTupleList1 = [(0,2),(3,9)]


div1 = mtf.Div(Field1)
ediv1 = mtf.Div(Field1,Point1)
print('Divergence expression is:\n',div1,'\nDivergence evaluation is:\n',ediv1)
curl1 = mtf.Curl(Field1,VarList=varList1)
ecurl1 = mtf.Curl(Field1,Point1,varList1)
print('Curl expression is:\n',curl1,'\nCurl evaluation is:\n',ecurl1)

#try:
#    potFun = findPotntialFnction(Field1,VarList=varList1)
#    print('A potential function is: \n',potFun)
#except:
#    print('Field:\n',Field1,'\nIs not conservative. Therefore a potential function does not exist for the field.')
#Oint = CountourInt(Field1,VarList=varList2,Hard=True,Eval=True,ParVarBounds=parBoundsTupleList1,ParVar=t,ParaField=parameterField,)
#print(Oint)
print(mtf.SurfArea(Surface1,varList2,ParVarList=parVarList1,ParaField=parameterField1,Eval=True,Easy=True,ParVarBounds=parBoundsTupleList1,Left=False))
#print(mtf.Magnitude(mtf.CrossProduct([-v*sp.sin(u),v*sp.cos(u),1],[sp.cos(u),sp.sin(u),0])))
#NormalPre = []
#for j in range(2):
#    NormalPre.append([sp.trigsimp(sp.diff([v*sp.cos(u),v*sp.sin(u),u][i],[u,v][j])) for i in range(3)])
#print(NormalPre)
#pmg.NumpyParaGraph(Field1,varTupleList1,parBoundsTupleList1)
#plt.show()

