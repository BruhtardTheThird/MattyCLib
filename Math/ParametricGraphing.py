import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import warnings as war
u,v,w = sp.symbols('u,v,w')
def NumpyParaCurve(ListFunction,EvalTuple,VarTupleList=None):
    if VarTupleList is None:
        VarTupleList = [(u,v)]
        war.warn('WARNING: Absence of VarList argument might break the output, as variables are set to [u,v]!',UserWarning)
    NumpyListFunction = []
    for i in ListFunction:
        NumpyListFunction.append(sp.lambdify(VarTupleList,i))
    EvalList = []
    for i in NumpyListFunction:
        EvalList.append(i(EvalTuple))
    return EvalList
def NumpyParaGraph(ListFunction,VarTupleList,ParBoundsTupleList):
    boundsList = []
    for j in ParBoundsTupleList:
        boundsList.append(np.linspace(*j,100))
    meshTuple = np.meshgrid(tuple(boundsList))
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    return ax.plot_surface(tuple(NumpyParaCurve(ListFunction,meshTuple,VarTupleList)),cmap='viridis')

