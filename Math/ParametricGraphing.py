import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import warnings as war
u,v,w = sp.symbols('u,v,w')
def NumpyParaCurve(ListFunction,EvalTuple,VarTupleList=None):
    """Returns a list of the evaluations of a parameterized function with variables specified by VarTupleList. Otherwise uses u and v.
    Essentially, changes a sympy parameterized function into one that can be accessed by numpy, which is why EvalTuple should be in numpy syntax.
    ListFunction should be either a sympy array of expressions or a list of sympy expressions. VarTupleList should be a list of a single tuple."""
    if VarTupleList is None:
        VarTupleList = [(u,v)]
        war.warn('WARNING: Absence of VarList argument might break the output, as variables are set to [u,v]!',UserWarning)
    NumpyListFunction = []
    for i in ListFunction:
            #creates a lambidified function for each function given, in terms of the variables given
        NumpyListFunction.append(sp.lambdify(VarTupleList,i))
    EvalList = []
    for i in NumpyListFunction:
            #evaluates each function at each point specified by the eval tuple argument
        EvalList.append(i(EvalTuple))
    return EvalList
def NumpyParaGraph(ListFunction,VarTupleList,ParBoundsTupleList):
    """Returns a surface plot for matplotlib.pyplot, where the surface is that created by the vector valued function, on the bounds given.
    ListFunction should be an arrayLike, consisting of sympy expressions. VarTupleList should be a list of a single tuple. 
    ParBoundsTupleList should be a list of tuples."""
    boundsList = []
    for j in ParBoundsTupleList:
        if j[0] == j[1]:
            raise TypeError('Expected 3-dimensional parameterization, too many variables.')
        boundsList.append(np.linspace(*j,100))#creates a list of arrays of 100 values between the bounds specified by ParBoundsTupleList
    meshTuple = np.meshgrid(*tuple(boundsList)) 
        #creates a tuple that allows a numpy function's variables to be evaluated on all points specified by the bounds
    ax = plt.figure().add_subplot(111,projection='3d') 
        #graphing stuff i dont really understand, but 'creates' a 3d graphing space
    return ax.plot_surface(*tuple(NumpyParaCurve(ListFunction,meshTuple,VarTupleList)),cmap='viridis') #plots a surface where x, y, and z are all a function of the NumpyParaCurve function
