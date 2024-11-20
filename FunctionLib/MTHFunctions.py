import sympy as sp
from sympy import diff
from sys import argv
import warnings
x, y, z, u, v, w = sp.symbols('x,y,z,u,v,w')
def InteriorProduct(Array1, Array2):
    if len(Array1) != len(Array2):
        raise ValueError('The arrays do not have the same length!')
    productExpr = 0
    for i in range(len(Array1)):
        productExpr = productExpr + Array1[i]*Array2[i]
    return productExpr
def CrossProduct(Array1, Array2):
    if len(Array1) != len(Array2):
        raise ValueError('The arrays do not have the same length!')
    if len(Array1) == 3:
        ProductList = [
            sp.trigsimp(Array1[1]*Array2[2]-Array1[2]*Array2[1]),
            sp.trigsimp(Array1[2]*Array2[0]-Array1[0]*Array2[2]),
            sp.trigsimp(Array1[0]*Array2[1]-Array1[1]*Array2[0])
            ]
        return ProductList
    if len(Array1)==2:
        ProductList = [
            0,
            0,
            sp.trigsimp(Array1[0]*Array2[1]-Array1[1]*Array2[0])
        ]
        return ProductList
    else:
        raise IndexError('Expected arrays of length 2 or 3, instead got arrays of length:',len(Array1))
def Magnitude(ArrayLike):
    MagExp = 0
    for i in ArrayLike:
        MagExp = MagExp + i**2
    ReturnExp = sp.sqrt(sp.simplify(MagExp))
    return ReturnExp
def Jacobian(x,y,z=None):
    """For a transformation T(x,y,z) = T(u,v,w), this function defines the Jacobian,
    taking x, y, and sometimes z as functions of u, v, and sometimes w."""
    matrix = sp.Matrix([[diff(x,u),diff(x,v)],[diff(y,u),diff(y, v)]])
    if z is None:
        return sp.trigsimp(matrix.det())
    else:
        matrix = matrix.col_insert(2, sp.Matrix([[diff(x,w)],[diff(y,w)]]))
        matrix = matrix.row_insert(2, sp.Matrix([[diff(z,u),diff(z,v),diff(z,w)]]))
        return sp.trigsimp(matrix.det())
def Div(Field,Point=None,VarList=None):
    """For a force vector field, F = <P,Q,R> and a point O = <x,y,z> (R&z = 0 for 2D) the divergence of F at point O is
    equal to the dot product of F and Del. Both arguments are input as lists of functions in terms of x, y, and z."""
    divExpr = sp.S.Zero
    if VarList is None:
        VarList = [x,y,z]
        warnings.warn('WARNING: Absence of VarList argument might break the output, as variables are set to [x,y,z]!',UserWarning)
    for f, s in zip(Field,VarList):
        divExpr = divExpr + diff(f,s)
    if Point == None:
        return divExpr
    else:
        for p,s in zip(Point,[x,y,z]):
            divExpr = divExpr.subs([(s,p)])
        return divExpr
def Curl(Field,Point=None,VarList=None):
    """For a force vector field, F = <P,Q,R> and a point O = <x,y,z> (R&z = 0 for 2D), the curl of F at point O is
    equal to the cross product of F and Del. Both arguments are input as lists of functions in terms of x, y, and z."""
    if VarList is None:
        VarList = [x,y,z]
        warnings.warn('WARNING: Absence of VarList argument might break the output, as variables are set to [x,y,z]!',UserWarning)
    else:
        gList = []
        for i in Field:
            iDiff = []
            for j in VarList:
                iDiff.append(sp.trigsimp(diff(i,j)))
            gList.append(iDiff)
            curlAray = sp.Array([gList])
    if len(curlAray) == 3:
        curlAray = sp.Array([sp.trigsimp(gList[2][1]-gList[1][2]),
                    sp.trigsimp(gList[0][2]-gList[2],[0]),
                    sp.trigsimp(gList[1][0]-gList[0][1])])
    else:
        curlAray = sp.Array([0,0,sp.trigsimp(gList[1][0]-gList[0][1])])
    if Point is None:
        return curlAray
    else:
        curlEval = []
        for i,j,k in zip(VarList,Point,curlAray):
            curlEval.append(k.subs(i,j))
        return curlEval
def findPotntialFnction(Field,VarList=None):
    if Curl(Field,VarList=VarList).tolist()!=[0]*len(Field):
        raise ValueError('Field:',Field,'Is not conservative. Therefore such a function does not exist.')
    if VarList is None:
        VarList = [x,y,z]
        warnings.warn('WARNING: Absence of VarList argument might break the output, as variables are set to [x,y,z]!',UserWarning)
    for f,s in zip(Field,VarList):
        for i in sp.expand(sp.integrate(f,s)).as_ordered_terms():
            if i not in VarList:
                VarList.append(i)
    potenExpr = sp.S.Zero
    for i in VarList:
        potenExpr = potenExpr + i
    return potenExpr
def ParametrizeExpr(Expression,MappedArray,VarList=None):
    """Returns an expression where the values of the mapping array have been substituted for the variables in the expression. If VarList is left blank,
    the function will substitute the values of the mapping array with the variables it 'sees' first. VarList should be input as a list of variables
    that are used in the original function."""
    if VarList != None:
        if len(MappedArray) != len(VarList):
            raise ValueError('The parametrization array:',MappedArray,'does not have the same dimension as the input:',VarList)
        for i,j in zip(VarList,MappedArray):
            Expression = Expression.subs(i,j)
        return Expression
    else:
        for i,j in zip(Expression.free_symbols,MappedArray):
            Expression = Expression.subs(i,j)
        return Expression
def CountourInt(Field,VarList,VarBounds=None,Eval=False,Hard=True,ParaField=None,ParVar=None,ParVars=None,ParVarBounds=None,):
    if Hard == True and ParVars == None:
        ParametrizedField = []
        for i in Field:
            ParametrizedField.append(i.subs(zip(VarList,ParaField)))
        integrand = 0
        for i in range(len(ParaField)):
            integrand = integrand + ParametrizedField[i]*sp.diff(ParaField[i],ParVar)
        hardIntegral = sp.Integral(integrand,(ParVar,ParVarBounds[0]))
        if Eval == False:
            print('Integrand for CounterInt =\n',integrand)
            return hardIntegral
        else:
            return hardIntegral.doit()
    else:
        if VarBounds==None:
            StokesIntegrand = (InteriorProduct(Curl(Field,VarList=VarList),[1,1,1])).subs(zip(VarList,ParaField))
            if ParaField == [sp.Symbol('r')*sp.cos(sp.Symbol('theta')),sp.Symbol('r')*sp.sin(sp.Symbol('theta'))] or ParaField == [sp.Symbol('r')*sp.cos(sp.Symbol('theta')),sp.Symbol('r')*sp.sin(sp.Symbol('theta')),sp.Symbol('z')]:
                StokesIntegrand = StokesIntegrand*sp.Symbol('r')
            StokesIntegral = sp.Integral(sp.simplify(StokesIntegrand),(ParVars[0],ParVarBounds[0]),(ParVars[1],ParVarBounds[1]))
            if Eval == False:
                print('Integrand for ContourInt =\n',StokesIntegrand)
                return StokesIntegral
            return StokesIntegral.doit()
        StokesIntegrand = InteriorProduct(Curl(Field,VarList=VarList),[1,1,1])
        StokesIntegral = sp.Integral(StokesIntegrand,(VarList[0],VarBounds[0]),(VarList[1],VarBounds[1]))
        if Eval == False:
            print('Integrand for ContourInt =\n',StokesIntegrand)
            return StokesIntegral
        return StokesIntegral.doit()
def SurfArea(Surface,VarList,ParaField=sp.Array|None,ParVarList=list|None,Eval=True,ParVarBounds=list|None,Easy=False,Left=True):
    if Easy == False:
        ParaSurf = []
        for i in Surface:
            ParaSurf.append(i.subs(zip(VarList,ParaField)))
        NormalPre = []
        for j in range(2):
            NormalPre.append([diff(ParaField[i],ParVarList[j]) for i in range(3)])
        Integrand = Magnitude(CrossProduct(NormalPre[0],NormalPre[1]))
    else:
        Integrand = sp.sqrt(diff(Surface[2],VarList[0])**2+diff(Surface[2],VarList[1])**2+1)
    if Left == True:
        Integral = sp.Integral(Integrand,(VarList[0],ParVarBounds[0]),(VarList[1],ParVarBounds[1]))
    else:
        Integral = sp.Integral(Integrand,(VarList[1],ParVarBounds[1]),(VarList[0],ParVarBounds[0]))
    if Eval == True:
        print('Surface Integral integrand is: \n',Integrand,'\nSurface Integral is:\n',Integral)
        EvalDIntegral = Integral.doit()
        return EvalDIntegral
    else:
        print('Surface Integral integrand is: \n',Integrand,'\nSurface Integral is:\n',Integral)
        return Integral

#print(Jacobian(u*sp.cos(v),u*sp.sin(v)))
#matrix1 = sp.Matrix([[diff(u*sp.sin(v)*sp.cos(w),u),diff(u*sp.sin(v)*sp.cos(w),v)],[diff(u*sp.sin(v)*sp.sin(w),u),diff(u*sp.sin(v)*sp.sin(w),v)]])
#matrix1 = matrix1.col_insert(2, sp.Matrix([[diff(u*sp.sin(v)*sp.cos(w),w)],[diff(u*sp.sin(v)*sp.sin(w),w)]]))
#print(matrix1)
#Testing with the matrix created by the function.
#print(Jacobian(u*sp.sin(v)*sp.cos(w),u*sp.sin(v)*sp.sin(w),u*sp.cos(v)))
assert Jacobian(u*sp.cos(v),u*sp.sin(v)) == u , 'The jacobian of a polar transformation should equal r = u.'
assert Jacobian(u*sp.sin(v)*sp.cos(w),u*sp.sin(v)*sp.sin(w),u*sp.cos(v)) == (u**2)*sp.sin(v), 'The jacobian of the mapping from cartesian to spherical co-ordinates should equal p^2sin(phi) = u^2sin(v)'
