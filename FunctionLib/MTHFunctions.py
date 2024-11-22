import sympy as sp
from sympy import diff
from sys import argv
import warnings
x, y, z, u, v, w = sp.symbols('x,y,z,u,v,w')
def InteriorProduct(Array1, Array2):
    if len(Array1) != len(Array2):
        raise ValueError('The arrays do not have the same length!')
    productExpr = sum([Array1[i]*Array2[i] for i in range(len(Array1))])
    return productExpr
def CrossProduct(Array1, Array2):
    if len(Array1) != len(Array2):
        raise ValueError('Arrays must have the same length.')

    if len(Array1) == 3:  # 3D
        return [
            sp.trigsimp(Array1[2]*Array2[1] - Array1[1]*Array2[2]),
            sp.trigsimp(Array1[0]*Array2[2] - Array1[2]*Array2[0]),
            sp.trigsimp(Array1[1]*Array2[0] - Array1[0]*Array2[1])
        ]
    elif len(Array1) == 2:  # 2D
        return [0, 0, sp.trigsimp(Array1[0]*Array2[1] - Array1[1]*Array2[0])]
    else:
        raise ValueError("Only 2D or 3D vectors are supported.")
def Magnitude(ArrayLike):
    ReturnExp = sp.sqrt(sp.simplify(sum([i**2 for i in ArrayLike])))
    return ReturnExp
def Jacobian(MappedField,VarList):
    """Returns the determinant of the Jacobian matrix for a given parametrization and variables."""
    if len(VarList) != len(MappedField):
        raise ValueError('MappedField and VarList must have the same length.')
    
    # Compute the Jacobian matrix
    jacobian_matrix = sp.Matrix([[diff(f, v) for v in VarList] for f in MappedField])
    return sp.trigsimp(jacobian_matrix.det())

def Div(Field,Point=None,VarList=None):
    """For a force vector field, F = <P,Q,R> and a point O = <x,y,z> (R&z = 0 for 2D) the divergence of F at point O is
    equal to the dot product of F and Del. Both arguments are input as lists of functions in terms of x, y, and z."""
    if VarList is None:
        VarList = [x,y,z]
        warnings.warn('WARNING: Absence of VarList argument might break the output, as variables are set to [x,y,z]!',UserWarning)
    divExpr = sum([diff(f,s) for f,s in zip(Field,VarList)])
    if Point == None:
        return divExpr
    else:
        for p,s in zip(Point,[x,y,z]):
            divExpr = divExpr.subs([(s,p)])
        return divExpr
def Curl(Field,Point=None,VarList=None):
    """Compute the curl of a vector field. Field is a sympy array or list, Point and VarList are lists."""
    VarList = VarList or [x, y, z]
    gList = [[diff(f, v) for v in VarList] for f in Field]
    
    if len(gList) == 3:
        curl_array = [
            sp.trigsimp(gList[2][1] - gList[1][2]),
            sp.trigsimp(gList[0][2] - gList[2][0]),
            sp.trigsimp(gList[1][0] - gList[0][1])
        ]
    else:
        curl_array = [0, 0, sp.trigsimp(gList[1][0] - gList[0][1])]

    if Point is None:
        return curl_array
    else:
        curlEval = [(k.subs(i,j)) for i,j,k in zip(VarList,Point,curl_array)]
        return curlEval
def FindPotentialFnction(Field,VarList=None):
    if Curl(Field,VarList=VarList).tolist()!=[0]*len(Field):
        raise ValueError('Field:',Field,'Is not conservative. Therefore such a function does not exist.')
    if VarList is None:
        VarList = [x,y,z]
        warnings.warn('WARNING: Absence of VarList argument might break the output, as variables are set to [x,y,z]!',UserWarning)
    NewList = []
    for f,s in zip(Field,VarList):
        for i in sp.expand(sp.integrate(f,s)).as_ordered_terms():
            if i not in NewList:
                NewList.append(i)
    potenExpr = sum([i for i in NewList])
    return potenExpr
def ParametrizeExpr(Expression,MappedArray,VarList=None):
    """Substitute variables in Expression using MappedArray or VarList."""
    if VarList is not None:
        if len(MappedArray) != len(VarList):
            raise ValueError('MappedArray and VarList must have the same length.')
        return Expression.subs(zip(VarList, MappedArray))
    return Expression.subs(zip(Expression.free_symbols, MappedArray))
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
                StokesIntegrand *= sp.Symbol('r')
            StokesIntegral = sp.Integral(sp.simplify(StokesIntegrand),(ParVars[0],ParVarBounds[0]),(ParVars[1],ParVarBounds[1]))
            if Eval:
                return StokesIntegral.doit()
            print('Integrand for ContourInt =\n',StokesIntegrand)
            return StokesIntegral
        StokesIntegrand = InteriorProduct(Curl(Field,VarList=VarList),[1,1,1])
        StokesIntegral = sp.Integral(StokesIntegrand,(VarList[0],VarBounds[0]),(VarList[1],VarBounds[1]))
        if Eval == False:
            print('Integrand for ContourInt =\n',StokesIntegrand)
            return StokesIntegral
        return StokesIntegral.doit()
def SurfArea(Surface,VarList,ParaField=list|None,ParVarList=list|None,Eval=True,ParVarBounds=list|None,Easy=False,Left=True,IntExpr=None):
    if Easy == False:
        ParaSurf = []
        for i in Surface:
            ParaSurf.append(i.subs(zip(VarList,ParaField)))
        NormalPre = []
        for j in range(2):
            NormalPre.append([sp.trigsimp(diff(ParaSurf[i],ParVarList[j])) for i in range(3)])
        Integrand = sp.factor(Magnitude(CrossProduct(NormalPre[0],NormalPre[1])))
        VarList = ParVarList
    else:
        Integrand = sp.sqrt(diff(Surface[2],VarList[0])**2+diff(Surface[2],VarList[1])**2+1)
    if IntExpr != None:
        Integrand = Integrand*IntExpr
        if ParVarList != None:
            VarList=ParVarList
    if Eval == True:
        if Left == True:
            Integral = sp.Integral(Integrand,(VarList[0],ParVarBounds[0]),(VarList[1],ParVarBounds[1]))
        else:
            Integral = sp.Integral(Integrand,(VarList[1],ParVarBounds[1]),(VarList[0],ParVarBounds[0]))
        print('Surface Integral integrand is: \n',Integrand,'\nSurface Integral is:\n',Integral)
        EvalDIntegral = Integral.doit()
        return EvalDIntegral
    else:
        return Integrand
def SurfFlux(Surface,Field,VarList,ParaVarList,Bounds=list|None,Eval=True,Left=True,Upward=False):
    ParaField = [f.subs(zip(VarList, Surface)) for f in Field]
    NormalDiffs = [[sp.trigsimp(diff(Surface[i], ParaVarList[j])) for i in range(3)] for j in range(2)] # gets the partial derivatives of the normal vector.
    Integrand = sp.factor(InteriorProduct(CrossProduct(NormalDiffs[0],NormalDiffs[1]),ParaField))
    if Upward:
        Integrand *= -1
    if Eval:
        if Left:
            Integral = sp.Integral(Integrand,(ParaVarList[0],Bounds[0]),(ParaVarList[1],Bounds[1]))
        else:
            Integral = sp.Integral(Integrand,(ParaVarList[1],Bounds[1]),(ParaVarList[0],Bounds[0]))
        print('Surface Integral integrand is: \n',Integrand,'\nSurface Integral is:\n',Integral)
        EvalDIntegral = Integral.doit()
        return EvalDIntegral
    else:
        return Integrand

#print(Jacobian(u*sp.cos(v),u*sp.sin(v)))
#matrix1 = sp.Matrix([[diff(u*sp.sin(v)*sp.cos(w),u),diff(u*sp.sin(v)*sp.cos(w),v)],[diff(u*sp.sin(v)*sp.sin(w),u),diff(u*sp.sin(v)*sp.sin(w),v)]])
#matrix1 = matrix1.col_insert(2, sp.Matrix([[diff(u*sp.sin(v)*sp.cos(w),w)],[diff(u*sp.sin(v)*sp.sin(w),w)]]))
#print(matrix1)
#Testing with the matrix created by the function.
#print(Jacobian(u*sp.sin(v)*sp.cos(w),u*sp.sin(v)*sp.sin(w),u*sp.cos(v)))
assert Jacobian([u*sp.cos(v),u*sp.sin(v)],[u,v]) == u , 'The jacobian of a polar transformation should equal r = u.'
assert Jacobian([u*sp.sin(v)*sp.cos(w),u*sp.sin(v)*sp.sin(w),u*sp.cos(v)],[u,v,w]) == (u**2)*sp.sin(v), 'The jacobian of the mapping from cartesian to spherical co-ordinates should equal p^2sin(phi) = u^2sin(v)'
