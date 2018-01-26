import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from sympy import *
from numpy import array
from sympy import log
from sympy.utilities.codegen import (InputArgument, OutputArgument,InOutArgument)
from sympy.printing import fcode
To,Fo = symbols('To,Fo')

def creepstrain(P,T1,F1,stress1,Ecp1,Ecs1,T2,F2,stress2,GetEcs,GetEcp,NumSteps):
# Integrates from F1 to F2 using NumSteps equally spaced dF steps for constant stress states i.e. stress1=stress2
    F = np.linspace(F1,F2,NumSteps+1)
    dF = (F2-F1)/NumSteps
    CS = []
    CP = []
    for i in range(NumSteps):
        Ecs2,dEcsdS = GetEcs(T1,0.,stress1,Ecs1,T2,dF,stress2,P)
        Ecp2,dEcpdS = GetEcp(T1,0.,stress1,Ecp1,T2,dF,stress2,P)
        
        if i == 0:
            Ecs2 = Ecs2/2
            dEcsdS = dEcsdS/2
            Ecp2 = Ecp2/2
            dEcpdS = dEcpdS/2
        
        Ecp1 = Ecp2
        Ecs1 = Ecs2
        
        CS.append(Ecs2[0])
        CP.append(Ecp2[0])
        
    return (F[0:-1]+F[1:])/2.,np.array(CP),np.array(CS)

def WRITE_FORTRAN_FUNCTION(NAMEFUNCTION,EXPRESSION,ARUGMENT_SEQUENCE):
    from sympy.utilities.codegen import codegen
    from sympy.parsing.sympy_parser import parse_expr
    from sympy import symbols
    exec(ARUGMENT_SEQUENCE+'=symbols("'+ARUGMENT_SEQUENCE+'")')
    result = codegen([NAMEFUNCTION, EXPRESSION], 'F95',header=False,argument_sequence=eval('['+ARUGMENT_SEQUENCE+']'))
    return result[0][1]

def WRITE_FORTRAN_SUBROUTINE(NAMESUBROUTINE,EXPRESSION,ARUGMENT_SEQUENCE):
    from sympy.utilities.codegen import codegen
    from sympy.parsing.sympy_parser import parse_expr
    from sympy import symbols
    exec(ARUGMENT_SEQUENCE+'=symbols("'+ARUGMENT_SEQUENCE+'")')
    result = codegen([NAMESUBROUTINE, EXPRESSION], 'F95',header=False,argument_sequence=eval('['+ARUGMENT_SEQUENCE+']'))
    return result[0][1]

def FORTRAN_FUNC(fid,name_subroutine,expression,P,Tmin,Tmax,Fmin,Fmax):
    from sympy.utilities.codegen import codegen
    from sympy.parsing.sympy_parser import parse_expr
    OUTPUT, To, Fo, T, F = symbols('OUTPUT To Fo T F')
    
    expression = str(eval(expression))
    #normalizeT = str(eval('(To-Tmin)/(Tmax-Tmin)'))
    #normalizeF = str(eval('(Fo-Fmin)/(Fmax-Fmin)'))    
    #expression = str(eval(expression).subs({T:normalizeT,F:normalizeF}))
    
    return expression


def FORTRAN_SUB(fid,name_subroutine,expression,P,Tmin,Tmax,Fmin,Fmax):
    from sympy.utilities.codegen import codegen
    from sympy.parsing.sympy_parser import parse_expr
    from sympy import Max
    OUTPUT, To, Fo, T, F = symbols('OUTPUT To Fo T F')
    
    normalizeT = str(eval('(To-Tmin)/(Tmax-Tmin)'))
    normalizeF = str(eval('(Fo-Fmin)/(Fmax-Fmin)'))    
    expression = str(eval(expression).subs({T:normalizeT,F:normalizeF}))
    
#    expr = Eq(OUTPUT, parse_expr(expression))
#    result = codegen([name_subroutine, [expr]], 'F95',header=False,argument_sequence=[To,Fo,OUTPUT])
    return expression

def symbols_from_dictionary(md):
# Construct Python Symbolic Symbols from Dictionary Keywords
    string1 = ''

    count = 1
    for key in md.keys():
        if count == 1:
            string1 = string1+str(key)
            count = 2
        else:
            string1 = string1+','+str(key)

    fullstring = string1+'='+'symbols("'+string1+'")'
    return fullstring

def loadColumn(filename,sheet,column):
# Load from FILENAME on an Excel Worksheet a specified Column of Data in addition to the Temperature and Fluence
# O1,O2,O3 = LOADCOLUMN(I1,I2,I3)
# 
# O1 Temperature
# O2 Fluence
# O3 Column specified by I3
#
# I1 Filename of the Excel Workbook
# I2 Name of the Worksheet
# I3 Name of the Column to be Loaded

    xls = pd.ExcelFile(filename)
    DATA = xls.parse(sheet)
    T = DATA['Temperature(C)']
    F = DATA['Fluence(dpa)']
    Y = DATA[column]
    if not(column == 'Orientation(withregardtoextrudingdirection)'):
        isData = ~np.isnan(Y)
        return np.array(T[isData],dtype=np.float64),np.array(F[isData],dtype=np.float64),np.array(Y[isData],dtype=np.float64)
    else:
        return np.array(T),np.array(F),Y

def residual1D(FORM,PARAMETERS,XDATA,YDATA):
    return FORM(XDATA,PARAMETERS) - YDATA    
    
def residual(FORM,PARAMETERS,XDATA,YDATA):
    X1 = XDATA[0]
    X2 = XDATA[1]
    return FORM(X1,X2,PARAMETERS) - YDATA

def normalize(X,Xmin,Xmax):
    return (X-Xmin)/(Xmax-Xmin)
    

def original(X,Xmin,Xmax):
    return X*(Xmax-Xmin)+Xmin

def normalize1(X):
    Xmin = X.min()
    Xmax = X.max()
    return (X-Xmin)/(Xmax-Xmin)+1,Xmin,Xmax

def fit1D(t,Y,FORM,PARAMETER0):

    def R(PARAMETERS,XDATA,YDATA):
        return residual1D(FORM,PARAMETERS,XDATA,YDATA)
    res_lsq = least_squares(R, PARAMETER0,loss='soft_l1', f_scale=0.1, args=(t,Y),method='trf', ftol=1e-08, xtol=1e-08, gtol=1e-08,max_nfev=10000)
    return res_lsq


def fit(t,f,Y,FORM,PARAMETER0):
#    t = normalize(t,t.min(),t.max())
#    f = normalize(f,f.min(),f.max())

    def R(PARAMETERS,XDATA,YDATA):
        return residual(FORM,PARAMETERS,XDATA,YDATA)
    res_lsq = least_squares(R, PARAMETER0,loss='soft_l1', f_scale=0.1, args=([t,f],Y),method='trf', ftol=1e-04, xtol=1e-04, gtol=1e-04,max_nfev=10000)
    return res_lsq
    
def axis3d(number):
	fig = plt.figure(number)
	ax = fig.add_subplot(111, projection='3d')
	return ax
	
def plot_fit(EVALF,T,F,ax):
    NT = 11
    NF = 11
    TT,FF = np.meshgrid(np.linspace(T.min(),T.max(),NT),np.linspace(F.min(),F.max(),NF))

    X1 = np.reshape(TT,(NT*NF))
    X2 = np.reshape(FF,(NT*NF))
    YFUNC = np.reshape(EVALF(X1,X2),(NT,NF))
    
    ax.plot_surface(TT,FF,YFUNC)
    
    return ax

def  evalsurface(EVALF,T,F):
    NT = 11
    NF = 11
    TT,FF = np.meshgrid(np.linspace(T.min(),T.max(),NT),np.linspace(F.min(),F.max(),NF))

    X1 = np.reshape(TT,(NT*NF))
    X2 = np.reshape(FF,(NT*NF))
    YFUNC = np.reshape(EVALF(X1,X2),(NT,NF))
    return TT,FF,YFUNC

def setconstantsub(Y,Name):
# O1,O2 = usnc.setconstantsub(I1,I2)
# Response (I2) is set as a constant function
#
# O1. Generated Python Function (T,F)
# O2. Generated Fortran Code with name specified by I4
#
# I1. Temperature,
# I2. Fluence,
# I3. Response to be fitted,
# I4. Name of the Fortran Subroutine to be generated as O3
# I5. Form of the function to be fitted with I6 number of parameters
# I6. Number of parameters in the form specified by I5

    CODE = (
'      subroutine {}(Value)\n'.format(Name)+
'!*************************************************************************\n'+
'! This subroutine constructs a subroutine that returns a constant\n'+
'! Output : Warning Message Printed to the Screen \n'+
'      implicit none\n'+
'      real*8 Value\n'+
'      Value = {}\n'.format(fcode(Y).strip())+
'      return\n'+
'      end\n')
   
    return CODE

def fitlogmodel(t,f,Y,Name,form,NumParams):
    T,F = symbols('T,F')  
    Tmin = t.min()
    Tmax = t.max()
    Fmin = 0.0
    Fmax = f.max()
    
    To = (np.array(t)-Tmin)/(Tmax-Tmin)
    Fo = (np.array(f)-Fmin)/(Fmax-Fmin)
    Yo = np.log10(np.array(Y))
    
    def EVAL(T,F,P):
        Expression = eval(form)
        return Expression

    result = fit(To,Fo,Yo,EVAL,[0.1]*NumParams)
    P = result['x']
    
    def MAKE_PYTHONFUNC(form,P0,Tmin,Tmax,Fmin,Fmax):
        def EVAL(To,Fo):
            T = (To-Tmin)/(Tmax-Tmin)
            F = (Fo-Fmin)/(Fmax-Fmin)
            P = P0
            Expression = eval(form)
            return 10**(Expression)
        return EVAL
    
    powerform = '10**'+'('+form+')'
    
    code = Get_Code(Name,fcode(eval(powerform),assign_to="Answer"),Tmin,Tmax,Fmin,Fmax)
    
    return MAKE_PYTHONFUNC(form,np.copy(P),Tmin,Tmax,Fmin,Fmax),result['x'],code

def fit1Dmodel(t,Y,Name,form,NumParams):
    T,F = symbols('T,F')  
    Tmin = t.min()
    Tmax = t.max()
    
    To = (np.array(t)-Tmin)/(Tmax-Tmin)
    Yo = np.array(Y)
    
    def EVAL(T,P):
        Expression = eval(form)
        return Expression

    result = fit1D(To,Yo,EVAL,[0.1]*NumParams)
    P = result['x']
    
    def MAKE_PYTHONFUNC(form,P0,Tmin,Tmax):
        def EVAL(To,Fo):
            T = (To-Tmin)/(Tmax-Tmin)
            P = P0
            Expression = eval(form)
            return Expression
        return EVAL
  
    code = Get_Code1D(Name,fcode(eval(form),assign_to="Answer"),Tmin,Tmax)

    
    return MAKE_PYTHONFUNC(form,np.copy(P),Tmin,Tmax),result['x'],code

def fit1Dlogmodel(t,Y,Name,form,NumParams):
    T,F = symbols('T,F')  
    Tmin = t.min()
    Tmax = t.max()
    
    To = (np.array(t)-Tmin)/(Tmax-Tmin)
    Yo = np.log10(np.array(Y))
    
    def EVAL(T,P):
        Expression = eval(form)
        return Expression

    result = fit1D(To,Yo,EVAL,[0.1]*NumParams)
    P = result['x']
    
    def MAKE_PYTHONFUNC(form,P0,Tmin,Tmax):
        def EVAL(To):
            T = (To-Tmin)/(Tmax-Tmin)
            P = P0
            Expression = eval(form)
            return 10**(Expression)
        return EVAL
    powerform = '10**'+'('+form+')'
    code = Get_Code1D(Name,fcode(eval(powerform),assign_to="Answer"),Tmin,Tmax)

    
    return MAKE_PYTHONFUNC(form,np.copy(P),Tmin,Tmax),result['x'],code




def fitmodel(t,f,Y,Name,form,NumParams):
    T,F = symbols('T,F')  
    Tmin = t.min()
    Tmax = t.max()
    Fmin = 0.0
    Fmax = f.max()
    
    To = (np.array(t)-Tmin)/(Tmax-Tmin)
    Fo = (np.array(f)-Fmin)/(Fmax-Fmin)
    Yo = np.array(Y)
    
    def EVAL(T,F,P):
        Expression = eval(form)
        return Expression

    result = fit(To,Fo,Yo,EVAL,[0.1]*NumParams)
    P = result['x']
    
    def MAKE_PYTHONFUNC(form,P0,Tmin,Tmax,Fmin,Fmax):
        def EVAL(To,Fo):
            T = (To-Tmin)/(Tmax-Tmin)
            F = (Fo-Fmin)/(Fmax-Fmin)
            P = P0
            Expression = eval(form)
            return Expression
        return EVAL
    
    code = Get_Code(Name,fcode(eval(form),assign_to="Answer"),Tmin,Tmax,Fmin,Fmax)

    
    return MAKE_PYTHONFUNC(form,np.copy(P),Tmin,Tmax,Fmin,Fmax),result['x'],code

def GetGet_invDelCode(Dir1,Dir2,Dir3,v1,v2,v3,G1,G2,G3):
    CODE = (
'      subroutine Get_invDel(T,F,iDel)\n'+
'!*************************************************************************\n'+
'! This subroutine computes the inverse elasticity tensor given the current\n'+
'! temperature and fluence\n'+
'! Input  : T    : Temperature\n'+
'!          F    : Fluence\n'+
'! Output : iDel : inverse elasticity tensor\n'+ 
'      implicit none\n'+
'      real*8 iDel(6,6)\n'+
'      real*8 T,F,c0,c2,c3,T0,b1,b2,b3,b4,b5,b6,b7,b8,b9\n'+
'      real*8 E10,E20,E30,nu12,nu23,nu31,G120,G230,G310\n'+
'      real*8 a10,a12,a22,a23,a20,a21,k1,EoverE0per,EoverE0par\n'+
'      real*8 Eper,Epar\n'+        
'\n'+
'      call GetE0par(Epar)\n'+
'      call GetEoverE0par(T,F,EoverE0par)\n'+
'      call GetE0per(Eper)\n'+
'      call GetEoverE0per(T,F,EoverE0per)\n'+
'      E10 = EoverE0{}*E{}\n'.format(Dir1,Dir1)+
'      E20 = EoverE0{}*E{}\n'.format(Dir2,Dir2)+
'      E30 = EoverE0{}*E{}\n'.format(Dir3,Dir3)+
'      nu12 = {}\n'.format(fcode(v1).strip())+
'      nu23 = {}\n'.format(fcode(v2).strip())+
'      nu31 = {}\n'.format(fcode(v3).strip())+
'      G120 = (EoverE0par+EoverE0per)/2.d0*{}\n'.format(fcode(G1).strip())+
'      G230 = (EoverE0par+EoverE0per)/2.d0*{}\n'.format(fcode(G2).strip())+
'      G310 = (EoverE0par+EoverE0per)/2.d0*{}\n'.format(fcode(G3).strip())+
'\n'+
'      iDel = 0.d0\n'+
'      iDel(1,1) = 1.d0/E10\n'+
'      iDel(2,2) = 1.d0/E20\n'+
'      iDel(3,3) = 1.d0/E30\n'+
'      iDel(4,4) = 1.d0/G120\n'+
'      iDel(5,5) = 1.d0/G310\n'+
'      iDel(6,6) = 1.d0/G230\n'+
'\n'+      
'      iDel(1,2) = -nu12/E20\n'+
'      iDel(1,3) = -nu31/E30\n'+
'      iDel(2,1) = -nu12/E10\n'+
'      iDel(2,3) = -nu31/E30\n'+
'      iDel(3,1) = -nu12/E10\n'+
'      iDel(3,2) = -nu31/E20\n'+
'\n'+
'      return\n'+
'      end\n')
    return CODE

def GetWarningCode(Tlow,Thigh,Flow,Fhigh):
    CODE = (
'      subroutine WarningCode(T,F)\n'+
'!*************************************************************************\n'+
'! This subroutine warns the user when extrapolation in T or F occurs\n'+
'! Input  : T   : Temperature\n'+
'!          F   : Fluence\n'+
'! Output : Warning Message Printed to the Screen \n'+
'      implicit none\n'+
'      real*8 T,F\n'+
'!      if (T.lt.{}) then \n'.format(fcode(Tlow).strip())+
"!          write(*,*) '*WARNING EXTRAPOLOTION* T BELOW Calibration Data'\n"+
"!          write(*,*)  'Temperature=', T\n"+
'!      endif\n'+
'!      if (T.gt.{}) then \n'.format(fcode(Thigh).strip())+
"!          write(*,*) '*WARNING EXTRAPOLOTION* T ABOVE Calibration Data'\n"+
"!          write(*,*)  'Temperature=', T\n"+
'!      endif\n'+
'!      if (F.lt.{}) then \n'.format(fcode(Flow).strip())+
"!          write(*,*) '*WARNING EXTRAPOLOTION* F BELOW Calibration Data'\n"+
"!	       write(*,*)  'Fluence=', F\n"+
'!      endif\n'+
'!      if (F.gt.{}) then \n'.format(fcode(Fhigh).strip())+
"!          write(*,*) '*WARNING EXTRAPOLOTION* F ABOVE Calibration Data'\n"+
"!	       write(*,*)  'Fluence=', F\n"+
'!      endif\n'+
'      return\n'+
'      end\n')
   
    return CODE


def GetGetTempPosTimeCode(Expression):
    CODE = (
'      subroutine GetTempPosTime(T,Coords,time)\n'+
'!************************************************************************\n'+
'! This subroutine computes the temperature as a function of time and position\n'+
'! Input  : Coords   : (x,y,z) coordinate in space\n'+
'!        : time     : time\n'+
'! Output : T        : Temperature\n'+
'\n'+
'      implicit none\n'+
'      real*8 T,Coords(3),time,X,Y,Z\n'+
'\n'+
'      X = Coords(1)\n'+
'      Y = Coords(2)\n'+
'      Z = Coords(3)\n'+        
'{}\n'.format(Expression)+
'      return\n'+
'      end\n'+
'!*************************************************************************\n')
    return CODE

def GetGetFluencePosTimeCode(Expression):
    CODE = (
'      subroutine GetFluencePosTime(F,Coords,time)\n'+
'!************************************************************************\n'+
'! This subroutine computes the temperature as a function of time and position\n'+
'! Input  : Coords   : (x,y,z) coordinate in space\n'+
'!        : time     : time\n'+
'! Output : F        : Fluence\n'+
'\n'+
'      implicit none\n'+
'      real*8 F,Coords(3),time,X,Y,Z\n'+
'\n'+
'      X = Coords(1)\n'+
'      Y = Coords(2)\n'+
'      Z = Coords(3)\n'+        
'{}\n'.format(Expression)+
'      return\n'+
'      end\n'+
'!*************************************************************************\n')
    return CODE


def Get_EwCODE(Term1,Term2,Term3):
    CODE = (
'              subroutine Get_Ew(T,F,Ew)\n'+
'!*************************************************************************\n'+
'! This subroutine computes the Wigner strain at the end of the time step\n'+
'! Input  : T   : Temperature\n'+
'!          F   : Fluence\n'+
'! Output : Eq  : Wigner strain\n'+
'      implicit none\n'+
'      real*8 T,F,Ew(6),dL_par,dL_per\n'+
'\n'+
'      call Get_Wigner_par(T,F,dL_par)\n'+
'      call Get_Wigner_per(T,F,dL_per)\n'+
#"      write(*,*) 'DL PAR:',dL_par,'| DL PER:',dL_per\n"+        
'      Ew = 0.d0\n'+
'      Ew(1) = dL_{}\n'.format(Term1.lower())+
'      Ew(2) = dL_{}\n'.format(Term2.lower())+
'      Ew(3) = dL_{}\n'.format(Term3.lower())+
'\n'+
'      return\n'+
'      end\n')
   
    return CODE

def Get_EthCODE(Term1,Term2,Term3,Ti):
    CODE = (
'      subroutine Get_Eth(T,F,Eth)\n'+
'!*************************************************************************\n'+
'! This subroutine computes the thermal strain at the end of the time step\n'+
'! Input  : T   : Temperature\n'+
'!          F   : Fluence\n'+
'! Output : Eth : Thermal strain\n'+
'      implicit none\n'+
'      real*8 T,F,Eth(6),CTE0(3),Ti,b1,b2,b3,b4,b5,b6\n'+
'      real*8 a0,a1,a2,CTE(3),Scale,CTE_par,CTE_per,CTE0_par,CTE0_per\n'+
'      call Get_CTEoCTE0par(T,F,CTE_par)\n'+
'      call Get_CTEoCTE0per(T,F,CTE_per)\n'+
'      call Get_CTE0par(T,CTE0_par)\n'+
'      call Get_CTE0per(T,CTE0_per)\n'+        
'\n'+
#"      write(*,*) 'CTE PAR:',CTE_par,'| CTE PER:',CTE_per\n"+
'      Eth = 0.d0\n'+
'      Eth(1) = CTE0_{}*CTE_{}*(T-{})\n'.format(Term1,Term1,fcode(Ti).strip())+
'      Eth(2) = CTE0_{}*CTE_{}*(T-{})\n'.format(Term2,Term2,fcode(Ti).strip())+
'      Eth(3) = CTE0_{}*CTE_{}*(T-{})\n'.format(Term3,Term3,fcode(Ti).strip())+
'      return\n'+
'      end\n')
   
    return CODE

def GetPrimaryRateCode(P,PrimaryCreepRate_Form,PrimaryCreepRate_Form_dEcp,PrimaryCreepRate_Form_dS):
    from sympy import symbols
    iDc_S, Ecp, iDc, T, F = symbols('iDc_S Ecp iDc T F')
    
    CODE = (
'!***********************************************************************\n'+
'      subroutine GetEcp_rate(T,F,stress,Ecp,Ecp_rate,dEcp_rate_dS,\n'+
'     .                       dEcp_rate_dE)\n'+
'!***********************************************************************\n'+
'! This subroutine computes the primary creep rate given temperature,\n'+
'! fluence and stress\n'+
'! Input  : T            : Temperature\n'+
'!          F            : Fluence\n'+
'!          stress       : Stress\n'+
'!          Ecp          : Primary creep strain\n'+
'! Output : Ecp_rate     : Primary creep strain rate (wrt fluence)\n'+
'!          dEcp_rate_dS : Derivative of primary creep strain rate wrt stress\n'+
'!          dEcp_rate_dE : Scalar derivative of primary creep strain rate\n'+
'!                         wrt primary creep strain\n'+
'      implicit none\n'+
'      integer i,j\n'+
'      real*8 T,F,stress(6),Ecp(6),Ecp_rate(6),alpha,G0,dEcp_rate_dE\n'+
'      real*8 iDel(6,6),iDc(6,6),iDc_S(6),EoverE0_ZF,dEcp_rate_dS(6,6)\n'+
'\n'
'      call Get_invDel(T,0.d0,iDel)\n'+
'      iDc = iDel\n'+
'      stress = stress\n'+        
'      iDc_S = 0.d0\n'+
'      do i=1,6\n'+
'         do j=1,6\n'+
'            iDc_S(i) = iDc_S(i) + iDc(i,j)*stress(j)\n'+
'         enddo\n'+
'      enddo\n'+
'      Ecp_rate = {}\n'.format(fcode(eval(PrimaryCreepRate_Form)).strip())+
'      dEcp_rate_dE = {}\n'.format(fcode(eval(PrimaryCreepRate_Form_dEcp)).strip())+
'      dEcp_rate_dS = {}\n'.format(fcode(eval(PrimaryCreepRate_Form_dS)).strip())+
'      return\n'+
'      end\n'+
'!************************************************************************\n')
    return CODE

def GetSecondaryRateCode(P,SecondaryCreepRate_Form,SecondaryCreepRate_Form_dS):
    from sympy import symbols
    iDc_S, Ecp, iDc, T, F = symbols('iDc_S Ecp iDc T F')


    CODE = (
'!************************************************************************\n'+
'      subroutine GetEcs_rate(T,F,stress,Ecs_rate,dEcsrate_dS)\n'+
'!************************************************************************\n'+
'! This subroutine computes the secondary creep rate, given the current\n'+
'! temperature, fluence and stress\n'+
'! Input :  T           : Temperature\n'+
'!          F           : Fluence\n'+
'!          stress      : Stress\n'+
'! Output : Ecs_rate    : Secondary creep rate\n'+
'!          dEcsrate_dS : Derivative of the seondary creep rate wrt stress\n'+
'      implicit none\n'+
'\n'
'      integer i,j\n'+
'      real*8 T,F,stress(6),Ecs_rate(6)\n'+
'      real*8 K,SC_a,SC_e,SC_k,beta,iDel(6,6),iDc(6,6),EoverE0_ZF\n'+
'      real*8 iDc_S(6),dEcsrate_dS(6,6)\n'+
'\n'+      
'      call Get_invDel(T,0.d0,iDel)\n'+
'      iDc = iDel\n'+
'      stress = stress\n'+     
'\n'+
'      iDc_S = 0.d0\n'+
'      do i=1,6\n'+
'         do j=1,6\n'+
'            iDc_S(i) = iDc_S(i) + iDc(i,j)*stress(j)\n'+
'         enddo\n'+
'      enddo\n'+
'\n'+
'      Ecs_rate    = {}\n'.format(fcode(eval(SecondaryCreepRate_Form)).strip())+
'      dEcsrate_dS = {}\n'.format(fcode(eval(SecondaryCreepRate_Form_dS)).strip())+
'      return\n'+
'      end\n')
    return CODE

def Get_Code(SubName,Expression,Tmin,Tmax,Fmin,Fmax):
    CODE = (
'      subroutine {}(To,Fo,Answer)\n'.format(SubName)+
'!*************************************************************************\n'+
'! This subroutine computes the wigner strain in the perpendicular direction\n'+
'! Input  : To      : Temperature\n'+
'!          Fo      : Fluence\n'+
'! Output : dL_per : Wigner strain in perpendicular direction\n'+
'      implicit none\n'+
'      integer i,j\n'+
'      real*8 To,Fo,Answer\n'+
'      real*8 T,F,Tmin,Tmax,Fmin,Fmax\n'+
'\n'+
'      Tmin = {}\n'.format(fcode(Tmin).strip())+
'      Tmax = {}\n'.format(fcode(Tmax).strip())+
'      Fmin = {}\n'.format(fcode(Fmin).strip())+
'      Fmax = {}\n'.format(fcode(Fmax).strip())+
'\n'+
'      T = (To-Tmin)/(Tmax-Tmin)\n'+
'      F = (Fo-Fmin)/(Fmax-Fmin)\n'+ 
'\n'+
'{}\n'.format(Expression)+
'      return\n'
'      end\n')
    return CODE

def Get_Code1D(SubName,Expression,Tmin,Tmax):
    CODE = (
'      subroutine {}(To,Answer)\n'.format(SubName)+
'!*************************************************************************\n'+
'! This subroutine computes the wigner strain in the perpendicular direction\n'+
'! Input  : To      : Temperature\n'+
'!          Fo      : Fluence\n'+
'! Output : dL_per : Wigner strain in perpendicular direction\n'+
'      implicit none\n'+
'      integer i,j\n'+
'      real*8 To,Answer\n'+
'      real*8 T,F,Tmin,Tmax\n'+
'\n'+
'      Tmin = {}\n'.format(fcode(Tmin).strip())+
'      Tmax = {}\n'.format(fcode(Tmax).strip())+
'\n'+
'      T = (To-Tmin)/(Tmax-Tmin)\n'+
'\n'+
'{}\n'.format(Expression)+
'      return\n'
'      end\n')
    return CODE
