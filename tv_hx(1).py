import numpy as np
import scipy.optimize
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import plotly.graph_objects as go

def ThrottleV(Tin,Po,P):
  Tref=75.0 #Reference Temperature (K)

  Z=1.0 #moles
  ZO2=0.21*Z #O2 moles (21% O2)
  ZN2=Z-ZO2 #N2 moles (79% N2)

  #Energy Balance in throttle valve
  #THout=120 #(K) initial T hot outlet value (HX), T inlet TV
  #t1=[]
  #while loop
  #while (abs(error2))>0.1:
    #Tin=THout #Temperature inlet of Throttle valve
  
  #TO DELETE ALL COMMENTS ABOVE (DOESNT MAKE SENSE). Or replace it with the right stuff
    
  def objective(X):
    L=X[0] #liquid moles
    xO2=X[1] #oxygen liquid molar fraction
    Tout=X[2] #Outlet T
    yO2=(ZO2-L*xO2)/(Z-L) #oxygen vapor molar fraction
    yN2=1.0-yO2 # nitrgogen vapor molar fraction
    TrO2=Tout/154.58 #Reduced Temperature O2
    TrN2=Tout/126.2 #Reduced temperature N2
    V=Z-L #Total vapour moles
    Psat=yO2*PsatO2(Tin)+yN2*PsatN2(Tin) #Saturation pressure of the mixture
    
    if Psat<=29.816: #?????
      E_B=V1*y1O2*CpO2V(Tin)*(Tin-Tref)+V1*y1N2*CpN2V(Tin)*(Tin-Tref)+L1*x1O2*CpO2Liq(Tin)*(Tin-Tref)+L1*x1N2*CpN2Liq(Tin)*(Tin-Tref)-L*xO2*(HvapO2(TrO2)+CpO2Liq(Tout)*(Tout-Tref))-L*(1-xO2)*(HvapN2(TrN2)+CpN2Liq(Tout)*(Tout-Tref))+60*28.9*V1-15*28.9*L#-(V)*yN2*(CpN2V(Tout)*(Tout-Tref))-(V)*yO2*(CpO2V(Tout)*(Tout-Tref))
    else:
      E_B=ZO2*CpO2V(Tin)*(Tin-Tref)+ZN2*CpN2V(Tin)*(Tin-Tref)-(V)*yN2*(CpN2V(Tout)*(Tout-Tref))-(V)*yO2*(CpO2V(Tout)*(Tout-Tref))
    return E_B
    
  def constraint1(X): #oxygen raoults law
    L=X[0] #liquid moles total
    xO2=X[1] #O2 moles
    Tout=X[2] #outlet valve temperature
    yO2=(ZO2-L*xO2)/(Z-L) #vapour moles O2
    yN2=1-yO2 #vapour moles N2
    return xO2*PsatO2(Tout)-yO2*P #Raoult's law

  def constraint2(X): #nitrogen raoults law
    L=X[0] #liquid moles total 
    xO2=X[1] #liquid moles oxygen
    Tout=X[2] #valve outlet temperature 
    yO2=(ZO2-L*xO2)/(Z-L) #oxygen vapour moles
    yN2=1-yO2 #nitrogen vapour moles
    xN2=1-xO2 #nitrogen liquid moles
    return xN2*PsatN2(Tout)-yN2*P #nitrogen raoults law

  def constraint3(X): #valve energy balance
    L=X[0] #liquid moles
    xO2=X[1] #oxygen liquid molar fraction
    Tout=X[2] #Outlet T
    yO2=(ZO2-L*xO2)/(Z-L) #oxygen vapor molar fraction
    yN2=1.0-yO2 # nitrgogen vapor molar fraction
    TrO2=Tout/154.58 #oxygen reduced temperature 
    TrN2=Tout/126.2 #nitrogen reduced temeprature
    V=Z-L #total vapour mples 
    Psat=yO2*PsatO2(Tin)+yN2*PsatN2(Tin) #mixture saturation pressure
    
    if Psat<=29.816: 
      E_B=V1*y1O2*CpO2V(Tin)*(Tin-Tref)+V1*y1N2*CpN2V(Tin)*(Tin-Tref)+L1*x1O2*CpO2Liq(Tin)*(Tin-Tref)+L1*x1N2*CpN2Liq(Tin)*(Tin-Tref)-L*xO2*(HvapO2(TrO2)+CpO2Liq(Tout)*(Tout-Tref))-L*(1-xO2)*(HvapN2(TrN2)+CpN2Liq(Tout)*(Tout-Tref))+60*28.9*V1-15*28.9*L##-(V)*yN2*(CpN2V(Tout)*(Tout-Tref))-(V)*yO2*(CpO2V(Tout)*(Tout-Tref))
    else:
      E_B=ZO2*CpO2V(Tin)*(Tin-Tref)+ZN2*CpN2V(Tin)*(Tin-Tref)-V*yN2*(CpN2V(Tout)*(Tout-Tref))-V*yO2*(CpO2V(Tout)*(Tout-Tref))
    return E_B
    
  #Initial values
  x0=[0.5,0.2,126.0] #liquid moles, oxygen liquid molar fraction, outlet temperature of valve
  #Boundaries
  bnds=((0,Z),(0,1),(54,154)) #setting up bounds for the max and min possible values
  #constraints
  con1={'type':'eq','fun':constraint1} #function 1, constraint 1
  con2={'type':'eq','fun':constraint2} #function 2, constraint 2
  con3={'type':'eq','fun':constraint3} #function 3, constraint 3
  con=[con1,con2,con3] #array of constraints
  solmin=minimize(objective,x0,method='SLSQP',bounds=bnds,constraints=con,options={'maxiter':10000}) #solving the objective function with the constraints that f1,f2,f3 must also simultaneously be solved (L,xO2,TOut)
  #print(solmin)
  #Results of the model
  L=solmin.x[0] 
  xO2=solmin.x[1]
  Tout=solmin.x[2]
  xN2=1-xO2
  V=Z-L
  TrO2=Tout/154.58
  TrN2=Tout/126.2
  yO2=(ZO2-L*xO2)/(Z-L) #oxygen vapor molar fraction
  yN2=1.0-yO2 # nitrgogen vapor molar fraction
  Psat=yO2*PsatO2(Tin)+yN2*PsatN2(Tin)
  if Psat<=29.816:
    #Enthalpy of liquid
    Enthalpy_L=L1*x1O2*CpO2Liq(Tin)*(Tin-Tref)+L1*x1N2*CpN2Liq(Tin)*(Tin-Tref)-L*xO2*(-HvapO2(TrO2)+CpO2Liq(Tout)*(Tout-Tref))-L*(1-xO2)*(-HvapN2(TrN2)+CpN2Liq(Tout)*(Tout-Tref))
    Enthalpy_V=V1*y1O2*CpO2V(Tin)*(Tin-Tref)+V1*y1N2*CpN2V(Tin)*(Tin-Tref)-(V)*yN2*(CpN2V(Tout)*(Tout-Tref))-(V)*yO2*(CpO2V(Tout)*(Tout-Tref))#
    
  else:
    
  #Enthalpy of Liquid
    Enthalpy_L=L*xO2*CpO2Liq(Tin)*(Tin-Tref)+L*xN2*CpN2Liq(Tin)*(Tin-Tref)-15*28.9*L-L*xO2*-HvapO2(TrO2)-L*xN2*-HvapN2(TrN2)#J
  #Enthalpy of vapor
    Enthalpy_V=ZO2*CpO2V(Tin)*(Tin-Tref)+ZN2*CpN2V(Tin)*(Tin-Tref)-V*yN2*(CpN2V(Tout)*(Tout-Tref))-V*yO2*(CpO2V(Tout)*(Tout-Tref))
  
  
  return [L,xO2,Tout,xN2,Enthalpy_L,Enthalpy_V]

def hx(TH):
      #Energy balance in heat exchanger
      #Hot fluid
      #Hin-Hout-Q=accumulation
      
      Psat=0.21*PsatO2(TH)+0.79*PsatN2(TH) #Saturation pressure (phase change)
      if Psat<=29.816:
        #Energy balance in HX when there's phase change
        def objective(X):
          L1=X[0] #liquid moles
          x1O2=X[1] #oxygen liquid molar fraction
          TH=X[2] #Outlet T
          y1O2=(ZO2-L1*x1O2)/(Z-L1) #oxygen vapor molar fraction
          y1N2=1.0-y1O2 # nitrogen vapor molar fraction
          TrO2=TH/154.58 #oxygen reduced temperature
          TrN2=TH/126.2 #nitrogen reduced temperature
          V1=Z-L1 #vapour moles
          
          Q=+V*yO2*CpO2V(TCin)*(TCin-Tref)+V*yN2*CpN2V(TCin)*(TCin-Tref)-V*yO2*CpO2V(TCout)*(TCout-Tref)-V*yN2*CpN2V(TCout)*(TCout-Tref)
          E_B=ZO2*CpO2V(THin)*(THin-Tref)+ZN2*CpN2V(THin)*(THin-Tref)-L1*x1O2*(HvapO2(TrO2)+CpO2Liq(TH)*(Tout-Tref))-L1*(1-x1O2)*(HvapN2(TrN2)+CpN2Liq(TH)*(TH-Tref))-(V1)*y1N2*(CpN2V(TH)*(TH-Tref))-(V1)*y1O2*(CpO2V(TH)*(TH-Tref))+60*28.9*Z-15*28.9*L1-V1*60*28.9-Q #
          return E_B
        def constraint1(X):
          L1=X[0]
          x1O2=X[1]
          TH=X[2]
          y1O2=(ZO2-L1*x1O2)/(Z-L1)
          y1N2=1-y1O2
          return x1O2*PsatO2(TH)-y1O2*Po
        def constraint2(X):
          L1=X[0]
          x1O2=X[1]
          TH=X[2]
          y1O2=(ZO2-L1*x1O2)/(Z-L1)
          y1N2=1-y1O2
          x1N2=1-x1O2
          return x1N2*PsatN2(TH)-y1N2*Po
        def constraint3(X):
          L1=X[0] #liquid moles
          x1O2=X[1] #oxygen liquid molar fraction
          TH=X[2] #Outlet T
          y1O2=(ZO2-L1*x1O2)/(Z-L1) #oxygen vapor molar fraction
          y1N2=1.0-y1O2 # nitrgogen vapor molar fraction
          TrO2=TH/154.58
          TrN2=TH/126.2
          V1=Z-L1
          
          
          Q=+V*yO2*CpO2V(TCin)*(TCin-Tref)+V*yN2*CpN2V(TCin)*(TCin-Tref)-V*yO2*CpO2V(TCout)*(TCout-Tref)-V*yN2*CpN2V(TCout)*(TCout-Tref)
          E_B=ZO2*CpO2V(THin)*(THin-Tref)+ZN2*CpN2V(THin)*(THin-Tref)-L1*x1O2*(HvapO2(TrO2)+CpO2Liq(TH)*(Tout-Tref))-L1*(1-x1O2)*(HvapN2(TrN2)+CpN2Liq(TH)*(TH-Tref))-(V1)*y1N2*(CpN2V(TH)*(TH-Tref))-(V1)*y1O2*(CpO2V(TH)*(TH-Tref))+60*28.9*Z-15*28.9*L1-V1*60*28.9-Q #
          return E_B
          
        #Initial values
        x0=[0.5,0.2,TH]
        #Boundaries
        bnds=((0,Z),(0,1),(54,154))
        #constraints
        con1={'type':'eq','fun':constraint1}
        con2={'type':'eq','fun':constraint2}
        con3={'type':'eq','fun':constraint3}
        con=[con1,con2,con3]
        solmin=minimize(objective,x0,method='SLSQP',bounds=bnds,constraints=con,options={'maxiter':1000})
        #print(solmin)
        #Results of the model
        L1=solmin.x[0]
        x1O2=solmin.x[1]
        TH=solmin.x[2]
        V1=Z-L1
        y1O2=(ZO2-L1*x1O2)/(Z-L1) #oxygen vapor molar fraction
        y1N2=1.0-y1O2 # nitrgogen vapor molar fraction
        x1N2=1-x1O2
        
      else:
        #initialize time
        t=0
        #initialize error
        error=100
        #while loop
        while (abs(error))>0.01:
          THout=TH #initial value for the iteration
          #it calculates the new temperature if there's no phase change
          mCpdTH_dt=(ZO2*CpO2V(THin)*(THin-Tref)+ZN2*CpN2V(THin)*(THin-Tref)-ZO2*CpO2V(THout)*(THout-Tref)-ZN2*CpN2V(THout)*(THout-Tref))+V*yO2*CpO2V(TCin)*(TCin-Tref)+V*yN2*CpN2V(TCin)*(TCin-Tref)-V*yO2*CpO2V(TCout)*(TCout-Tref)-V*yN2*CpN2V(TCout)*(TCout-Tref)
          TH=mCpdTH_dt*dt/(ZO2*CpO2V((THout+THin)/2)+ZN2*CpN2V((THout+THin)/2))+THout
          L1=0
          V1=1
          y1O2=0.21
          y1N2=0.79
          x1N2=0
          x1O2=0
          #time iteration 
          t=t+dt
          #calculates the new error
          error=abs(THout-TH)
      return [TH,L1,V1,y1O2,y1N2,x1N2,x1O2]

def CpO2Liq(T):
    #The following properties were calculated using Perry's Handbook (2018) formulas
    CpO2L=(175430.0-6152.3*T+113.92*T**2-0.92382*T**3+0.0027963*T**4)/1000 #J/mol*K oxygen liquid heat capacity
    return CpO2L
def CpN2Liq(T):
  #The following properties were calculated using Perry's Handbook (2018) formulas
  CpN2L=(281970.0-12281.0*T+248.0*T**2-2.2182*T**3+0.0074902*T**4)/1000 #J/mol*K nitrogen liquid heat capacity
  return CpN2L

def CpO2V(Tout):
  #The following properties were calculated using Perry's Handbook (2018) formulas
  CpO2Vap=(0.29103*1.0e5+0.10040*1.0e5*(2.5265*1.0e3/Tout/np.sinh(2.5265*1.0e3/Tout))**2+0.09356*1.0e5*(1153.8/Tout/np.cosh(1153.8/Tout))**2)/1000 #J/mol*K oxygen vapor heat capacity
  return CpO2Vap

def CpN2V(Tout):
  #The following properties were calculated using Perry's Handbook (2018) formulas
  CpN2Vap=(0.29105*1.0e5+0.08615*1.0e5*(1.7016*1.0e3/Tout/np.sinh(1.7016*1.0e3/Tout))**2+0.00103*1.0e5*(909.79/Tout/np.cosh(909.79/Tout))**2)/1000 #J/mol*K nitrogen vapor heat capacity  
  return CpN2Vap
def HvapO2(TrO2):
  #The following properties were calculated using Perry's Handbook (2018) formulas
  HvaporO2=(0.9008*1e7*(1.0-TrO2)**(0.4542-0.4096*TrO2+0.3183*TrO2**2.0))/1000 #J/mol oxygen heat of vaporization
  return HvaporO2
def HvapN2(TrN2):
  #The following properties were calculated using Perry's Handbook (2018) formulas
  HvaporN2=(0.74905*1e7*(1.0-TrN2)**(0.40106-0.317*TrN2+0.27343*TrN2**2.0))/1000 #J/mol nitrogen heat of vaporization
  return HvaporN2
def PsatO2(Tout):
  Psat_O2=np.exp(51.245-1200.2/Tout-6.4361*np.log(Tout)+0.028405*Tout)/1e5 #From Perry's Chemical Engineer Handbook (2018)
  return Psat_O2
def PsatN2(Tout):
  Psat_N2=np.exp(58.282-1084.1/Tout-8.3144*np.log(Tout)+0.044127*Tout)/1e5 #From Perry's Chemical Engineer Handbook (2018)
  return Psat_N2

THin=303 #K
P=2 #bar
Po=30
THout=140 #(K) initial T hot outlet value (HX), T inlet TV
Tref=75.0 #K
Tpinch=0 #
Z=1.0 #moles
ZO2=0.21*Z #moles
ZN2=Z-ZO2 #moles
error2=100
t1=[]
while (abs(error2))>0.1:
  Tin=THout
  [L,xO2,Tout,xN2,Enthalpy_L,Enthalpy_V]=ThrottleV(Tin,Po,P)
  
  TCin=Tout
  TCout=THin-Tpinch
  dt=0.1
  V=Z-L
  yO2=(ZO2-L*xO2)/(Z-L) #oxygen vapor molar fraction
  yN2=1.0-yO2 # nitrgogen vapor molar fraction
  [THout,L1,V1,y1O2,y1N2,x1N2,x1O2]=hx(Tin)
  #t1.append(t)
  error2=abs(THout-Tin)
  
  #ttot=sum(t1)
  #Plot table
for Tout in np.arange(60,300,0.001):

  Psatout=0.21*PsatO2(Tout)+0.79*PsatN2(Tout)
  if round(Psatout,3)==P:
    
    break
  T1=Tout
fig1 = go.Figure(data=[go.Table(header=dict(values=['Outlet Pressure (bar)','Outlet Temperature of TV (K)','Hot Outlet temperature of condenser (K)','Pinch T (K)', 'Liquid moles','liquid molar fraction O2','liquid molar fraction N2','vapor molar fraction O2','vapor molar fraction N2','Liquid enthalpy (J)','Vapor enthalpy (J)']),
                  cells=dict(values=[['%.2f'%P],['%.2f'%T1],['%.2f'%THout],['%.2f'%Tpinch],['%.5f'%L], ['%.3f'%xO2],['%.3f'%xN2],['%.3f'%yO2],['%.3f'%yN2],['%.2f'%Enthalpy_L],['%.2f'%Enthalpy_V]]))
                      ])
fig1.show()
  #Plot P vs T
fig = plt.figure(figsize=[5, 5])
ax = plt.axes() 
Temp=np.linspace(Tin,Tout,100)
Pressure=np.linspace(Po,P,100)
ax.plot(Temp, Pressure, 'b')
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Pressure (bar)')
ax.legend()