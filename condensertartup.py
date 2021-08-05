import numpy as np
import scipy.optimize
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import plotly.graph_objects as go


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