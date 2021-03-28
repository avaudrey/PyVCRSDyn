# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import CoolProp.CoolProp as cp
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


class duct:
    """
    A duct is the geometrical part of a fluid flow. Different ducts can be put
    together to form a heat exchanger, each of them possibly housing a different
    type of fluid flow, a single phase or a two phase one.
    """
    def __init__(self,length,cross_sectional_area,wet_perimeter):
        # Length of the concerned duct, in meter
        self.length = length
        # Cross sectional area of the duct, in square meter
        self.cross_sectional_area=cross_sectional_area
        # Wet perimeter, in meter
        self.wet_perimeter=wet_perimeter
  
    def volume(self):
        """Inner volume of the tube, in cubic meter."""
        return self.length*self.cross_sectional_area
  
    def heat_exchange_area(self):
        """Heat exchange area, corresponding to its inner surface, in square
        meter."""   
        return self.wet_perimeter*self.length

class exchanger:
    def __init__(self,wall_thickness,average_heat_capacity,mass,wall_thermal_conductivity):
        self.wall_thickness=wall_thickness
        self.average_heat_capacity=average_heat_capacity
        self.mass=mass
        self.wall_thermal_conductivity=wall_thermal_conductivity

    def total_heat_capacity(self):
        return self.mass*self.average_heat_capacity
  
  #def conductive_thermal_resistance(self):


class fluid:
    """
    Fluid involved in single or two-phase flows.
    """
    def __init__(self,fluid_name,pressure,temperature=None):
        # XXX : Why did you affect the value "None" to the temperature?
        # Name of the concerned fluid, as a string
        self.fluid_name=fluid_name
         # Values of the fluid pressure and temperature, in Pa and K,
        # respectively
        self.pressure=pressure
        self.temperature=temperature
  
    def Rho_v(self):
        # Density of the fluid as saturated vapor, in kg/m3
        return cp.PropsSI('D',\
                          'P',self.pressure,\
                          'Q',1,\
                          self.fluid_name)
  
    def Rho_l(self):
        # Density of the fluid as saturated liquid, in kg/m3
        return cp.PropsSI('D',\
                          'P',self.pressure,\
                          'Q',0,\
                          self.fluid_name)

    def Rho(self):
        # Actual density of the fluid, as a function of its pressure and
        # temperature, in kg/m3
        return cp.PropsSI('D',\
                          'P',self.pressure,\
                          'T',self.temperature,\
                          self.fluid_name)

    def h_l(self):
        # Mass specific enthalpy of the fluid, in J/kg, when in a saturated
        # liquid phase
        return cp.PropsSI('H',\
                          'P',self.pressure,\
                          'Q',0,\
                          self.fluid_name)

    def h_v(self):
        # Mass specific enthalpy of the fluid, in J/kg, when in a saturated
        # vapor phase
        return cp.PropsSI('H',\
                          'P',self.pressure,\
                          'Q',1,\
                          self.fluid_name)

    def h(self):
        # Actual mass specific enthalpy of the fluid, as a function of its
        # pressure and temperature, in J/kg
        return cp.PropsSI('H',\
                          'P',self.pressure,\
                          'T',self.temperature,\
                          self.fluid_name)
  
    def Prandtl(self):
        # Prandtl number of the fluid
        return cp.PropsSI('PRANDTL',\
                          'P',self.pressure,\
                          'T',self.temperature,\
                          self.fluid_name)

    def conductivity(self):
        # Conductivity of the fluid
        return cp.PropsSI('CONDUCTIVITY',\
                          'P',self.pressure,\
                          'T',self.temperature,\
                          self.fluid_name)

  ########################## Partial derivatives ############################

    def dRhol_dP(self):
        # Partial derivative of saturated liquid density
        # respect to pressure
        return cp.PropsSI('d(Dmass)/d(P)|sigma',\
                          'P',self.pressure,\
                          'Q',0,\
                          self.fluid_name)
  
    def dRhov_dP(self):
        # Partial derivative of saturated vapor density
        # respect to pressure
        return cp.PropsSI('d(Dmass)/d(P)|sigma',\
                          'P',self.pressure,\
                          'Q',1,\
                          self.fluid_name)

    def dRho_dP(self,enthalpy=None):
        # Partial derivative of density
        # respect to pressure, when is function of
        # Pressure and enthalpy
        return cp.PropsSI('d(Dmass)/d(P)|Hmass',\
                          'P',self.pressure,\
                          'H',enthalpy,\
                          self.fluid_name)
  
    def dRho_dH(self,enthalpy=None):
        # Partial derivative of density
        # respect to enthalpy, when is function of
        # Pressure and enthalpy
        return cp.PropsSI('d(Dmass)/d(Hmass)|P',\
                          'P',self.pressure,\
                          'H',enthalpy,\
                          self.fluid_name)

    def dHl_dP(self):
        # Partial derivative of saturated liquid enthalpy
        # respect to pressure
        return cp.PropsSI('d(Hmass)/d(P)|sigma',\
                          'P',self.pressure,\
                          'Q',0,\
                          self.fluid_name)
  
    def dHv_dP(self):
        # Partial derivative of saturated vapor enthalpy
        # respect to pressure
        return cp.PropsSI('d(Hmass)/d(P)|sigma',\
                          'P',self.pressure,\
                          'Q',1,\
                          self.fluid_name)

  
class Heat_transfer:
    def __init__(self, alpha_in, alpha_out, Heat_Transfer_length,  Tc_sortie=None, Tc_entree=None, Tf_sortie=None, Tf_entree=None):
        self.alpha_in = alpha_in
        self.alpha_out = alpha_out
        self.Heat_Transfer_length = Heat_Transfer_length
        self.Tc_sortie = Tc_sortie
        self.Tc_entree = Tc_entree
        self.Tf_sortie = Tf_sortie
        self.Tf_entree = Tf_entree
  
    def UA(self): # Ver si requiero poner lo de la conductividad y esas cosas locas
        diametro = (4* ducto_total.cross_sectional_area / (np.pi))**0.5 
        Area =  np.pi * diametro * ducto_total.length
        return 1 / ( 1/(self.alpha_in * Area * self.Heat_Transfer_length / ducto_total.length) + 1/(self.alpha_out * Area * self.Heat_Transfer_length / ducto_total.length))
  
    def DMLT(self):
        dt1 = self.Tc_entree - self.Tf_sortie
        dt2 = self.Tc_sortie - self.Tf_entree
        return (dt1 - dt2) / (np.log(dt1/dt2))
  
    def Heat(self):
        return self.DMLT()*self.UA()

fluido1 = 'R134a'
fluido2 = 'WATER'
# Datos requeridos por el programa: 
d_in = 0.015
d_ext = 0.025
ducto_total = duct(5.277137578892268,np.pi * d_in**2/4,1) # el wet_perimeter no lo he usado 
length_l_petit= 4.3806496889366455
m_water = 1.0
cp_water = 4178.740883004221
h_4 = 256409.24455736773 # este valor es una constante??? Pues no!!!!
Aex = np.pi * d_in * ducto_total.length
Tcaliente_in = 40 + 273.15
m_ref = 0.27034739892528964
# mdot_4 = m_ref
# mdot_1 = m_ref
m_ine = m_water
m_oe = m_water
p_ext = 300_000
h_hi = cp.PropsSI('H', 'P', p_ext,'T',Tcaliente_in,fluido2)

alpha_SH_in = 1744.3426527250208
alpha_TP_in = 11011.778404079098
alpha_SH_out = 8109.623096189451
alpha_TP_out = alpha_SH_out

tau =0.3 

class transient:
  
    def __init__(self,alpha_SH, alpha_TP, Total_Lenght_L, Lenght_l,Void_fraction):
        self.alpha_SH = alpha_SH
        self.alpha_TP = alpha_TP
        self.Total_Lenght_L = Total_Lenght_L
        self.Lenght_l = Lenght_l
        self.Void_fraction = Void_fraction
  
    def model(self,Z,t,mdot_1): 
        H_g = ( fluid(fluido1,Z[0]).h_v() + Z[1] ) / 2 
        T_g = cp.PropsSI('T','P',Z[0],'H',H_g,fluido1)
        T_4 = cp.PropsSI('T','P',Z[0],'Q',0,fluido1)
        Rho_g = ( fluid(fluido1,Z[0],T_g).Rho() )
        
        Tf_out = cp.PropsSI('T','P', Z[0], 'H', Z[1], fluido1)
        Tc_out = cp.PropsSI('T','P',Z[3],'H',Z[4],fluido2)
        Tc_l = T_4 + ( Tc_out - T_4 )/np.exp(-Heat_transfer(alpha_TP_in,alpha_TP_out,Z[2]).UA() * Z[2]/(ducto_total.length * m_water * cp_water) )

        T_ex =  T_4*Z[2]/ducto_total.length + (Tcaliente_in - T_4)*(m_water * cp_water)/(Heat_transfer(alpha_TP_in,alpha_TP_out,Z[2]).UA())*(1-np.exp( (-Heat_transfer(alpha_TP_in,alpha_TP_out,Z[2]).UA()*Z[2]) / (ducto_total.length*m_water*cp_water) )) \
            + T_g * (ducto_total.length - Z[2])/ducto_total.length - (Tc_l - T_g) * (m_water * cp_water)/(Heat_transfer(alpha_SH_in,alpha_SH_out,Z[2]).UA())*(np.exp( (-Heat_transfer(alpha_SH_in,alpha_SH_out,Z[2]).UA()*Z[2]) / (ducto_total.length*m_water*cp_water) ) - np.exp( (-Heat_transfer(alpha_SH_in,alpha_SH_out,Z[2]).UA() / (m_water*cp_water) ))) # debo definir la manera de determinar la temperatura promedio del fluido externo
        h_ex = cp.PropsSI('H','P', Z[3], 'T', T_ex, fluido2)  

        Q_l = Heat_transfer(alpha_TP_in,alpha_TP_out,Z[2],Tc_out,Tc_l,T_4,T_4).Heat() 
        Q_Ll = Heat_transfer(alpha_SH_in,alpha_SH_out,ducto_total.length-Z[2],Tc_l,Tcaliente_in,Tf_out,T_4).Heat()
        Q_total = Q_l + Q_Ll

        z11 = ( (self.Void_fraction * fluid(fluido1,Z[0]).dRhov_dP()+(1-self.Void_fraction)*fluid(fluido1,Z[0]).dRhol_dP())) * ducto_total.cross_sectional_area * Z[2] \
           + ducto_total.cross_sectional_area * (ducto_total.length-Z[2]) * (fluid(fluido1,Z[0]).dRho_dP(H_g) + 0.5 * fluid(fluido1,Z[0]).dRho_dH(H_g) *fluid(fluido1,Z[0]).dHv_dP())    
        z12 = ducto_total.cross_sectional_area*(ducto_total.length-Z[2])*0.5*fluid(fluido1,Z[0]).dRho_dH(H_g)
        z13 = (self.Void_fraction*fluid(fluido1,Z[0]).Rho_v()+(1-self.Void_fraction) * fluid(fluido1,Z[0]).Rho_l()-Rho_g)*ducto_total.cross_sectional_area  
        z24 = Aex * ducto_total.length * fluid(fluido2,Z[3],T_ex).dRho_dP(h_ex) # Ver lo de si entalpia o temperatura
        z25 = Aex * ducto_total.length * 0.5 * fluid(fluido2,Z[3],T_ex).dRho_dH(h_ex) # Ver lo de si entalpia o temperatura
        z34 = Aex * ducto_total.length * (h_ex * fluid(fluido2,Z[3],T_ex).dRho_dP(h_ex) - 1)
        z35 = 0.5 * Aex * ducto_total.length * (h_ex * fluid(fluido2,Z[3],T_ex).dRho_dH(h_ex) + fluid(fluido2,Z[3],T_ex).Rho())
        z41 = ( self.Void_fraction * fluid(fluido1,Z[0]).dHv_dP() * fluid(fluido1,Z[0]).Rho_v() + (1-self.Void_fraction) * fluid(fluido1,Z[0]).dRhol_dP()*(fluid(fluido1,Z[0]).h_l()-fluid(fluido1,Z[0]).h_v())+(1-self.Void_fraction)*fluid(fluido1,Z[0]).dHl_dP()*fluid(fluido1, Z[0]).Rho_l()-1)*ducto_total.cross_sectional_area*Z[2]   
        z43 = ( (1-self.Void_fraction) * fluid(fluido1,Z[0]).Rho_l() * (fluid(fluido1,Z[0]).h_l()-fluid(fluido1,Z[0]).h_v())-Z[0])*ducto_total.cross_sectional_area 
        z51 = ( fluid(fluido1,Z[0]).dRho_dP(H_g) * (H_g - fluid(fluido1,Z[0]).h_v()) + fluid(fluido1,Z[0]).dRho_dH(H_g) * 0.5 * fluid(fluido1,Z[0]).dHv_dP() * (H_g-fluid(fluido1,Z[0]).h_v()) \
           + 0.5 * Rho_g * fluid(fluido1,Z[0]).dHv_dP() - 1 ) * ducto_total.cross_sectional_area * (ducto_total.length - Z[2]) 
        z52 = ( 0.5 * fluid(fluido1,Z[0]).dRho_dP(H_g) * (H_g - fluid(fluido1,Z[0]).h_v()) + 0.5 * Rho_g) * ducto_total.cross_sectional_area * (ducto_total.length - Z[2])
        z53 = ( Z[0] + Rho_g * ( fluid(fluido1,Z[0]).h_v() - H_g ) ) * ducto_total.cross_sectional_area
        f1 = Z[5]-mdot_1
        f2 = m_ine-m_oe
        f3 = -Q_total + m_ine * h_hi - m_oe * Z[4] #Revisar la convencion de signos de todas las ecuaciones :(
        f4 = Q_l + Z[5] * (h_4 - fluid(fluido1,Z[0]).h_v())
        f5 = Q_Ll + mdot_1 * ( fluid(fluido1,Z[0]).h_v() - Z[1] )

        dPdt = f1*z43*z52/(z11*z43*z52 + z12*z41*z53 - z12*z43*z51 - z13*z41*z52) \
            + f4*(z12*z53 - z13*z52)/(z11*z43*z52 + z12*z41*z53 - z12*z43*z51 - z13*z41*z52) \
            - f5*z12*z43/(z11*z43*z52 + z12*z41*z53 - z12*z43*z51 - z13*z41*z52)

        dH_outdt = f1*(z41*z53-z43*z51)/(z11*z43*z52 + z12*z41*z53 - z12*z43*z51 - z13*z41*z52) \
            + f4*(-z11*z53 + z13*z51)/(z11*z43*z52 + z12*z41*z53 - z12*z43*z51 - z13*z41*z52) \
                + f5*(z11*z43-z13*z41)/(z11*z43*z52 + z12*z41*z53 - z12*z43*z51 - z13*z41*z52)

        dldt = -f1*(z41*z52)/(z11*z43*z52 + z12*z41*z53 - z12*z43*z51 - z13*z41*z52) \
            + f4*(z11*z52 - z12*z51)/(z11*z43*z52 + z12*z41*z53 - z12*z43*z51 - z13*z41*z52) \
                + f5*(z12*z41)/(z11*z43*z52 + z12*z41*z53 - z12*z43*z51 - z13*z41*z52)

        dPexdt = (f2*z35)/(z24*z35-z25*z34) - (f3*z25)/(z24*z35-z25*z34)

        dh_hodt= (-f2*z34)/(z24*z35-z25*z34) + (f3*z24)/(z24*z35-z25*z34)

        dmdot_4dt = (mdot_1 - Z[5])/tau 

        return [dPdt,dH_outdt,dldt, dPexdt, dh_hodt, dmdot_4dt]
  
  # initial condition
p0 = 243342.36987785524
h_out0 = 404367.02095548273
l_0 = 4.3806496889366455
Pex0= 300_000
h_ho0= 127792.33926518963
mdot_40 = m_ref

z0=[p0, h_out0, l_0, Pex0, h_ho0, mdot_40]

# Number of points
n =141

# time points
t = np.linspace(0,(n-1)/10,n)

# step input
mdot_1 = np.zeros(n)
mdot_1[0:]=m_ref
# change to 0.3 at time = 5.0
mdot_1[45:] = 0.27

test=transient(1744.307640564927, 7757.236902536498,0.8808487889318191 + 4.421219527097692,1.0973653460950592, 0.41)

# store solution
P = np.empty_like(t)
H_out = np.empty_like(t)
l = np.empty_like(t)
Pex= np.empty_like(t)
h_ho = np.empty_like(t)
mdot_4 = np.empty_like(t)

# record initial conditions
P[0] = z0[0]
H_out[0] = z0[1]
l[0] =z0[2]
Pex[0]= z0[3]
h_ho[0] =z0[4]
mdot_4[0] =z0[5]

# solve ODE
for i in range(1,n):
    # span for next time step
    tspan = [t[i-1],t[i]]
    # solve for next step
    Z = odeint(test.model,z0,tspan,args=(mdot_1[i],))
    # store solution for plotting
    P[i] = Z[1][0]
    H_out[i] = Z[1][1]
    l[i] = Z[1][2]
    Pex[i] = Z[1][3]
    h_ho[i] = Z[1][4]
    mdot_4[i] = Z[1][5]
    # next initial condition
    z0 = Z[1]

# plot results
"""plt.plot(t,mdot_1,'g:',label='u(t)')"""
plt.plot(t,P,'b-',label='x(t)')
"""plt.plot(t,H_out,'r--',label='y(t)')
plt.plot(t,l)
plt.ylabel('values')
plt.xlabel('time')
plt.legend(loc='best')"""
plt.show()
