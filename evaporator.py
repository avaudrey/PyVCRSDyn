#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Python package dedicated to calculations dealing with the dynamic behaviour of
evaporator involved in vapor-compression refrigeration systems.

Author: 
"""

import numpy as np

# If CoolProp is installed, it will be used for the calculation of physical
# properties of water
try:
    import CoolProp.CoolProp as cp
    is_coolprop_present = True
except ImportError:
    is_coolprop_present = False

class Duct:
    """
    A duct is the geometrical part of a fluid flow. Different ducts can be put
    together to form a heat exchanger, each of them possibly housing a different
    type of fluid flow, a single phase or a two phase one.
    """
    def __init__(self, length , cross_sectional_area, wet_perimeter):
        # Length of the concerned duct, in meter
        self.length = length
        # Cross sectional area of the duct, in square meter
        self.cross_sectional_area = cross_sectional_area
        # Wet perimeter, in meter
        self.wet_perimeter = wet_perimeter
    def volume(self):
        """Inner volume of the tube, in cubic meter."""
        return self.length*self.cross_sectional_area
    def heat_exchange_area(self):
        """Heat exchange area, corresponding to its inner surface, in square
        meter."""
        return self.wet_perimeter*self.length

class HeatExchanger(Duct):
    """
    Heat exchanger composed by one or several duct housing fluid flows.
    """
    def __init__(self, wall_thickness, average_specific_heat_capacity, mass,\
                 wall_thermal_conductivity):
        # Thickness of the wall separating two different ducts, in meter.
        self.wall_thickness = wall_thickness
        # Average heat capacity of the heat exchanger, in J/(kg.K).
        self.average_specific_heat_capacity = average_specific_heat_capacity
        # Total mass of the heat exchanger, in kilogram.
        self.mass = mass
        # Thermal conductivity of the wall separating two different ducts.
        self.wall_thermal_conductivity = wall_thermal_conductivity
    def total_heat_capacity(self):
        """ Total heat capacity of the heat exchanger, in J/K."""
        return self.mass*self.average_heat_capacity
    def conductive_thermal_resistance(self):
        # TODO : to finish.
        pass

class Fluid:
    """
    Fluid involved in single or two-phase flows.
    """
    def __init__(self, fluid_name, pressure, temperature=None):
        # XXX : Why did you affect the value "None" to the temperature?
        # Name of the concerned fluid, as a string
        self.name = fluid_name
        # Values of the fluid pressure and temperature, in Pa and K,
        # respectively
        self.pressure = pressure
        self.temperature = temperature
    def saturated_vapor_density(self):
        # Density of the fluid as saturated vapor, in kg/m3
        return cp.PropsSI('D',\
                          'P', self.pressure,\
                          'Q', 1.0,\
                          self.name)
    def saturated_liquid_density(self):
        # Density of the fluid as saturated liquid, in kg/m3
        return cp.PropsSI('D',\
                          'P', self.pressure,\
                          'Q', 0.0,\
                          self.name)
    def density(self):
        # Actual density of the fluid, as a function of its pressure and
        # temperature, in kg/m3
        return cp.PropsSI('D',\
                          'P', self.pressure,\
                          'T', self.temperature,\
                          self.name)
    def saturated_liquid_specific_enthalpy(self):
        # Mass specific enthalpy of the fluid, in J/kg, when in a saturated
        # liquid phase
        return cp.PropsSI('H',\
                          'P', self.pressure,\
                          'Q', 0.0 ,\
                          self.name)
    def saturated_vapor_specific_enthalpy(self):
        # Mass specific enthalpy of the fluid, in J/kg, when in a saturated
        # vapor phase
        return cp.PropsSI('H',\
                          'P', self.pressure,\
                          'Q', 1.0,\
                          self.name)
    def specific_enthalpy(self):
        # Actual mass specific enthalpy of the fluid, as a function of its
        # pressure and temperature, in J/kg
        return cp.PropsSI('H',\
                          'P', self.pressure,\
                          'T', self.temperature,\
                          self.name)
    def Prandtl_number(self):
        # Prandtl number of the fluid
        return cp.PropsSI('PRANDTL',\
                          'P', self.pressure,\
                          'T', self.temperature,\
                          self.name)
    def conductivity(self):
        # TODO : to finish
        #return cp.PropsSI()
        pass
    ########################## Partial derivatives ############################
    def dRhol_dP(self):
        return cp.PropsSI('d(Dmass)/d(P)|sigma','P',self.pressure,'Q',0,self.which_fluid)
    def dRhov_dP(self):
        return cp.PropsSI('d(Dmass)/d(P)|sigma','P',self.pressure,'Q',1,self.which_fluid)
    def dHl_dP(self):
        return cp.PropsSI('d(Hmass)/d(P)|sigma','P',self.pressure,'Q',0,self.which_fluid)
    def dHv_dP(self):
        return cp.PropsSI('d(Hmass)/d(P)|sigma','P',self.pressure,'Q',1,self.which_fluid)

R134a=fluid('R134a',250_000)
h_v=R134a.dRhol_dP()
print(h_v)

# class two_phase:
#   def __init__(self,length):
#     self.length=length
  
#   def Nusselt_number(self):
#     Nu=
#     return Nu

#   def Reynolds(self):
#     Re=
#     return Re 

# ##### Realizando la programacion de lo encontrado en papers #########

# class evap:
#   Z11=
#   z12=
#   z21=
#   z22=
#   z23=
#   z31=
#   z32=
#   z33=
#   z44=
#   z51=
#   z55=

# class steady_state:
#   def __init__(self.heat_transfer_Area_in,self.heat_transfer_Area_out,self.Temperature_cold_fluid_in,self.Temperature_cold_fluid_out,self.Temperature_hot_fluid_in,self.Temperature_hot_fluid_out):
#     self.heat_transfer_Area_in = heat_transfer_Area_in
#     self.heat_transfer_Area_out = heat_transfer_Area_out
#     self.Temperature_cold_fluid_in = Temperature_cold_fluid_in
#     self.Temperature_cold_fluid_out = Temperature_cold_fluid_out
#     self.Temperature_hot_fluid_in = Temperature_hot_fluid_in
#     self.Temperature_hot_fluid_out = Temperature_hot_fluid_out
#     self.Diameter_in =
     
# def L_length(self):
#   return

"""  Class transient is commented 
class transient:
  
  def __init__(self,alpha_in, alpha_out, Total_Lenght_L, Lenght_l):
    self.alpha_in = alpha_in
    self.alpha_out = alpha_out
    self.Total_Lenght_L = Total_Lenght_L
    self.Lenght_l = Lenght_l
  
  def Matrix(self):
    z11=
    z12=
    z21=
    z22=
    z23=
    z31=
    z32=
    z33=
    z44=
    z51=
    z55=
"""