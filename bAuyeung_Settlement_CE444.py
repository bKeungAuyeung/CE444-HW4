from math import log10
import numpy as np

class Foundation:
    def __init__(self, B, L, D, H, surcharge = 'Na',thickness = 'Na', point_load = 'Na', gamma = 'Na', metric = True):
        self.B = B
        self.L = L
        self.D = D
        self.H = H
        self.thickness = thickness
        self.gamma = gamma
        if surcharge == 'Na':
            self.load = (point_load / (B*L)) + (gamma*thickness)
        else:
            self.load = surcharge
        if metric == False: # if the system units are in ft
            self.unit_system = 'ft, psf, pcf' 
            self.Pa = 2117 # psf
            self.refL = 3.28 # ft
            self.gamma_w = 62.4 # pcf
        else:
            self.unit_system = 'm, kPa, kN/m^3'
            self.Pa = 100 # kPa
            self.refL = 1 # m
            self.gamma_w = 9.81 # kN/m^3

    def influence_depth(self):
        self.influence_depth = self.refL * ((self.B / self.refL)**0.79)

    def sand(self, gamma_soil, gwt, avgN60, excavation = False):
        self.gamma_soil = gamma_soil # match units specified by the foundation
        self.gwt = gwt
        self.avgN60 = avgN60
        if excavation == False:
            self.excavation = 'No'
        else:
            self.excavation = excavation # excavation must be in specified units

    def clay(self, gamma_soil, gwt, Su, PI, K, e0, OCR = 1, excavation = False):
        self.gamma_soil = gamma_soil # match units specified by the foundation
        self.gwt = gwt
        self.Su = Su
        self.PI = PI
        self.K = K # obtain K from the Duncan and Buchignani (1987) graph
        self.Eu = K * Su
        self.OCR_clay = OCR
        self.e0 = e0
        if excavation == False:
            self.excavation = 'No'
        else:
            self.excavation = excavation # excavation must be in specified units

    def settlement_sand(self, time, deep_soil = False, dynamic = False):
        if deep_soil == True:
            self.layer_thickness_factor = 1
        else:
            H_infl_ratio = self.H / self.influence_depth
            self.layer_thickness_factor = H_infl_ratio * (2 - H_infl_ratio)
        
        self.comp_index = 1.71 / (self.avgN60**1.4)

        print(f'Time of Settlement Considered: {time} yrs')
        if time == 0:
            self.time_factor = 1
        elif dynamic == True:
            R3yr, Rt = 0.7, 0.8
            self.time_factor = 1 + R3yr + (Rt*log10(time/3)) # time must be in yrs
        else:
            R3yr, Rt = 0.3, 0.2
            self.time_factor = 1 + R3yr + (Rt*log10(time/3)) # time must be in yrs

        LB_ratio = self.L / self.B
        self.shape_factor = ((1.25*LB_ratio)/(LB_ratio + 0.25))**2

        if self.gwt == 0:
            self.preconsol_stress = (self.gamma_soil - self.gamma_w)*self.excavation
        elif 0 < self.gwt < self.excavation:
            self.preconsol_stress = (self.gamma_soil * self.excavation) - ((self.gamma_w) * ( self.excavation - self.gwt))
        else:
            self.preconsol_stress = (self.gamma_soil * self.excavation)

        if self.preconsol_stress > self.load:
            #calculate settlement as OCOC
            print('OCOC Settlement Case')
            self.settlement = (self.refL) * (0.1*self.shape_factor*self.layer_thickness_factor*self.time_factor*self.comp_index*((self.B/self.refL)**0.7)) * (self.load / (3*self.Pa))
        elif self.preconsol_stress < self.load:
            #calculate settlement as OCNC
            print('OCNC Settlement Case')
            self.settlement = (self.refL) * (0.1*self.shape_factor*self.layer_thickness_factor*self.time_factor*self.comp_index*((self.B/self.refL)**0.7)) * ((self.load - (2/3*self.preconsol_stress)) / self.Pa)            
        else:
            #calculate settlement as NC
            print('NC Settlement Case')
            self.settlement = (self.refL )* (0.1*self.shape_factor*self.layer_thickness_factor*self.time_factor*self.comp_index*((self.B/self.refL)**0.7)) * (self.load / self.Pa)

    def settlement_clay(self, points: np.array, Cc = 'Na', Cr = 'Na', layer_thickness = 'Na', net = True, I0 = 'Na', I1 = 'Na', skempA = 'Na', skempAlpha = 'Na'):

        if self.gwt == 0:
            self.excav_load = (self.gamma_soil - self.gamma_w)*self.excavation
        elif 0 < self.gwt < self.excavation:
            self.excav_load = (self.gamma_soil * self.excavation) - ((self.gamma_w) * ( self.excavation - self.gwt))
        else:
            self.excav_load = (self.gamma_soil * self.excavation)
        
        if net == True:
            self.load = self.load - self.excav_load
        else:
            self.load = self.load + self.excav_load

        self.veoStress_clay = np.zeros(len(points))

        for index, value in enumerate(points):
            if self.gwt == 0:
                self.veoStress_clay[index] = (self.gamma_soil - self.gamma_w)*value
            elif 0 < self.gwt < value:
                self.veoStress_clay[index] = (self.gamma_soil * value) - ((self.gamma_w) * ( value - self.gwt))
            else:
                self.veoStress_clay[index] = (self.gamma_soil * value)
        
        self.w_i = I0 * I1  * ((self.load * self.B)/self.Eu) #immediate settlement in clay
        
        self.two_to_one_delta = ((self.load * self.B * self.L)/((self.B + (points - self.D))*(self.L + (points - self.D))))

        self.vefStress_clay = self.veoStress_clay + self.two_to_one_delta

        self.preconsol_stress = self.OCR_clay * self.veoStress_clay

        self.w_c = np.zeros(len(self.preconsol_stress))
        for index in range(len(self.preconsol_stress)):
            if self.OCR_clay == 1 or self.veoStress_clay[index] == self.preconsol_stress[index]:   
                print(f'Sublayer {index+1}: NC Settlement Case')
                self.w_c[index] = (layer_thickness/(1 + self.e0)) * Cc * log10(self.vefStress_clay[index]/self.veoStress_clay[index])
            elif self.veoStress_clay[index] < self.vefStress_clay[index] < self.preconsol_stress[index]:
                print(f'Sublayer {index+1}: OCOC Settlement Case')
                self.w_c[index] = (layer_thickness/(1 + self.e0)) * Cr * (log10(self.vefStress_clay[index]/self.veoStress_clay[index]))
            else:
                print(f'Sublayer {index+1}: OCNC Settlement Case')
                self.w_c[index] = (layer_thickness/(1 + self.e0)) * (((Cr * (log10(self.preconsol_stress[index]/self.veoStress_clay[index]))) + ((Cc * (log10(self.vefStress_clay[index]/self.preconsol_stress[index]))))))

        self.skempMu = (skempA + skempAlpha * (1-skempA))
        self.w_c_3d = self.skempMu * self.w_c
        self.settlement = np.sum(self.w_c_3d) + self.w_i


    def display_foundation(self):
        print(f'Unit System: {self.unit_system}')
        print(f'Foundation Short Side: {self.B}')
        print(f'Foundation Long Side: {self.L}')
        print(f'Foundation Depth: {self.D}')
        print(f'Foundation Thickness: {self.thickness}')
        print(f'Stress considered: {self.load}')
        print(f'Depth from bottom of Foundation to Stiff Layer: {self.H}')
        print(f'Unit Weight of Foundation: {self.gamma}')
        print('')

    def display_sand(self):
        print(f'Unit System: {self.unit_system}')
        print(f'Unit Weight of Soil: {self.gamma_soil}')
        print(f'Average N60: {self.avgN60}')
        print(f'Groundwater table depth: {self.gwt}')
        print(f'Excavation is present to depth {self.excavation}')
        print('')

    def display_clay(self):
        print(f'Unit System: {self.unit_system}')
        print(f'Unit Weight of Soil: {self.gamma_soil}')
        print(f'Groundwater table depth: {self.gwt}')
        print(f'Excavation is present to depth {self.excavation}')
        print(f'Su: {self.Su}')
        print(f'PI: {self.PI}')
        print(f'OCR: {self.OCR_clay}')
        print(f'K: {self.K}')
        print(f'e0: {self.e0}')
        print('')

    def display_factors_sand(self): # display the results
        print(f'Unit System: {self.unit_system}')
        print(f'Influence Depth: {self.influence_depth}')
        print(f'Layer Thickness Factor: {self.layer_thickness_factor}')
        print(f'Compressibility Index: {self.comp_index}')
        print(f'Shape Factor: {self.shape_factor}')
        print(f'Time Factor: {self.time_factor}')
        print(f'Preconsolidation Stress: {self.preconsol_stress}')
        print('')

    def display_factors_clay(self):
        print(f'Unit System: {self.unit_system}')
        print(f'Influence Depth: {self.influence_depth}')
        print(f'Initial Vertical Effective Stress: {self.veoStress_clay}')
        print(f'Preconsolidation Stress: {self.preconsol_stress}')
        print(f'Additional Vertical Effective Stress (via 2:1 Method): {self.two_to_one_delta}')
        print(f'Final Vertical Effective Stress: {self.vefStress_clay}')
        print(f'Settlement in Clay Layers (1D): {self.w_c}')
        print(f"Skempton's μ: {self.skempMu}")
        print('')

    def display_results_sand(self):
        print(f'Unit System: {self.unit_system}')
        print(f'Stress considered: {self.load}')
        print(f'Settlement: {self.settlement}')
        print('')

    def display_results_clay(self):
        print(f'Unit System: {self.unit_system}')
        print(f'Stress considered: {self.load}')
        print(f'Immediate Settlement: {self.w_i}')
        print(f'Settlement in Clay Layers (3D): {self.w_c_3d}')
        print(f'SUM(3D): {np.sum(self.w_c_3d)}')
        print(f'Immediate + SUM(3D): {self.settlement}')
        print('')

