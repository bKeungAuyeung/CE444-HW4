from math import log10

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

    def soil(self, gamma_soil, gwt, avgN60, excavation = False):
        self.gamma_soil = gamma_soil # match units specified by the foundation
        self.gwt = gwt
        self.avgN60 = avgN60
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

    def settlement_clay(self, time, depth_to_base, I0 = -1, I1 = -1, OCR = 1):
        # 

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

    def display_soil(self):
        print(f'Unit System: {self.unit_system}')
        print(f'Unit Weight of Soil: {self.gamma_soil}')
        print(f'Average N60: {self.avgN60}')
        print(f'Unit Weight of Water: {self.gamma_w}')
        print(f'Groundwater table depth: {self.gwt}')
        print(f'Excavation is present to depth {self.excavation}')
        print('')

    def display_factors(self): # display the results
        print(f'Unit System: {self.unit_system}')
        print(f'Influence Depth: {self.influence_depth}')
        print(f'Layer Thickness Factor: {self.layer_thickness_factor}')
        print(f'Compressibility Index: {self.comp_index}')
        print(f'Shape Factor: {self.shape_factor}')
        print(f'Time Factor: {self.time_factor}')
        print(f'Preconsolidation Stress: {self.preconsol_stress}')
        print('')


    def display_results(self):
        print(f'Unit System: {self.unit_system}')
        print(f'Stress considered: {self.load}')
        print(f'Settlement: {self.settlement}')
        print('')

