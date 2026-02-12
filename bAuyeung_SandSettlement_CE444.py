from math import log10

class Foundation:
    def __init__(self, B, L, D, H, thickness, load, gamma, metric = True):
        self.B = B
        self.L = L
        self.D = D
        self.H = H
        self.thickness = thickness
        self.load = load 
        self.gamma = gamma
        if metric == False: # if the system units are in ft
            self.Pa = 2117 # psf
            self.refL = 3.28 # ft
            self.gamma_w = 62.4 # pcf
        else:
            self.Pa = 100 # kPa
            self.refL = 1 # m
            self.gamma_w = 9.81 # kN/m^3


    def soil(self, gamma_soil, gwt, layer_depth, avgN60, excavation = False):
        self.gamma_soil = gamma_soil # match units specified by the foundation
        self.gwt = gwt
        self.layer_depth = layer_depth
        self.avgN60 = avgN60
        if excavation == True:
            self.excavation = excavation
        else:
            self.excavation = -1 

    def influence_depth(self):
        self.influence_depth = self.refL * ((self.B / self.refL)**0.79)

    def settlement(self, time, deep = False, dynamic = False):
        if deep == True:
            self.layer_thickness_factor = 1
        else:
            H_infl_ratio = self.H / self.influence_depth
            self.layer_thickness_factor = H_infl_ratio * (2 - H_infl_ratio)
        
        self.comp_index = 1.71 / (self.avgN60**1.4)

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

        



