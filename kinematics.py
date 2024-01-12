import numpy as np

class Component():
    def __init__(self,x,y,maj,min,pos,flux,date,mjd,year,component_number=-1):
        self.x = x
        self.y = y
        self.mjd = mjd
        self.maj = maj
        self.min = min
        self.pos = pos
        self.flux = flux
        self.date = date
        self.mjd = mjd
        self.year = year
        self.component_number=component_number

    def assign_component_number(self,number):
        self.component_number=number

class ComponentCollection():
    def __init__(self,components=[]):
        self.components=components
    def get_speed(self):
        pass



