#!/usr/bin/env python
# coding: utf-8

import numpy as np
from particles import *

class Cell_Structure:

    def __init__(self, box_length:float, r_cutoff:float, r_skin:float, part_num:int, num_types:int, particles:list[Particle]):
        self.box_l = box_length
        self.r_cut = r_cutoff
        self.r_skin = r_skin
        self.part_num = part_num
        self.num_type = num_types
        self.particles = np.array(particles) 

    def build_verlet_list(self, part_num):
        self.verlet_list = [ [] for _ in range(part_num)]

        for i in range(part_num):
            for j in range(i+1,part_num):

                rij = np.linalg.norm(self.particles[i].pos - self.particles[j].pos)
                
                rij = rij - self.box_l  * np.round(rij/self.box_l)   # minimum image convention
                
                if(rij < self.r_skin + self.r_cut):
                    self.verlet_list[i].append(j)
                    self.verlet_list[j].append(i)
        
        return self.verlet_list

    def update_verlet(self, max_disp = None):
        if(max_disp == None):
            max_disp = self.r_cut/2

        if(max_disp > self.r_skin):
            self.verlet_list = self.build_verlet_list(self.part_num)
            return self.verlet_list    
        
        return None


