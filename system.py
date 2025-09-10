#!/usr/bin/env python
# coding: utf-8

import numpy as np
from particles import *
from non_bonded_interactions import *
from cell_structure import *


class System:
    def __init__(self, particles:list[Particle], Interactions:Non_bonded_Interactions, cell_structure:Cell_Structure, KT: float = 1.0):
        self.part = particles
        self.interactions = Interactions
        self.cell = cell_structure
        self.KT = KT

        self.verlet_list = self.cell.build_verlet_list(self.cell.part_num)

        self.force = np.zeros((len(self.part),3))
        self.torque = np.zeros((len(self.part),3))


    def calculate_force_ij(self, p1, p2, lambda_ratios, b):
        #  define the forces due to non-bonded interactions, noise, external field and activity
        exp_inter = self.interactions.exp_potential(p1,p2,lambda_ratios, b)
        fij = exp_inter.force_1_due_2()

        print(f"fij = {fij}")

        return fij


    def calculate_torque_ij(self,p1,p2,lambda_ratios, b):
        # define the torque due to non-bonded interactions, noise, external field and activity
        exp_inter = self.interactions.exp_potential(p1,p2,lambda_ratios, b)
        tij = exp_inter.torque_1_due_2()
        print(f"tij = {tij}")


        return tij


    def force_calculations(self, verlet_list):
        for i in range(len(verlet_list)):
            for j in verlet_list[i]:
                force_ij = self.calculate_force_ij(self.part[i], self.part[j], lambda_ratios=[1/3, 1/3, 1/3], b=1.0)
                self.force[i] += force_ij + self.part[i].ext_f 


    def torque_calculations(self, verlet_list):
        for i in range(len(verlet_list)):
            for j in verlet_list[i]:
                torque_ij = self.calculate_torque_ij(self.part[i], self.part[j],lambda_ratios=[1/3, 1/3, 1/3], b=1.0)
                self.torque[i] += torque_ij + self.part[i].ext_torque
        


## TODO : Compute the stresses at time of calculation  for forces
