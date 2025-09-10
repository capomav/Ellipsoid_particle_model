#!/usr/bin/env python
# coding: utf-8

import numpy as np
from system import *
from Quaternions import *
from Quaternions import Quaternion, rotation_about_axis
from particles import *
from non_bonded_interactions import *
from cell_structure import *


class Integrators:
   def __init__(self, time_step, system:System):
      self.dt = time_step
      self.system = system

   
   cnt = 0
   
   def step(self, quaternions, gamma_trans, gamma_rot):
      # default euler scheme is implemented

      noise_trans_mag = np.sqrt(2*self.system.KT*gamma_trans/ self.dt)
      noise_rot_mag = np.sqrt(2*self.system.KT*gamma_rot/ self.dt)
      noise_trans = noise_trans_mag * np.random.normal(0,1, size=np.shape(self.system.force))
      noise_rot =  noise_rot_mag * np.random.normal(0,1, size=np.shape(self.system.force))

      max_disp = 0.0

      for i in range(len(self.system.part)):
         
         self.system.force_calculations(self.system.verlet_list)
         self.system.torque_calculations(self.system.verlet_list)
          
         self.system.part[i].vel = (1/gamma_trans) * (self.system.force[i] + noise_trans[i]) 
         displacement = self.system.part[i].vel * self.dt
         self.system.part[i].pos += displacement

         max_disp = max(max_disp, float(np.linalg.norm(displacement)) )
         
         self.system.part[i].pos = self.system.part[i].pos - self.system.cell.box_l * np.round( self.system.part[i].pos / self.system.cell.box_l)  #  periodic boundary condition1

         self.system.part[i].ang_vel = (1/gamma_rot) * (self.system.torque[i] + noise_rot[i]) 

         dtheta = np.linalg.norm(self.system.part[i].ang_vel) * self.dt
          
         if(self.cnt == 0):
            u1_new = rotation_about_axis(dtheta,self.system.part[i].ang_vel, self.system.part[i].u1)
            u2_new = rotation_about_axis(dtheta,self.system.part[i].ang_vel, self.system.part[i].u2)         

         insta_ang_vel = self.system.part[i].ang_vel   
         q_ang_vel = Quaternion(0,insta_ang_vel[0],insta_ang_vel[1],insta_ang_vel[2])

         quaternions[i] += (0.5 *self.dt) * quaternions[i] * q_ang_vel

         q_u1 = Quaternion(0,self.system.part[i].u1[0],self.system.part[i].u1[1],self.system.part[i].u1[2])
         q_u2 = Quaternion(0,self.system.part[i].u2[0],self.system.part[i].u2[1],self.system.part[i].u2[2])

         q_u1_new = quaternions[i] * q_u1 * quaternions[i].conjugate()
         q_u2_new = quaternions[i] * q_u2 * quaternions[i].conjugate() 

         u1_new = q_u1_new.get_director()
         u2_new = q_u2_new.get_director()

         print(f"u1_new=  {u1_new}")
         print(f"u2_new=  {u2_new}")
         print("\n")
         #print(np.shape(u1_new))
         #print(u2_new)

         #u1_new = Quaternion.get_director(q_u1_new)
         #u2_new = Quaternion.get_director(q_u2_new)

         self.system.part[i].u1 = u1_new / np.linalg.norm(u1_new)
         self.system.part[i].u2 = u2_new / np.linalg.norm(u2_new)
      
      print("completed step = ", self.cnt)
      self.cnt += 1
      
      self.system.cell.update_verlet(max_disp) 

