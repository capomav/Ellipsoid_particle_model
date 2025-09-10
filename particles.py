#!/usr/bin/env python
# coding: utf-8
import numpy as np

class Particle :
    # Particle class with necessary arguments for the particle properties are initilaised here.
    # id,particle_type(species) position, velocity, angular velocity and internal orientational axes (u1,u2). 

    counter = 0    #  counts the total number of particles created. If a particle is deleted this counter decreases. 
    cnt = 0 # counts the total times particle instance is created to number the ids. If the particle is deleted this counter doesn't decrease.

    def __init__(self,  id = None, par_type: int = 0, pos: np.ndarray = np.random.rand(1,3), vel : np.ndarray = np.array([0.0,0.0,0.0]), \
                 ang_vel : np.ndarray = np.array([0,0,0.0]), mass : float = 1.0, u1 : np.ndarray = np.array([1.0,0,0]), u2 : np.ndarray = np.array([0,1.0,0]), \
                ext_f: list[float] = [0.0,0.0,0.0], ext_torque: list[float] = [0.0,0.0,0.0], cnt: int = cnt, counter: int = counter) :
        
        self.id = id if id is not None else Particle.counter
        self.pos = np.random.rand(3) if pos is None else np.array(pos).flatten()
        self.vel = np.zeros(3) if vel is None else np.array(vel).flatten()
        self.ang_vel = np.zeros(3) if ang_vel is None else np.array(ang_vel).flatten()
        self.u1 = np.array([1.0, 0.0, 0.0]) if u1 is None else np.array(u1).flatten()
        self.u2 = np.array([0.0, 1.0, 0.0]) if u2 is None else np.array(u2).flatten()
        self.ext_f = np.zeros(3) if ext_f is None else np.array(ext_f).flatten()
        self.ext_torque = np.zeros(3) if ext_torque is None else np.array(ext_torque).flatten()

        
        '''
        if(isinstance(id,int)):
            self.id = id
        else:
            self.id = cnt    
        
        if(isinstance(par_type,int)):
            self.par_type = par_type

        if(isinstance(mass,int)):
            self.m = mass
        
        if(len(pos) == 3):
            self.pos = pos

        if(len(vel) == 3 and isinstance(vel[:], float)):
            self.vel = vel   

        if(len(ang_vel) == 3 and isinstance(vel[:], float)):
            self.ang_vel = ang_vel    
        
        if(len(u1) == 3):
            self.u1 = u1
        
        if(len(u2) == 3):
            self.u2 = u2
        '''

        counter += 1
        cnt += 1 

    # defining the getter functions to retrieve the values of particle instance 
      
    def get_particle_count(self):
        # return the total number of particles created 
        return Particle.counter

    
    def get_id(self):
        # return the particle id of the current particle
        return self.id
    

    def get_particle_id(self, part):
        # return the particle id of the particle passed as an argument
        return part.get_id()
    

    def get_particle_velocity(self):
        # return the velocity array of the particle that calls the function
        return self.vel 
    
    
    def get_particle_pos(self):
        # returns the postion of the particle that calls the function
        return self.pos
    
    
    def __del__(self):
        # deletes the particle and reduces the total count
        if(Particle.counter >1):
            Particle.counter -= 1

