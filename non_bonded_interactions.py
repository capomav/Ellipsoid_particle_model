#!/usr/bin/env python
# coding: utf-8

import numpy as np
from particles import *

class Non_bonded_Interactions:
        
        class LJ:
            def __init__(self, epsilon:float= 1.0, sigma:float = 1.0, r_cut:float = 1.0, shift:float = 0.0, offset:float = 0.0):
                self.eps = epsilon
                self.sig = sigma
                self.r_cut = r_cut
                self.shift = shift
                self.offset = offset
            

            def lj_potential(self, dist):
                frac = self.sig/(dist - self.offset)
                if(self.shift == "auto"):
                    frac_c = self.r_cut/(dist - self.offset) 
                    sft = - frac_c**12 + frac_c**6
                else:
                    sft = self.shift
                lj_pot = 4*self.eps * ((frac**12 - frac**6) + sft)
                return lj_pot
            

            def lj_force(self, dist):
                frac = self.sig/(dist-self.offset)
                force = 48*self.sig*(frac**12 * (1/dist) - 0.5* frac**6 * (1/dist) )
                sft = 48*self.sig*( (self.sig**12 / self.r_cut**13)  - 0.5*(self.sig **6 / self.r_cut**7) )
                force -= sft
                return force
            
        
        class WCA:
            def __init__(self, epsilon:float= 1.0, sigma:float = 1.0, r_cut:float = 1.0):
                self.eps = epsilon
                self.sig = sigma
                self.r_cut = r_cut

            def wca_potential(self,dist):
                if(dist < self.r_cut):
                    frac6 = (self.sig/dist)**6
                    return 4.0*self.eps*(frac6**2 - frac6 + 0.25)
                else:
                    return 0.0
                
            def wca_force(self,dist):
                if(dist < self.r_cut):
                    frac6 = (self.sig/dist)**6
                    return 48*self.eps*frac6*(frac6 - 0.5) / (dist*dist)
                return 0.0            
        
        class Gay_Berne:
            def __init__(self, epsilon_0:float = 1.0, sigma_0:float = 1.0, kappa:float = 1.0 , kappa_prime:float = 1.0, mu : float = 1.0, nu:float = 1.0):
                self.eps0 = epsilon_0
                self.sig0 = sigma_0
                self.kap = kappa
                self.kappr = kappa_prime
                self.mu = mu
                self.nu = nu 

                # return the chi and chi prime respectively
                self.chi = (self.kap**2 - 1)/(self.kap **2 + 1)
                self.chi_prime = (1 - self.kappr**(1/self.mu) )/(1 + self.kappr ** (1/self.mu))                

            def sigma(self,u1, u2, r12):
                return self.sig0 *( 1 - 0.5* self.chi * (  (np.dot(r12,u1) + np.dot(r12,u2))**2 / (1 + self.chi * np.dot(u1,u2)) +  (np.dot(r12,u1) - np.dot(r12,u2))**2 / (1 - self.chi * np.dot(u1,u2))  ))** -0.5

            def epsilon(self, u1, u2, r12):

                eps1 = (1 - ( self.chi* np.dot(u1,u2) )**2 )**-0.5
                eps2 =  1 - 0.5* self.chi_prime * (  (np.dot(r12,u1) + np.dot(r12,u2))**2 / (1 + self.chi_prime * np.dot(u1,u2)) +  (np.dot(r12,u1) - np.dot(r12,u2))**2 / (1 - self.chi_prime * np.dot(u1,u2)))                
                eps = self.eps0 * eps1**self.nu * eps2 ** self.mu
                return eps
            
                  
            def gb_potential(self, u1, u2, r12):
                sig_cal = ( self.sig0 / (np.linalg.norm(r12) - self.sigma(u1,u2,r12) + self.sig0) )**6
                pot = 4*self.epsilon(u1, u2, r12) * ( sig_cal**2 - sig_cal )
                return pot
            
        

        class exp_potential:
            def __init__(self, part1:Particle ,part2:Particle, lambda_ratios: list = [1/3,1/3,1/3], b:float = 1.0):

                
                self.part1 = part1
                self.part2 = part2
                self.r12 = part2.pos - part1.pos 
                self.lambda_ratios = lambda_ratios
                self.b = b

            '''
            p1 = np.array(self.part1.u1)
            p2 = np.array(self.part2.u1)
            q1 = np.array(self.part1.u2)
            q2 = np.array(self.part2.u2)
            b = self.b
            r12 = np.array(self.r12)
            '''

            def S(self,r12,p1,q1,p2,q2):

                S1 = np.dot(p1,r12) * np.dot(p2,r12) / np.linalg.norm(r12)**2 
                S2 = np.dot(p1,q1) * np.dot(p2,q2)
                S3 = np.dot(q1,r12) * np.dot(q2,r12) / np.linalg.norm(r12)**2 

                S = [S1,S2,S3]
                return S     
            
            def potential(self,r12, p1,q1,p2,q2,b):
                V = np.exp(-np.linalg.norm(r12)) - self.S(r12,p1,q1,p2,q2)*np.exp(-np.linalg.norm(r12) /b )
                return V 


            def force_1_due_2(self):
                # total force from elife paper
                
                p1 = np.array(self.part1.u1)
                p2 = np.array(self.part2.u1)
                q1 = np.array(self.part1.u2)
                q2 = np.array(self.part2.u2)
                b = self.b
                r12 = np.array(self.r12)
                r12_mag = np.linalg.norm(r12)
                r12_u = r12/r12_mag 


                gam =  np.exp(-r12_mag*(b-1)/b) - np.sum(self.S(r12,p1,q1,p2,q2)) /b  + 2 * (1/r12_mag) * ( self.lambda_ratios[0]*np.dot(r12_u,p1)* np.dot(r12_u,p2) + self.lambda_ratios[2]*np.dot(r12_u,q1)*np.dot(r12_u,q2) )

                force_tot = np.exp(-r12_mag/b) * ( gam * r12_u- self.lambda_ratios[0] * (1/r12_mag) *  (  np.dot(r12_u,p2)*p1 + np.dot(r12_u,p1)*p2 ) - self.lambda_ratios[2]* (1/np.linalg.norm(r12) ) * (np.dot(r12_u,q2)*q1 + np.dot(r12_u,q2)*q1) )      

                '''
                def f_S1(self, ):    
                    return f_S1
                
                def f_S2(self, ):
                    return f_S2
                
                def f_S3(self, ):
                    return f_S3 
                '''             

                return force_tot
            

            def torque_1_due_2(self):

                p1 = np.array(self.part1.u1)
                p2 = np.array(self.part2.u1)
                q1 = np.array(self.part1.u2)
                q2 = np.array(self.part2.u2)
                b = self.b
                r12 = np.array(self.r12)
                r12_mag = np.linalg.norm(r12)
                r12_u = r12/r12_mag
                
                S = self.S(r12,p1,q1,p2,q2)

                dV_dp1 = np.exp(r12_mag/b) * (  self.lambda_ratios[0] * ( S[0]* p1 - p2 + np.dot(r12_u,p2)*r12_u ) + self.lambda_ratios[1] * ( S[1]*p1 - np.dot(q1,q2)*p2 + np.dot(q1,p2)*q2 )  )
                
                dV_dq1 = np.exp(r12_mag/b) * (  self.lambda_ratios[1] * ( S[1] * q1 - np.dot(p1,p2)*q2 + np.dot(p1,q2)*p2 ) + self.lambda_ratios[2]* ( S[2] * q1 - q2 + np.dot(r12_u,q2) * r12_u )  ) 

                torque_p =  - np.cross(p1, dV_dp1)
            
                torque_q = - np.cross(q1, dV_dq1)

                torque_1_due_2_tot = torque_p + torque_q

                return torque_1_due_2_tot
            
        '''
        if(len(forces)!= 0):
            if("LJ" in forces):
                lj = LJ()
    
            if("WCA" in forces):
                wca = WCA()
    
            if("Gay-Berne" in forces):
                gb = Gay_Berne()
                
            else:
                raise ValueError("The argument 'name' doesn't match pre-defined interactions or provide a incorrect input")
        else:
            print("\n No non-bonded interactions present in the system \n")
        '''


