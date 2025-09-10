import numpy as np
from system import *
from Quaternions import *
from particles import *
from non_bonded_interactions import *
from cell_structure import *
from integrator import *

part_num = 2
box_l = 15
r_skin = 2.5
r_cut = 3.5
num_types = 1
time_step = 0.01
N_steps = 10

pos = [ np.random.random(size= (3,1)) * box_l for _ in range(part_num) ]
#u_random1= np.array([t /np.linalg.norm(t) for _ in range(part_num) for t in [np.random.random(size= (3,1))] ])
#u_random2= [t /np.linalg.norm(t) for _ in range(part_num) for t in [np.random.random(size= (3,1))] ]

t = np.random.random(size= (part_num,3))
u_random1 = t / np.linalg.norm(t, axis=1, keepdims=True)

t2 = np.random.random(size= (part_num,3))
u_random2 = t2 / np.linalg.norm(t2, axis=1, keepdims=True)

print(f"pos = {pos}")
print(f"u1 = {u_random1}")
print(f"u1 = {u_random2}\n\n")


particles = [Particle(pos = pos[i].flatten(), u1 = np.asarray(u_random1[i]), u2= np.asarray(u_random2[i]) ) for i in range(len(pos))]

interact = Non_bonded_Interactions()
cell = Cell_Structure(box_length = box_l, r_cutoff = r_cut, r_skin = r_skin, part_num = part_num, num_types = num_types, particles= particles)
system = System(particles, interact, cell, 1.0)
integ = Integrators(time_step, system)

quaternions = [quat_initial(u_random2[i]) for i in range(part_num)]

for _ in  range(N_steps):
    integ.step(quaternions, gamma_trans = 1.0,gamma_rot = 1.0)
