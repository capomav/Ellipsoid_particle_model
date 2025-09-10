#!/usr/bin/env python
# coding: utf-8

import numpy as np

# ## Defining the quaternions and its algebraic rules
# 
# A quaternion is an expression of the form  :<br>
# a + b <b>i</b>  + c <b>j</b> + d <b>k</b> <br>
# where a, b, c, d, are real numbers, and i, j, k, are symbols that can be interpreted as unit-vectors pointing along the three spatial axes. 
# 
# Set of quaternions is a 4-dimensional vector space over the real number with {1,<b>i</b>,<b>j</b>,<b>k</b>} as basis vector. It can also be represented as [a,<b>v</b>] where a is scalar quantity and <b>v</b> is a 3-dimensional vector.
# 
# ### Multiplication of basis elements : 
# 
# Multiplication of the basis elemnt is defined as follows : <br>
# 
# 1 is a multiplicative identity <br>
# $ 1 \hat{i} = \hat{i} 1 = \hat{i} $, same for other components ... <br>
# 
# product of other basis elements are <br>
# $ \hat{i}^{2} = \hat{j}^{2} = \hat{k}^{2} = -1 $ <br>
# $ \hat{i}\hat{j} = -\hat{j}\hat{i} = \hat{k} $, $\hat{j}\hat{k} = -\hat{k}\hat{j} =  \hat{i}$,   $\hat{k}\hat{i} = - \hat{i}\hat{k} = \hat{j} $ <br>
# 
# Combining these rules we get, <br>
# $ \hat{i}\hat{j}\hat{k} = -1 $
# 
# ### Scalar and Hamilton product
# 
# 
# Multiplying a quaternion with a scalar <br>
# 
# $ \lambda( a + b \hat{i} + c \hat{j} + d \hat{k}) = \lambda a + (\lambda b)\hat{i}  + (\lambda c)\hat{j} + (\lambda d)\hat{k} $ <br>
# 
# 
# Hamilton product for two quaternions is defined as <br>
# 
# For two quaternions 
# 
# $ q_{1} =  a_{1} + b_{1} \hat{i} + c_{1} \hat{j} + d_{1} \hat{k}$ and $q_{2} = a_{2} + b_{2} \hat{i} + c_{2} \hat{j} + d_{2} \hat{k} $ <br>
# 
# $ q_{1}q_{2} $ is determined by the product of their basis elements and distributive law  which gives the following expression<br>
# 


class Quaternion :

    def __init__(self, a,b,c,d) : 
        self.a = float(a)
        self.b = float(b)
        self.c = float(c)
        self.d = float(d)

    
    def __str__(self):
        return f"{self.a} + {self.b}i + {self.c}j + {self.d}k"

    
    def __add__(self,other):
        if not isinstance(other, Quaternion):
            raise TypeError(f"Unsupported operand type(s) for +: 'Quaternion' and '{type(other).__name__}'")
        return Quaternion(self.a + other.a, self.b + other.b, self.c + other.c, self.d + other.d) 


    
    def __sub__(self,other):
        if not isinstance(other, Quaternion):
            raise TypeError(f"Unsupported operand type(s) for -: 'Quaternion' and '{type(other).__name__}'")
        return Quaternion(self.a - other.a, self.b - other.b, self.c - other.c, self.d - other.d) 

    
    def __mul__(self,other):
        if not isinstance(other, Quaternion):
            if(type(other).__name__ == float or int):
                return Quaternion(self.a * other, self.b * other, self.c * other, self.d * other)
            raise TypeError(f"Unsupported operand type(s) for *: 'Quaternion' and '{type(other).__name__}'")
        a = self.a * other.a - (self.b * other.b + self.c * other.c + self.d * other.d)
        b = self.a * other.b + self.b * other.a + self.c * other.d - self.d * other.c
        c = self.a * other.c - self.b * other.d + self.c * other.a + self.d * other.b
        d = self.a * other.d + self.b * other.c - self.c * other.b + self.d * other.a

        return Quaternion(a,b,c,d)
    
    def __rmul__(self,other):
        if not isinstance(other, Quaternion):
            if(type(other).__name__ == float or int):
                return Quaternion(self.a * other, self.b * other, self.c * other, self.d * other)
            raise TypeError(f"Unsupported operand type(s) for *: 'Quaternion' and '{type(other).__name__}'")
        a = self.a * other.a - (self.b * other.b + self.c * other.c + self.d * other.d)
        b = self.a * other.b + self.b * other.a + self.c * other.d - self.d * other.c
        c = self.a * other.c - self.b * other.d + self.c * other.a + self.d * other.b
        d = self.a * other.d + self.b * other.c - self.c * other.b + self.d * other.a

        return Quaternion(a,b,c,d)

    
    def conjugate(self):
        return Quaternion(self.a, -self.b, -self.c, -self.d)

    def norm(self):
        return np.sqrt(self.a*self.a + self.b*self.b + self.c*self.c +self.d*self.d)

    
    def __truediv__(self,other):
        if not isinstance(other, Quaternion):
            if(type(other).__name__ == float or int):
                return Quaternion(self.a / other, self.b / other, self.c / other, self.d / other)
            raise ErrorType(f"Unsupported operand type(s) for *: 'Quaternion' and '{type(other).__name__}'")
       
        sq_norm = np.square(self.norm())
        prod = self * self.conjugate()
        
        return Quaternion (prod.a/sq_norm, prod.b/sq_norm, prod.c/sq_norm, prod.d/sq_norm)
        
        
    def __iadd__(self, other):
        if not isinstance(other, Quaternion):
            raise TypeError(f"unsupported operand type(s) for +=: 'Quaternion' and '{type(other).__name__}'")
        self.a += other.a
        self.b += other.b
        self.c += other.c
        self.d += other.d
        return self

    
    def __isub__(self, other):
        if not isinstance(other, Quaternion):
            raise TypeError(f"unsupported operand type(s) for -=: 'Quaternion' and '{type(other).__name__}'")
        self.a -= other.a
        self.b -= other.b
        self.c -= other.c
        self.d -= other.d
        return self

    
    def __imul__(self, other):
        if not isinstance(other, Quaternion):
            raise TypeError(f"unsupported operand type(s) for *=: 'Quaternion' and '{type(other).__name__}'")
        self.a = self.a * other.a - self.b * other.b - self.c * other.c - self.d * other.d
        self.b = self.a * other.b + self.b * other.a + self.c * other.d - self.d * other.c
        self.c = self.a * other.c - self.b * other.d + self.b * other.a + self.d * other.b
        self.d = self.a * other.d + self.b * other.c - self.c * other.b + self.d * other.a
        return self

    
    def __itruediv__(self, other):
        if not isinstance(other, Quaternion):
            if(type(other).__name__ == float or int):
                self.a = self.a / other
                self.b = self.b / other
                self.c = self.c / other
                self.d = self.d / other
                return self
            raise TypeError(f"unsupported operand type(s) for /=: 'Quaternion' and '{type(other).__name__}'")
        
        sq_norm = np.square(self.norm())
        prod = self * self.conjugate()
        
        self.a = prod.a / sq_norm
        self.b = prod.b / sq_norm
        self.c = prod.c / sq_norm
        self.d = prod.d / sq_norm
        
        return self


    def qdot(self, q1, q2):
        if(isinstance(q1, Quaternion) and isinstance(q2, Quaternion)):
            qdot = q1.a * q2.a + q1.b * q2.b + q1.c * q2.c + q1.d * q2.d

            return qdot
        else:
            raise TypeError("Provide quaternion inputs")



    def get_director(self):
        # returns the director vector i.e axis of rotation 

        q0 = self.a
        q1 = self.b
        q2 = self.c
        q3 = self.d

        ux = 2*(q1*q3 + q0*q2)
        uy = 2*(q2*q3 - q0*q1)
        uz = q0*q0 - q1*q1 * q2*q2 + q3*q3

        return [ux,uy,uz]


    def get_rotation_mat(self):
        # returns the matrix for rotating a 3d vector 


        q0 = self.a
        q1 = self.b
        q2 = self.c
        q3 = self.d
        
        rot_mat = np.zeros((3,3))
        rot_mat[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3
        rot_mat[0][1] = 2*(q1*q2 + q0*q3)
        rot_mat[0][2] = 2*(q1*q3 - q0*q2)

        rot_mat[1][0] = 2*(q1*q2 - q0*q3)
        rot_mat[1][1] = q0*q0 + q2*q2 - q1*q1 - q3*q3 
        rot_mat[1][2] = 2*(q2*q3 + q0*q1)

        rot_mat[2][0] = 2*(q1*q3 + q0*q2)
        rot_mat[2][1] = 2*(q2*q3 - q0*q1)
        rot_mat[2][2] = q0*q0 -q1*q1 - q2*q2 + q3*q3
        
        return rot_mat


def rotation_about_axis(theta, axis, vec):
    # rotates a vector (vec) about the director vector (axis) by an angle (theta)

    if(len(axis) == 3 and len(vec) == 3):
        axis = axis/np.linalg.norm(axis) 
        a = np.cos(theta/2)
        b = np.sin(theta/2)*axis[0]
        c = np.sin(theta/2)*axis[1]
        d = np.sin(theta/2)*axis[2]

        q_rot = Quaternion(a,b,c,d)
        q_rot = q_rot * (1/q_rot.norm())
        q_axis = Quaternion(0, axis[0], axis[1], axis[2])
        q_vec = Quaternion(0, vec[0], vec[1], vec[2])

        q_vec_rot = q_rot * q_vec * q_rot.conjugate()

        vec_rot = np.array([q_vec_rot.b, q_vec_rot.c, q_vec_rot.d])
        
        #print(q_rot.get_rotation_mat())
        
        return vec_rot
    
    else: 
        if(len(axis)!=3) : raise ValueError(f"size of 'axis' should be (3,) but was given : {len(axis)}")
        if(len(vec)!=3) : raise ValueError(f"size of 'vec' should be (3,) but was given : {len(vec)}")
        raise    


def quat_initial(n_cap):
    theta_0 = np.arccos(np.dot( np.array([0.,0.,1.]), n_cap))
    w_0 = np.cross(np.array([0.,0.,1.]), n_cap)

    q_0 = np.cos(theta_0/2)
    q_1 = np.sin(theta_0/2)* w_0[0]
    q_2 = np.sin(theta_0/2)* w_0[1]
    q_3 = np.sin(theta_0/2)* w_0[2] 
    
    quat_0 = Quaternion(q_0, q_1, q_2,q_3)

    return quat_0