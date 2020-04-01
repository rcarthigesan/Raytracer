# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 15:43:16 2016

@author: R Carthigesan

Module implementing a 3-D optical ray tracer, which can be used to model the behaviour of some simple optical systems.
"""
import numpy as np
import numpy.linalg as npl


class Ray:
    """
    Class representing optical ray. Create instance at starting point point_init and direction dir_init using Ray(point_init, dir_init)
    """
    
    def __init__(self, dir_init, point_init=np.array([0.0, 0.0, 0.0])):
        """
        Initialises ray: if initial point not provided, starts at origin. Starting point and direction are given as NumPy arrays of size 3, representing Cartesian 3 vectors.
        """
        if dir_init.size != 3:
            raise Exception("Direction not 3-vector")
        elif point_init.size != 3:
            raise Exception("Point not 3-vector")
        else:
            self.dir_init = dir_init.astype(np.float)
            self.dirlist = [self.dir_init]
            self.point_init = point_init.astype(np.float)
            self.pointlist = [self.point_init]
        
    def pcurrent(self):
        """Returns current point of ray."""
        return self.pointlist[-1]
        
    def dcurrent(self):
        """Returns current direction of ray."""
        return self.dirlist[-1]
        
    def append(self,d,p):
        """Adds new direction d and point p  to ray."""
        if d.size != 3:
            raise Exception("New direction not 3-vector")
        elif p.size != 3:
            raise Exception("New point not 3-vector")
        else:
            self.dirlist.append(d.astype(np.float))
            self.pointlist.append(p.astype(np.float))
        
    def vertices(self):
        """Returns list of all points along ray."""
        return self.pointlist

        
class OpticalElement:
    
  def propagate_ray(self, ray):
    """Propagate a ray through the optical element."""
    raise NotImplementedError()

    
class SphericalRefraction(OpticalElement):
    """
    Class representing a spherical refracting surface centred on the optical axis. z_0 is the intercept of the surface with the axis; curv, the curvature, is the reciprocal of the radius of curvature; n_1 and n_2 are the refractive indices either side of the surface; ap_rad is the aperture radius - the maximum extent of the surface from the optical axis.
    """
    def __init__(self, z_0, curv, n_1, n_2, ap_rad):
        self.z_0 = z_0
        self.curv = curv
        self.n_1 = n_1
        self.n_2 = n_2
        self.ap_rad = ap_rad
        
    def intercept(self,ray):
        """Calculates the first valid intercept of a ray with the spherical surface."""
        P=ray.pcurrent()
        R=1.0/self.curv
        centrecurve = (np.array([0.0,0.0,self.z_0 + R])).astype(np.float) #centre of curvature
        k_hat=(ray.dcurrent())/(npl.norm(ray.dcurrent()))
        r = np.subtract(P, centrecurve)
        Q = np.add(P, np.dot(r, k_hat))
        if npl.norm(np.subtract(Q, centrecurve)) > abs(R):
            return None
        else:
            l_plus = -1 * (np.dot(r, k_hat)) + np.sqrt((np.dot(r, k_hat))**2 - ((npl.norm(r))**2 - R**2))
            l_minus = -1 * (np.dot(r, k_hat)) - np.sqrt((np.dot(r, k_hat))**2 - ((npl.norm(r))**2 - R**2))
            intersect_plus = np.add(P, l_plus * k_hat)
            intersect_minus = np.add(P, l_minus * k_hat)
            return intersect_plus, intersect_minus
            if R>0 and P[-1]<self.z_0:
                if np.sqrt((intersect_minus[0])**2 + (intersect_minus[1])**2) > abs(R):
                    return None
                else:
                    return intersect_minus
            elif R<0 and P[-1]>self.z_0:
                if np.sqrt((intersect_minus[0])**2 + (intersect_minus[1])**2) > abs(R):
                    return None
                else:
                    return intersect_minus
            elif R>0 and P[-1]>(-self.z_0):
                if np.sqrt((intersect_plus[0])**2 + (intersect_plus[1])**2) > abs(R):
                    return None
                else:
                    return intersect_plus
            elif R<0 and P[-1]<(-self.z_0):
                if np.sqrt((intersect_plus[0])**2 + (intersect_plus[1])**2) > abs(R):
                    return None
                else:
                    return intersect_plus
            elif np.sign(R) == np.sign(k_hat[-1]):
                return None
            
            
        
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        