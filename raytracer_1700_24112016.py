# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 15:43:16 2016

@author: R Carthigesan

Module implementing a 3-D optical ray tracer, which can be used to model the behaviour of some simple optical systems.
"""
from numpy import array

class Ray:
    """
    Class representing optical ray. Create instance at starting point point_init and direction dir_init using Ray(point_init, dir_init)
    """
    
    def __init__(self, dir_init=array([]), point_init=array([0.0, 0.0, 0.0])):
        """
        Initialises ray: if initial point not provided, starts at origin.
        """
        if dir_init.size != 3:
            raise Exception("Direction not 3-vector")
        if point_init.size != 3:
            raise Exception("Point not 3-vector")
        self.dir_init = dir_init
        self.dirlist = [dir_init]
        self.point_init = point_init
        self.pointlist = [point_init]
        
    def pcurrent(self):
        """
        Returns current point of ray.
        """
        return self.pointlist[-1]
        
    def dcurrent(self):
        """
        Returns current direction of ray.
        """
        return self.dirlist[-1]
        
    def append(self,d,p):
        """
        Adds new direction d and point p  to ray.
        """
        if d.size != 3:
            raise Exception("New direction not 3-vector")
        if p.size != 3:
            raise Exception("New point not 3-vector")
        self.dirlist.append(d)
        self.pointlist.append(p)
        
    def vertices(self):
        """
        Returns list of all points along ray.
        """
        return self.pointlist
        
class OpticalElement:
    
  def propagate_ray(self, ray):
    "propagate a ray through the optical element"
    raise NotImplementedError()
    
class SphericalRefraction(OpticalElement):
    """
    Class representing a spherical refracting surface
    """
    def __init__(self, z_0, curv, n_1, n_2, ap_rad):
        self.z_0 = z_0
        self.curv = curv
        self.n_1 = n_1
        self.n_2 = n_2
        self.ap_rad = ap_rad