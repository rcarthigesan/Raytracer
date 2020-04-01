# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 15:43:16 2016

@author: R Carthigesan

Module implementing a 3-D optical ray tracer, which can be used to model the behaviour of some simple optical systems.
"""
import numpy as np
import numpy.linalg as npl
import genpolar as gp
import matplotlib.pyplot as pl


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
    Class representing a spherical refracting surface centred on the optical axis. z_0 is the intercept of the surface with the axis; curv, the curvature, is the reciprocal of the radius of curvature; n_1 and n_2 are the refractive indices inside and outside the surface; ap_rad is the aperture radius - the maximum extent of the surface from the optical axis.
    """
    def __init__(self, z_0, curv, n_1, n_2, ap_rad):
        self.z_0 = z_0
        self.curv = curv
        self.n_1 = n_1
        self.n_2 = n_2
        self.ap_rad = ap_rad
        
    def intercept(self,ray):
        """Calculates the first valid intercept of a ray with the spherical surface."""
        P = ray.pcurrent()
        k_hat=(ray.dcurrent())/(npl.norm(ray.dcurrent()))
        if self.curv == 0:
            if k_hat[-1] == 0 and P[-1] != self.z_0:
                return None
            elif P[-1] < self.z_0:
                if np.sign(k_hat[-1]) > 0:
                    Lambda = (self.z_0-P[2])/(k_hat[2])
                    intersect = np.add(P, Lambda*k_hat)
                    if np.sqrt((intersect[0])**2 + (intersect[1])**2) > self.ap_rad: 
                        return None
                    else:
                        return intersect
                else:
                    return None
            elif P[-1] > self.z_0:
                if np.sign(k_hat[-1]) < 0:
                    Lambda = (self.z_0-P[2])/(k_hat[2])
                    intersect = np.add(P, Lambda*k_hat)
                    if np.sqrt((intersect[0])**2 + (intersect[1])**2) > self.ap_rad: 
                        return None
                    else:
                        return intersect
                else:
                    return None
            else:
                return P
        else:
            R = 1.0/self.curv
            centrecurve = (np.array([0.0,0.0,self.z_0 + R])).astype(np.float)
            r = np.subtract(P, centrecurve)
            if ((np.dot(r, k_hat))**2 - ((npl.norm(r))**2 - R**2)) < 0:
                return None
            else:
                l_plus = -1 * (np.dot(r, k_hat)) + np.sqrt((np.dot(r, k_hat))**2 - ((npl.norm(r))**2 - R**2))
                l_minus = -1 * (np.dot(r, k_hat)) - np.sqrt((np.dot(r, k_hat))**2 - ((npl.norm(r))**2 - R**2))
                intersect_plus = np.add(P, l_plus * k_hat)
                intersect_minus = np.add(P, l_minus * k_hat)
                if np.sign(k_hat[-1]) == np.sign(R) and npl.norm(r) > R:
                    if np.sqrt((intersect_minus[0])**2 + (intersect_minus[1])**2) > self.ap_rad:
                        return None
                    else:
                        return intersect_minus
                elif np.sign(k_hat[-1]) != np.sign(R) and npl.norm(r) > R:
                    if np.sqrt((intersect_plus[0])**2 + (intersect_plus[1])**2) > self.ap_rad:
                        return None
                    else:
                        return intersect_plus
                elif np.sign(R) == np.sign(k_hat[-1]):
                    return None
                else:
                    if np.sqrt((intersect_minus[0])**2 + (intersect_minus[1])**2) > self.ap_rad:
                        return None
                    else:
                        return intersect_minus

    def refract(self, ray, intercept):
        """Calculates new direction vector of ray due to refraction at lens surface"""
        intercept = intercept.astype(np.float)
        r_i = (np.subtract(intercept, ray.pcurrent())) / (npl.norm(np.subtract(intercept, ray.pcurrent())))
        if self.curv != 0:
            centrecurve = (np.array([0.0,0.0,self.z_0 + (1/self.curv)])).astype(np.float)
            normal_out = (np.subtract(intercept, centrecurve)) / (npl.norm(np.subtract(intercept, centrecurve)))
            normal_in = (np.subtract(centrecurve, intercept)) / (npl.norm(np.subtract(centrecurve, intercept)))
            if np.dot(normal_out, r_i) < 0:
                normal = normal_out
                c_1 = np.dot((-1*normal), r_i)
                eta = (self.n_1/self.n_2)
            else:
                normal = normal_in
                c_1 = np.dot((-1*normal), r_i)
                eta = (self.n_2/self.n_1)
            if np.sqrt(1 - (eta**2) * (1 - (c_1**2))) < 0: #total internal reflection
                return None
            else:
                c_2 = np.sqrt(1 - (eta**2) * (1 - (c_1**2)))
                r_r = eta*r_i + (eta * c_1 - c_2)*normal
                return r_r
        else:
            normal_left = np.array([0,0,-1])
            normal_right = np.array([0,0,1])
            if np.dot((-1*normal_left), r_i) > 0:
                normal = normal_left
                c_1 = np.dot((-1*normal), r_i)
                eta = (self.n_1/self.n_2)
            else:
                normal = normal_right
                c_1 = np.dot((-1*normal), r_i)
                eta = (self.n_2/self.n_1)
            if np.sqrt(1 - (eta**2) * (1 - (c_1**2))) < 0: #total internal reflection
                return None
            else:
                c_2 = np.sqrt(1 - (eta**2) * (1 - (c_1**2)))
                r_r = eta*r_i + (eta * c_1 - c_2)*normal
                return r_r
        
    def propagate_ray(self, ray):
        """Uses the intercept and refract method to append the intercept and new direction to the ray as it intercepts a lens"""
        intercept = self.intercept(ray)
        if intercept is None:
            return "ray terminates: no intercept"
        else:
            newdir = self.refract(ray, intercept)
            if newdir is None:
                return "ray terminates: total internal reflection"
            else:
                ray.append(newdir, intercept)

class OutputPlane(OpticalElement):
    """
    Class representing the output plane. Initialise by defining a normal vector n and some point in the plane k, both as NumPy 3-arrays.
    """
    def __init__(self, normal, point):
        if normal.size != 3:
            raise Exception("Normal not 3-vector")
        elif point.size != 3:
            raise Exception("Point not 3-vector")
        self.normal = normal
        self.point = point
    
    def intercept(self, ray):
        """Calculate intercept of ray with output plane"""
        P = ray.pcurrent()
        m = self.point
        n = self.normal
        k = ray.dcurrent()
        Lambda = (np.dot(np.subtract(m, P), n)) / (np.dot(k, n))
        return P + Lambda*k
    
    def propagate_ray(self, ray):
        """Propagates ray to output plane. Appends point of intersection with plane to ray, and appends new direction as zero-vector"""
        intercept = self.intercept(ray)
        if intercept is None :
            return "ray terminates: no intercept with output plane"
        else:
            ray.append(np.array([0,0,0]), intercept)
            
class SphericalReflection(OpticalElement):
    """
    Class representing a spherical reflecting surface centred on the optical axis. z_0 is the intercept of the surface with the axis; curv, the curvature, is the reciprocal of the radius of curvature; n_1 and n_2 are the refractive indices inside and outside the surface; ap_rad is the aperture radius - the maximum extent of the surface from the optical axis.
    """
    def __init__(self, z_0, curv, n_1, n_2, ap_rad):
        self.z_0 = z_0
        self.curv = curv
        self.n_1 = n_1
        self.n_2 = n_2
        self.ap_rad = ap_rad
        
    def intercept(self,ray):
        """Calculates the first valid intercept of a ray with the spherical surface."""
        P = ray.pcurrent()
        k_hat=(ray.dcurrent())/(npl.norm(ray.dcurrent()))
        if self.curv == 0:
            if k_hat[-1] == 0 and P[-1] != self.z_0:
                return None
            elif P[-1] < self.z_0:
                if np.sign(k_hat[-1]) > 0:
                    Lambda = (self.z_0-P[2])/(k_hat[2])
                    intersect = np.add(P, Lambda*k_hat)
                    if np.sqrt((intersect[0])**2 + (intersect[1])**2) > self.ap_rad: 
                        return None
                    else:
                        return intersect
                else:
                    return None
            elif P[-1] > self.z_0:
                if np.sign(k_hat[-1]) < 0:
                    Lambda = (self.z_0-P[2])/(k_hat[2])
                    intersect = np.add(P, Lambda*k_hat)
                    if np.sqrt((intersect[0])**2 + (intersect[1])**2) > self.ap_rad: 
                        return None
                    else:
                        return intersect
                else:
                    return None
            else:
                return P
        else:
            R = 1.0/self.curv
            centrecurve = (np.array([0.0,0.0,self.z_0 + R])).astype(np.float)
            r = np.subtract(P, centrecurve)
            if ((np.dot(r, k_hat))**2 - ((npl.norm(r))**2 - R**2)) < 0:
                return None
            else:
                l_plus = -1 * (np.dot(r, k_hat)) + np.sqrt((np.dot(r, k_hat))**2 - ((npl.norm(r))**2 - R**2))
                l_minus = -1 * (np.dot(r, k_hat)) - np.sqrt((np.dot(r, k_hat))**2 - ((npl.norm(r))**2 - R**2))
                intersect_plus = np.add(P, l_plus * k_hat)
                intersect_minus = np.add(P, l_minus * k_hat)
                if np.sign(k_hat[-1]) == np.sign(R) and npl.norm(r) > R:
                    if np.sqrt((intersect_minus[0])**2 + (intersect_minus[1])**2) > self.ap_rad:
                        return None
                    else:
                        return intersect_minus
                elif np.sign(k_hat[-1]) != np.sign(R) and npl.norm(r) > R:
                    if np.sqrt((intersect_plus[0])**2 + (intersect_plus[1])**2) > self.ap_rad:
                        return None
                    else:
                        return intersect_plus
                elif np.sign(R) == np.sign(k_hat[-1]):
                    return None
                else:
                    if np.sqrt((intersect_minus[0])**2 + (intersect_minus[1])**2) > self.ap_rad:
                        return None
                    else:
                        return intersect_minus

    def reflect(self, ray, intercept):
        """Calculates new direction vector of ray due to reflection at surface"""
        intercept = intercept.astype(np.float)
        r_i = (np.subtract(intercept, ray.pcurrent())) / (npl.norm(np.subtract(intercept, ray.pcurrent())))
        if self.curv != 0:
            centrecurve = (np.array([0.0,0.0,self.z_0 + (1/self.curv)])).astype(np.float)
            normal_out = (np.subtract(intercept, centrecurve)) / (npl.norm(np.subtract(intercept, centrecurve)))
            normal_in = (np.subtract(centrecurve, intercept)) / (npl.norm(np.subtract(centrecurve, intercept)))
            if np.dot((-1*normal_out), r_i) > 0:
                normal = normal_out
            else:
                normal = normal_in
            r_r = r_i - 2*np.dot(r_i, normal)*normal
            return r_r
        else:
            normal_left = np.array([0,0,-1])
            normal_right = np.array([0,0,1])
            if np.dot((-1*normal_left), r_i) > 0:
                normal = normal_left
            else:
                normal = normal_right
            r_r = r_i - 2*np.dot(r_i, normal)*normal
            return r_r
        
    def propagate_ray(self, ray):
        """Uses the intercept and reflect method to append the intercept and new direction to the ray as it intercepts a reflecting object"""
        intercept = self.intercept(ray)
        if intercept is None:
            return "ray terminates: no intercept"
        else:
            newdir = self.reflect(ray, intercept)
            ray.append(newdir, intercept)
            
def ray_bundle_propagator(objects, n, rmax, m):
    """Propagates a bundle of rays starting at the x-y plane, evenly distributed over a disc of radius rmax, in n rings, with m points in the first ring, through a sequence of objects and output plane outplane. Also plots a ray diagram and calculates the rms spot radius at the output plane."""
    bundlegen = gp.rtuniform(n, rmax, m)
    bundlestartpts = [k for k in bundlegen]
    bundleendpts_x = []
    bundleendpts_y = []
    
    """Propagate ray bundle and plot ray diagram"""
    for j in range(len(bundlestartpts)):
        r = bundlestartpts[j][0]
        theta = bundlestartpts[j][1]
        ray = Ray(np.array([0,0,1]), np.array([r*np.cos(theta),r*np.sin(theta),0]))
        for object in objects:
            object.propagate_ray(ray)
        bundleendpts_x.append(ray.pcurrent()[0])
        bundleendpts_y.append(ray.pcurrent()[1])
        y_values_ray_plot = []
        z_values_ray_plot = []
        for i in range(len(ray.vertices())):
            point = ray.vertices()[i]
            y_values_ray_plot.append(point[1])
            z_values_ray_plot.append(point[2])
            pl.plot(z_values_ray_plot, y_values_ray_plot, 'blue')
    
    """Calculate RMS spot radius at output"""
    squareradii = []
    for i in range(len(bundleendpts_x)-1):
        squareradii.append(bundleendpts_x[i]**2 + bundleendpts_y[i]**2)
    rms_spot_radius = np.sqrt(sum(squareradii)/len(squareradii))
        
    return bundleendpts_x, bundleendpts_y, rms_spot_radius