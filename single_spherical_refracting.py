# -*- coding: utf-8 -*-
"""
Created on Fri Dec 02 14:34:38 2016

@author: rc1515
"""

"""
Script that propagates a bundle of rays through a single spherical refracting surface, as in "getting started" section. Saves three figures: a ray diagram; a spot diagram of the ray source; a spot diagram at the output plane.
"""
import numpy as np
import raytracer as rt
import matplotlib.pyplot as pl
import genpolar as gp

reload(rt)
reload(gp)

lens = rt.SphericalRefraction(100, 0.03, 1.0, 1.5168, 1.0/0.03)
outplane = rt.OutputPlane(np.array([0,0,-1]), np.array([0,0,197.833]))

bundleendpts_x, bundleendpts_y, rms_spot_radius = rt.ray_bundle_propagator([lens, outplane], 8, 5, 6)

"""Plot Ray Diagram"""
pl.grid()
pl.xlabel("z/mm")
pl.ylabel("y/mm")
pl.axis([0, 200, -5.2, 5.2])
pl.title('Ray Diagram')
pl.savefig('ray_diagram_single_spherical_refract.png', dpi=300, bbox_inches='tight')
pl.clf()

"""Plot input spot diagram"""
for r, t in gp.rtuniform(8, 2.5, 6):
    pl.plot(r * np.cos(t), r * np.sin(t), 'bo')
pl.title('Input')
pl.xlabel("x/mm")
pl.ylabel("y/mm")
pl.axis('square')
pl.grid()
pl.savefig('input_single_spherical_refract.png', dpi=300, bbox_inches='tight')
pl.clf()

"""Plot output spot diagram"""
pl.plot(bundleendpts_x, bundleendpts_y, 'bo')
pl.title('Output')
pl.xlabel("x/mm")
pl.ylabel("y/mm")
pl.axis('square')
pl.grid()
#pl.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#pl.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
pl.savefig('output_single_spherical_refract.png', dpi=300, bbox_inches='tight')
pl.clf()