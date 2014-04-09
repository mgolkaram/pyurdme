#!/usr/bin/env python
""" pyURDME model file for the annihilation cylinder 3D example. """

import os
import pyurdme
import dolfin

import matplotlib.pyplot as plt
import numpy

# Global Constants
MAX_X_DIM = 5.0
MIN_X_DIM = -5.0
TOL = 1e-1 


class Edge1(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] > (MAX_X_DIM - TOL)

class Edge2(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] < (MIN_X_DIM + TOL)

class cylinderDemo3D(pyurdme.URDMEModel):
    def __init__(self, model_name="cylinder_demo3d"):
        pyurdme.URDMEModel.__init__(self, model_name)

        # System constants
        D_const = 0.1
        k_react = pyurdme.Parameter(name="k_react", expression=1)
        k_creat = pyurdme.Parameter(name="k_creat", expression=200)
        self.addParameter([k_react, k_creat])
        
        # Define Species
        A = pyurdme.Species(name="A", diffusion_constant=D_const)
        B = pyurdme.Species(name="B", diffusion_constant=D_const)
        self.addSpecies([A, B])
        
        # Define Geometry
        pt1 = dolfin.Point(MAX_X_DIM, 0, 0)
        pt2 = dolfin.Point(MIN_X_DIM, 0, 0)
        cylinder = dolfin.Cylinder(pt1, pt2, 1.0)
        self.mesh = pyurdme.URDMEMesh(mesh=dolfin.Mesh(cylinder, 32))
        
        # Define Subdomains
        subdomains = dolfin.MeshFunction("size_t", self.mesh, self.mesh.topology().dim()-1)
        subdomains.set_all(1)
        # Mark the boundary points
        Edge1().mark(subdomains,2)
        Edge2().mark(subdomains,3)
        self.addSubDomain(subdomains)
        
        # Define Reactions
        R1 = pyurdme.Reaction(name="R1", reactants=None, products={A:1}, rate=k_creat, restrict_to=2)
        R2 = pyurdme.Reaction(name="R2", reactants=None, products={B:1}, rate=k_creat, restrict_to=3)
        R3 = pyurdme.Reaction(name="R3", reactants={A:1, B:1}, products=None, rate=k_react)
        self.addReaction([R1, R2, R3])

        # Define simulation timespan
        self.timespan(range(200))


if __name__ == "__main__":
    
    model = cylinderDemo3D()
    result = pyurdme.urdme(model, report_level=1)
    
    # This line here dumps the state of A at all timepoints to Paraview comaptible output (VTK). The trajectory
    # is written to a folder "Aout", where each snapshot is stored in a separate file. To open the "movie",
    # just open Aout/trajectory.pvd, then you can animate etc.
    if not os.path.isdir('Aout'):
        print "Writing species 'A' to folder 'Aout'"
        result.toVTK(species='A',folder_name="Aout")
    if not os.path.isdir('Bout'):
        print "Writing species 'B' to folder 'Bout'"
        result.toVTK(species='B',folder_name="Bout")


    # Plot of the time-average spatial concentration.
    x_vals = model.mesh.coordinates()[:, 0]
    A_vals = numpy.mean(result.getSpecies("A", concentration=True), axis=0)
    B_vals = numpy.mean(result.getSpecies("B", concentration=True), axis=0)
    plt.plot(x_vals,A_vals,'.r',x_vals,B_vals,'.b')
    plt.legend(['A', 'B'])
    plt.show()
