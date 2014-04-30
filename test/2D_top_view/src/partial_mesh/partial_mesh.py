###
###     2D COMPLETE MESH (obstacle + object)
###

# USER DEFINED PARAMETERS
mesh_max_size = 0.2
refined_mesh_max_size = 0.1


import sys
import os
import salome
import math

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)
sys.path.insert( 0, os.getcwd() )

###
### GEOMETRY
###

import GEOM
from salome.geom import geomBuilder
import SALOMEDS

# create new geometry object
geompy = geomBuilder.New(theStudy)

# create vertices
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

# vertex are number counterclock-wise, starting from the bottom left
domain_vertex_1 = geompy.MakeVertex (0, 0, 0);
domain_vertex_2 = geompy.MakeVertex (10, 0, 0);
domain_vertex_3 = geompy.MakeVertex (10, 7, 0);
domain_vertex_4 = geompy.MakeVertex (0, 7, 0);
circle_center   = geompy.MakeVertex (3.5, 3.5, 0);
circle_radius   = 0.5;


# create lines
domain_line_1 = geompy.MakeLineTwoPnt (domain_vertex_1, domain_vertex_2);
domain_line_2 = geompy.MakeLineTwoPnt (domain_vertex_2, domain_vertex_3);
domain_line_3 = geompy.MakeLineTwoPnt (domain_vertex_3, domain_vertex_4);
domain_line_4 = geompy.MakeLineTwoPnt (domain_vertex_4, domain_vertex_1);
circumference = geompy.MakeCircle (circle_center, None, circle_radius)

# create faces
domain_face = geompy.MakeFaceWires ([domain_line_1, domain_line_2, domain_line_3, domain_line_4, circumference], 1)

# create groups
refinement_group = geompy.CreateGroup (domain_face, geompy.ShapeType["EDGE"])
geompy.UnionList (refinement_group, [circumference])

## add everything to what gui will show
geompy.addToStudy (O, 'O')
geompy.addToStudy (OX, 'OX')
geompy.addToStudy (OY, 'OY')
geompy.addToStudy (OZ, 'OZ')
geompy.addToStudy (domain_vertex_1, 'domain_vertex_1')
geompy.addToStudy (domain_vertex_2, 'domain_vertex_2')
geompy.addToStudy (domain_vertex_3, 'domain_vertex_3')
geompy.addToStudy (domain_vertex_4, 'domain_vertex_4')
geompy.addToStudy (circle_center, 'circle_center')
geompy.addToStudy (domain_line_1, 'domain_line_1')
geompy.addToStudy (domain_line_2, 'domain_line_2')
geompy.addToStudy (domain_line_3, 'domain_line_3')
geompy.addToStudy (domain_line_4, 'domain_line_4')
geompy.addToStudy (circumference, 'circumference')
geompy.addToStudy (domain_face, 'domain_face')
geompy.addToStudyInFather (domain_face, refinement_group, 'refinement_group')

####
#### MESH
####

#import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(theStudy)
from salome.StdMeshers import StdMeshersBuilder
from salome.NETGENPlugin import NETGENPluginBuilder

# create new mesh object
mesh = smesh.Mesh(domain_face)

# global mesh parameters
NETGEN_2D = mesh.Triangle (algo=smeshBuilder.NETGEN_1D2D)
mesh_parameters = NETGEN_2D.Parameters()
mesh_parameters.SetMaxSize (mesh_max_size)
mesh_parameters.SetSecondOrder (0)
mesh_parameters.SetOptimize (1)
mesh_parameters.SetFineness (3)
mesh_parameters.SetMinSize (0)
mesh_parameters.SetQuadAllowed (0)

# refine mesh
submesh_algorithm = mesh.Segment (geom=refinement_group)
submesh_algorithm.MaxSize (refined_mesh_max_size)

# compute mesh
isDone = mesh.Compute()

### set object names
smesh.SetName(mesh.GetMesh(), 'mesh')
smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN_2D')
smesh.SetName(mesh_parameters, 'mesh_parameters')

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(1)

# export mesh in unv format
mesh.ExportUNV (os.getcwd () + '/partial_mesh.unv')

