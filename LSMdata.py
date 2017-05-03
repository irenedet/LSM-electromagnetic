from netgen.csg import *
from ngsolve import *
from ngsolve.internal import *
from netgen.meshing import SetTestoutFile
from ffcode import ExportFFMesh, SetupFFMesh, findffvec
from mygeometry_1 import brick_geometry, sphere_geometry
from aux_functions import *
import os

## Load Physical, Geometrical, and Computational parameters
from parameters import *

################################################################################################################
#      Definition of the geometry
################################################################################################################
if geom == 1:
    geometry = brick_geometry(Rminus, Rplus, Rext, Rpml, delta, hsample, hmax)
    ngmesh = geometry.GenerateMesh()

if geom == 2:
    geometry = sphere_geometry(Rminus, Rplus, Rext, Rpml, Rother, c, delta, hsample, hmax)
    ngmesh = MakeGeometry.GenerateMesh(maxh=hmax)
    
    # This is for the boundary layer
    # 4 = number of boundary conditions:
    for i in range(6):
        print("BC: ",ngmesh.GetBCName(i))
        if ngmesh.GetBCName(i) == "crack":
            crackIndex = i+1

    #the crack should go into domain o_plus, there are 3 materials defined
    for i in range(1,6):
        print("Material: ",ngmesh.GetMaterial(i))
        if ngmesh.GetMaterial(i) == "oplus":
            crackDomIndex = i+1
    ngmesh = MakeGeometry().GenerateMesh(maxh=hmax)
        
    # create boundary layer with: surface index, thickness, index of domain into which the layer will go, new material name of layer
    ngmesh.BoundaryLayer(crackIndex,delta,crackDomIndex,"olayer")

mesh = Mesh(ngmesh)
print(mesh.GetBoundaries())

if (plot_id ==1):
    Draw(mesh)
    viewoptions.clipping.enable = 1
    viewoptions.clipping.nx = -1
    visoptions.clipsolution = "scal"

nv=specialcf.normal(mesh.dim)
collection='interface'

# Finite element spaces involved in the variational formulations:
Vext = HCurl(mesh, order=2, complex=True, definedon=mesh.Materials("air")+mesh.Materials("pml"), dirichlet="outer")
Vplus  = HCurl(mesh, order=2, complex=True, definedon=mesh.Materials("oplus")+mesh.Materials("olayer")+mesh.Materials("ominus"))
V=FESpace([Vext,Vplus])

# Define test and trial functions
uext,uplus = V.TrialFunction()
vext,vplus = V.TestFunction()

################################################################################################################
#     Sesquilinear forms & assemblement of the FEM matrix using a Nitsche's method on the 
#     exterior boundary of the inhomogeneity ("interface").
################################################################################################################
# Material properties:
nu_b, eps_b = materials(mu1,mu2,mu2,eps1,eps2,eps2,mesh) #Background media material functions
nu_l, eps_l = materials(mu1,mu2,mu0,eps1,eps2,eps0,mesh) #Defective media material functions

# FEM Matrices
a_b = Nitsches_transmission(alpha_pml,Rpml,nu_b,eps_b,k,mu2,nv,mesh,V,gamma,hmax,uext,uplus,vext,vplus) #Background domain
a_l = Nitsches_transmission(alpha_pml,Rpml,nu_l,eps_l,k,mu2,nv,mesh,V,gamma,hmax,uext,uplus,vext,vplus) #Defective domain

#################################################################################################################
#    Data for the LSM (Linear Sampling Method) : Far-field data and right-hand-side computation
#################################################################################################################
print('Saving to '+data_dir)
if not os.path.exists(data_dir):
    os.makedirs(data_dir)


## Mesh of unit sphere gives directions for incoming waves/measurements for FF pattern
ffgeometry = CSGeometry()
ffgeometry.Add(Sphere(Pnt(0,0,0),1).bc('ff'))
ffnetgenMesh = ffgeometry.GenerateMesh(maxh=FFh)

FFP,bfaces,FFpoints = SetupFFMesh(ffnetgenMesh)
ExportFFMesh(FFP,bfaces,data_dir+'/ffmesh.msh')
ffnetgenMesh.Save(data_dir+'/ffmesh.vol') # save for imepdance calculation
N=len(FFP)
## Set sampling points for the LSM
Nsample,SSpoints,SSpointsp = Sample_Points(ngmesh,mesh,data_dir,Rminus) 

#Cross product
Cross = lambda u,v: CoefficientFunction((u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])) # cross product
Ocross=lambda a,b: (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-b[0]*a[1]) # just gives a tuple
nopml = CoefficientFunction ([0 if mat=='pml' else 1 for mat in mesh.GetMaterials() ]) #indicator function for the physical domain
##################################################################################################################
#     Loop for FF and rhs computation, when sending all the different incident plane waves
##################################################################################################################

np=0
#fff=open(data_dir+'/FFP.txt','w')
#print(k,file=fff)
#print(eps0,file=fff)
#print(delta,file=fff)
#collection='interface'
rhs=open(data_dir+'/RHS.txt','w')
print(Nsample,file=rhs)

for d_direc in FFP:
    dv=d_direc.p
    dvm = (-dv[0],-dv[1],-dv[2]) #necessary to later compute RHS of the FF equation
    ndv=sqrt(dv[0]**2+dv[1]**2+dv[2]**2)
    if ndv>1.-0.00001:    
        pv0,pv1 = Get_polarizations(dv,N)

        #Definition of GridFunctions for ff computation
        u_b0 = GridFunction(V)
        u_l0 = GridFunction(V) 
        Es0 = GridFunction(Vext,"Es0")

        # Define solution functions
        u_b1 = GridFunction(V)# Background problem solution
        u_l1 = GridFunction(V)
        Es1=GridFunction(Vext,"Es0")

        #Definition of GridFunctions for rhs evaluation
        u_rhs0 = GridFunction(V)
        u_rhs1 = GridFunction(V)

        # Incident field for the first polarization pv0
        f0, Einext0 = NitschesPlaneWaveSource(k,dv,pv0,nv,mesh,Vext,V,gamma,hmax,mu2,vext,vplus)
        f1, Einext1 = NitschesPlaneWaveSource(k,dv,pv1,nv,mesh,Vext,V,gamma,hmax,mu2,vext,vplus)
        fm0, Einextm0 = NitschesPlaneWaveSource(k,dvm,pv0,nv,mesh,Vext,V,gamma,hmax,mu2,vext,vplus)
        fm1, Einextm1 = NitschesPlaneWaveSource(k,dvm,pv1,nv,mesh,Vext,V,gamma,hmax,mu2,vext,vplus)
        # Define solution functions
        with TaskManager():
            u_b0.vec.data += a_b.mat.Inverse(V.FreeDofs()) * f0.vec
            u_l0.vec.data += a_l.mat.Inverse(V.FreeDofs()) * f0.vec
            u_rhs0.vec.data += a_b.mat.Inverse(V.FreeDofs()) * fm0.vec
            u_b1.vec.data += a_b.mat.Inverse(V.FreeDofs()) * f1.vec
            u_l1.vec.data += a_l.mat.Inverse(V.FreeDofs()) * f1.vec
            u_rhs1.vec.data += a_b.mat.Inverse(V.FreeDofs()) * fm1.vec
        # Scattered field for the defect alone - To construct the FF operator in the LSM method
        Es0.Set(u_l0.components[0]-u_b0.components[0])     
        Es1.Set(u_l1.components[0]-u_b1.components[0])
        u_RHS0 = CoefficientFunction(u_rhs0.components[1])#to be able to evaluate at sample points in SSpoints
        u_RHS1 = CoefficientFunction(u_rhs1.components[1])#to be able to evaluate at sample points in SSpoints

        write_in_rhs(np,rhs,dv,pv0,0,SSpoints,u_RHS0)
        write_in_rhs(np,rhs,dv,pv1,1,SSpoints,u_RHS1)

        if (plot_id ==1):
            Draw(CoefficientFunction(u_rhs1.components[1]) +nopml*CoefficientFunction(u_rhs1.components[0])+nopml*Einextm1,mesh,"E_rhs1")
            Redraw()
################## FFcomputation #####################
        #print('starting ffvec')
        #FF0=findffvec(k,FFP,Es0,order,mesh,nv,collection,Vext)
        #FF1=findffvec(k,FFP,Es1,order,mesh,nv,collection,Vext)
        #print("done")
        ##xxx  # stop
        #write_in_fff(np,fff,dv,0,pv0,FF0)
        #write_in_fff(np,fff,dv,1,pv1,FF1)
        #np=np+1
        #print('MS: Done ',np,' incident fields out of ',len(FFP))
rhs.close()
#fff.close()

if (plot_id ==1):
    viewoptions.clipping.enable = 1
    viewoptions.clipping.nx = -1
    visoptions.clipsolution = "scal"



