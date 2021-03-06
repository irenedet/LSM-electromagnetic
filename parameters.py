from ngsolve import *
from netgen.csg import *
from ngsolve.internal import *

#Geometrical parameters
k = 3. # Wave number
geom = 4 # 1-->Brick with defective side, 2-->Sphere with circular layer, 3-->Brick with circular layer
Rext = 2.7 #Radius of external circle
Rpml = 2. #Radius of PML
alpha_pml = 0.6 # absorbing coefficient of PML
Rplus =  1.3 #Radius of Oplus
Rminus = 0.6 #Radius of Ominus
delta = 0.05 #Layer thickness

if (geom == 2):
    #Geometry of the circular delamination
    c = (2*delta+delta*delta)/(2+2*delta-sqrt(3)) #center of other
    Rother = 1-c+delta # Radius of other

#Physical parameters
mu2 = 1. #mu in oplus
mu1 = 1 #1.5 #mu in ominus
mu0 = 1 #2 #mu in the layer

eps2 = 2 #+ 1J*2.01 #epsilon in oplus
eps1 = 3.5 #+ 1J*3.51#2 + 1J*2.01 #epsilon in ominus
eps0 = 1 #+ 1J*1.01 #epsilon in the layer

#Nitsche's parameter
gamma=-1.e3*1J#check this

#FE parameters
order = 2 #order of the integration in the error computation
hmax=0.2
# Surface refinement to define points on the interface to be sampled
hsample = 0.05
# Far field pattern computation parameters
pi=4.*atan(1.)
FFh=0.4
data_dir='Data_'+str(k)+'_'+str(eps0)

# Visualization parameters 
plot_id = 1 # 1 = to plot the solutions

print('Meshing parameter',2*3.141/k/hmax*order)
