from ngsolve import *
from netgen.csg import *
from ngsolve.internal import *


def brick2_geometry(Rminus, Rplus, Rext, Rpml, delta, hsample, hmax):
    geometry = CSGeometry()
    o_minus = (OrthoBrick(Pnt(-Rminus,-Rminus,-Rminus),Pnt(Rminus,Rminus,Rminus)).maxh(0.1))
    
    geometry.Add (o_minus.mat("ominus").maxh(0.3))    
    #geometry.Add ((box * pl1 * pl2).mat("olayer").maxh(hmax),bcmod=[(pl1,"crack"),(box,"sides"),(pl2,"top")])

    #slices = [2**(-i) for i in reversed(range(1,6))]
    #geometry.CloseSurfaces(pl1,pl2)#,slices)
    
    return geometry
def brick_geometry(Rminus, Rplus, Rext, Rpml, delta, hsample, hmax):
    geometry = CSGeometry()
    o_ext = (Sphere(Pnt(0,0,0), Rext)).bc("outer")
    pml = Sphere(Pnt(0,0,0),Rpml)
    o_plus = Sphere(Pnt(0,0,0), Rplus).bc("interface")
    #This is to define the two close surfaces for the thin layer:
    box = OrthoBrick(Pnt(-Rminus,-Rminus,-Rminus),Pnt(Rminus,Rminus,Rminus+delta))
    pl1 = Plane(Pnt(0,0,Rminus),Vec(0,0,-1)).bc("crack")
    pl2 = Plane(Pnt(0,0,Rminus+delta),Vec(0,0,1))
    o_minus = (box - pl1).maxh(hsample)
    
    geometry.Add (o_minus.mat("ominus").maxh(hmax))    
    geometry.Add ((o_ext - pml).mat("pml"))
    geometry.Add ((pml-o_plus).mat("air"))
    geometry.Add ((o_plus-box).mat("oplus").maxh(hmax))
    geometry.Add ((box * pl1 * pl2).mat("olayer").maxh(hmax))

    #slices = [2**(-i) for i in reversed(range(1,6))]
    geometry.CloseSurfaces(pl1,pl2)#,slices)
    
    return geometry

def irreg_geometry(Rminus, Rplus, Rext, Rpml, delta, hsample, hmax):
    geometry = CSGeometry()
    o_ext = (Sphere(Pnt(0,0,0), Rext)).bc("outer")
    pml = Sphere(Pnt(0,0,0),Rpml)
    o_plus = Sphere(Pnt(0,0,0), Rplus + 0.1 ).bc("interface")

    box = OrthoBrick(Pnt(-Rminus,-Rminus,-Rminus),Pnt(Rminus,Rminus,Rminus+delta)).maxh(hsample)
    pl1 = Plane(Pnt(0,0,Rminus),Vec(0,0,-1)).bc("crack")
    pl2 = Plane(Pnt(0,0,Rminus+delta),Vec(0,0,1))
    small_cylinder = Cylinder(Pnt(0.4,0.4,1),Pnt(0.4,0.4,-1),0.4)
    big_cylinder = Cylinder(Pnt(1,1,1),Pnt(-1,-1,-1),0.6).maxh(hsample)
    layer = (box * pl1 * pl2).maxh(hsample)
        
    box_minus = ((box-pl1) * big_cylinder).maxh(hsample).bc("sample")
    box_plus = ((box-pl1) - big_cylinder).maxh(hsample)
    layer = (box * pl1 * pl2)

    geometry.Add ((o_ext - pml).mat("pml"))
    geometry.Add ((pml-o_plus).mat("air"))
    geometry.Add ((o_plus-box).mat("oplus").maxh(hmax))
    geometry.Add ((layer * small_cylinder).mat("olayer").maxh(hsample))
    geometry.Add ((layer - small_cylinder).mat("oplus").maxh(hsample))
    geometry.Add (box_plus.mat("oplus").maxh(hsample))
    geometry.Add (box_minus.mat("ominus").maxh(hsample))
    
    return geometry

def half_sphere(Rminus, Rplus, Rext, Rpml, delta, hsample, hmax):
    geometry = CSGeometry()        

    o_ext = (Sphere(Pnt(0,0,0), Rext)).bc("outer")
    pml = Sphere(Pnt(0,0,0),Rpml)
    o_plus = Sphere(Pnt(0,0,0), Rplus+.2).bc("interface")

    #This is to define the two close surfaces for the thin layer:
    circle = Sphere(Pnt(0,0,0),Rminus)
    small_cylinder = Cylinder(Pnt(0,0,1),Pnt(0,0,-1),0.3)

    pl1 = Plane(Pnt(0,0,0),Vec(0,0,-1))
    pl2 = Plane(Pnt(0,0,delta),Vec(0,0,1))
    o_minus_with_layer = (circle*pl2)
    o_minus = (circle - pl1).maxh(hsample).bc("sample")
    layer = (circle * pl1 * pl2)

    geometry.Add ((layer * small_cylinder).mat("olayer").maxh(hsample))
    geometry.Add ((layer - small_cylinder).mat("oplus").maxh(hsample))

    slices = [2**(-i) for i in reversed(range(1,6))]
    geometry.CloseSurfaces(pl1,pl2,slices)
    
    geometry.Add ((o_ext - pml).mat("pml"))
    geometry.Add ((pml-o_plus).mat("air"))
    geometry.Add ((o_plus-o_minus_with_layer).mat("oplus").maxh(hmax))
    geometry.Add ((o_minus).mat("ominus").maxh(hsample))
    
    return geometry

def sphere_geometry(Rminus, Rplus, Rext, Rpml, Rother, c, delta, hsample, hmax):
    geometry = CSGeometry()
    o_ext = (Sphere(Pnt(0,0,0), Rext)).bc("outer")
    pml = Sphere(Pnt(0,0,0),Rpml)
    o_plus = Sphere(Pnt(0,0,0), Rplus).bc("interface")

    #This is to define the two close surfaces for the thin layer:
    
    o_minus = (Sphere(Pnt(0,0,0), Rminus)).maxh(hmax).mat("ominus").maxh(hsample)
    other = Sphere(Pnt(0,c,0),Rother)
    withCrack = (o_minus * other)
    withoutCrack = (o_minus - other)

    geometry.Add ((o_ext - pml).mat("air"))
    geometry.Add ((pml-o_plus).mat("air"))
    geometry.Add ((o_plus-o_minus-other).mat("oplus"))
    geometry.Add ((other-o_minus).mat("olayer"))
    geometry.Add (withCrack,bcmod=[(o_minus,"crack")])
    geometry.Add (withoutCrack,bcmod=[(o_minus,"nocrack")])

    return geometry
