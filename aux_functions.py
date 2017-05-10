from ngsolve import *
from netgen.csg import *
from ngsolve.internal import *
import numpy as np
from time import time

#Cross product
Cross = lambda u,v: CoefficientFunction((u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])) # cross product
Ocross=lambda a,b: (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-b[0]*a[1]) # just gives a tuple

def materials(mu1,mu2,mu0,eps1,eps2,eps0,mesh):
    mur = {"ominus" : mu1 , "oplus" : mu2 , "olayer" : mu0 , "air" : 1, "pml" : 1 }
    nu = CoefficientFunction ([ 1/(mur[mat]) for mat in mesh.GetMaterials() ])

    epsilonr = {"ominus" : eps1 , "oplus" : eps2 , "olayer" : eps0 , "air" : 1, "pml" : 1 }
    epsilon = CoefficientFunction ([ epsilonr[mat] for mat in mesh.GetMaterials() ])
    return nu,epsilon

def Nitsches_transmission(alpha_pml,Rpml,nu,eps,k,mu2,nv,mesh,V,gamma,hmax,uext,uplus,vext,vplus):
    # u and v refer to trial and test-functions in the definition of forms below
    SetPMLParameters(rad=Rpml,alpha=alpha_pml)
    a = BilinearForm(V, symmetric=True, flags={"printelmat":True})
    a.components[0] += BFI("PML_curlcurledge",coef=1)
    a.components[0] += BFI("PML_massedge",coef=-k*k)
    a += SymbolicBFI(nu*curl(uplus)*curl(vplus)   - k*k*eps*uplus*vplus)


    # Nietsche's method for the transmission between the scattered field in the exterior and the total field inside the obstacle
    a += SymbolicBFI(0.5*( (1./mu2)*curl(uplus).Trace() + curl(uext).Trace()) * Cross(nv,vext.Trace()-vplus.Trace() ) ,BND,definedon=mesh.Boundaries("interface"))
    a += SymbolicBFI( Cross(nv,uext.Trace()-uplus.Trace()) * 0.5*( (1./mu2)*curl(vplus).Trace() + curl(vext).Trace() ) ,BND,definedon=mesh.Boundaries("interface"))
    a += SymbolicBFI( gamma/hmax*k* (uext.Trace() - uplus.Trace()) * (vext.Trace() - vplus.Trace()),BND,definedon=mesh.Boundaries("interface"))
    with TaskManager():
        a.Assemble()
        Minv = a.mat.Inverse(V.FreeDofs(),inverse='pardiso')
    return a, Minv

def NitschesPlaneWaveSource(k,dv,pv,nv,mesh,Vext,V,gamma,hmax,mu2,vext,vplus):
    #Incident field
    Ein= exp(k*1.J*(dv[0]*x+dv[1]*y+dv[2]*z))*CoefficientFunction(pv)
    curlEin = k*1.J*Cross(dv,pv)*exp(k*1.J*(dv[0]*x+dv[1]*y+dv[2]*z))
    Einext=GridFunction(Vext,"gEin")
    Einext.Set(Ein)

    f = LinearForm(V)
    f += SymbolicLFI( curlEin * 0.5 * (Cross(nv,vplus.Trace() + vext.Trace() )),BND,definedon=mesh.Boundaries("interface"))
    f += SymbolicLFI(Ein * 0.5*Cross(nv,(1./mu2)*curl(vplus).Trace()+curl(vext).Trace() ),BND,definedon=mesh.Boundaries("interface"))
    f += SymbolicLFI(-gamma/hmax*k* Ein * (vext.Trace()-vplus.Trace()),BND,definedon=mesh.Boundaries("interface"))
    with TaskManager():
        f.Assemble()
    
    return f,Einext

def Sample_Points(ngmesh,mesh,data_dir,Rminus,geom):
    #Sampling points for the LSM
    #SampleP=ngmesh.Points()
    SSpoints=[]# mip points where to evaluate the background sols
    SSpointsp=[]# simple list of sampling points
    if (geom == 1):
        spfile=open(data_dir+'/Spoints.txt','w')
        for p in ngmesh.Points():
            zsamp = p.p
            if(zsamp[0]*zsamp[0]==Rminus*Rminus)&(zsamp[1]*zsamp[1]<=Rminus*Rminus)&(zsamp[2]*zsamp[2]<=Rminus*Rminus):
                if (zsamp not in SSpointsp):
                    SSpointsp.append(zsamp)
                    mip=mesh(zsamp[0],zsamp[1],zsamp[2])
                    SSpoints.append(mip)
                    print(zsamp[0],zsamp[1],zsamp[2],file=spfile)
            if(zsamp[0]*zsamp[0]<=Rminus*Rminus)&(zsamp[1]*zsamp[1]==Rminus*Rminus)&(zsamp[2]*zsamp[2]<=Rminus*Rminus):
                if (zsamp not in SSpointsp):
                    SSpointsp.append(zsamp)
                    mip=mesh(zsamp[0],zsamp[1],zsamp[2])
                    SSpoints.append(mip)
                    print(zsamp[0],zsamp[1],zsamp[2],file=spfile)
            if(zsamp[0]*zsamp[0]<=Rminus*Rminus)&(zsamp[1]*zsamp[1]<=Rminus*Rminus)&(zsamp[2]*zsamp[2]==Rminus*Rminus):
                if (zsamp not in SSpointsp):
                    SSpointsp.append(zsamp)
                    mip=mesh(zsamp[0],zsamp[1],zsamp[2])
                    SSpoints.append(mip)
                    print(zsamp[0],zsamp[1],zsamp[2],file=spfile)
        Nsample = len(SSpoints)
        spfile.close()
    if (geom == 4):
        spfile=open(data_dir+'/Spoints.txt','w')
        for p in ngmesh.Points():
            zsamp = p.p
            if(zsamp[0]*zsamp[0]+zsamp[1]*zsamp[1]<=Rminus*Rminus)&(zsamp[2]==0):
                if (zsamp not in SSpointsp):
                    SSpointsp.append(zsamp)
                    mip=mesh(zsamp[0],zsamp[1],zsamp[2])
                    SSpoints.append(mip)
                    print(zsamp[0],zsamp[1],zsamp[2],file=spfile)
            if((zsamp[0]*zsamp[0]+zsamp[1]*zsamp[1]+zsamp[2]*zsamp[2])==Rminus*Rminus)&(zsamp[2]<0):
                if (zsamp not in SSpointsp):
                    SSpointsp.append(zsamp)
                    mip=mesh(zsamp[0],zsamp[1],zsamp[2])
                    SSpoints.append(mip)
                    print(zsamp[0],zsamp[1],zsamp[2],file=spfile)
        Nsample = len(SSpoints)
        spfile.close()
    
    if (geom == 3)|(geom == 2):
        print('The Sample_Points function has to be defined for this geometry')
    return Nsample,SSpoints,SSpointsp

def Get_polarizations1(dv):
    dvp1=(dv[2]-dv[1],dv[0]-dv[2],dv[1]-dv[0])
    pv0=Ocross(dv,dvp1)
    npv=sqrt(pv0[0]**2+pv0[1]**2+pv0[2]**2)
    if abs(npv)==0:
            print('Error: pv=0 ')
    pv0=(pv0[0]/npv ,pv0[1]/npv, pv0[2]/npv)
    # make sure there are no orthogonality errors
    pvdv=pv0[0]*dv[0]+pv0[1]*dv[1]+pv0[2]*dv[2]
    if abs(pvdv)>1.e-15:
        print('d is not orthogonal to p')
        xxstop
    pv1=Ocross(dv,pv0)
    npv=sqrt(pv1[0]**2+pv1[1]**2+pv1[2]**2)
    if abs(npv)==0:
        print('Error: pv=0 ')
        xxstop
    pv1=(pv1[0]/npv ,pv1[1]/npv, pv1[2]/npv)
    return pv0,pv1

def Get_polarizations(dv,np,FFpoints,FFP):
    dvp1=FFpoints[(np+1)%len(FFP)] # just to obtain a non-colinear
    # direction (this fails if FFpoints[(np+1)%len(FFP)]=-dv)
    pv0=Ocross(dv,dvp1)
    npv=sqrt(pv0[0]**2+pv0[1]**2+pv0[2]**2)
    if abs(npv)==0:
        dvp1=(dv[2]-dv[1],dv[0]-dv[2],dv[1]-dv[0])
        pv0=Ocross(dv,dvp1)
        npv=sqrt(pv0[0]**2+pv0[1]**2+pv0[2]**2)
        if abs(npv)==0:
            print('Error: pv=0 ')
            xxxstop
    pv0=(pv0[0]/npv ,pv0[1]/npv, pv0[2]/npv)
    # make sure there are no orthogonality errors
    pvdv=pv0[0]*dv[0]+pv0[1]*dv[1]+pv0[2]*dv[2]
    if abs(pvdv)>1.e-15:
        print('d is not orthogonal to p')
        xxstop
    pv1=Ocross(dv,pv0)
    npv=sqrt(pv1[0]**2+pv1[1]**2+pv1[2]**2)
    if abs(npv)==0:
        print('Error: pv=0 ')
        xxstop
    pv1=(pv1[0]/npv ,pv1[1]/npv, pv1[2]/npv)
    return pv0,pv1

def write_in_rhs(np,rhs,dv,pv,pol,SSpoints,u_RHS):    
    print(np, pol, file=rhs)
    print(dv[0], dv[1], dv[2], file=rhs)
    print(pv[0], pv[1], pv[2], file=rhs)
    for mip in SSpoints:
        j=u_RHS(mip)#evaluate at sampling points
        u1=j[0]
        u2=j[1]
        u3=j[2]
        print("%12.12e %12.12e %12.12e %12.12e %12.12e %12.12e"%(
            u1.real,u1.imag,u2.real,u2.imag,u3.real,u3.imag),file=rhs)

def write_in_fff(np,fff,dv,pol,pv,FF):
    print(np, pol, file=fff)
    print(dv[0], dv[1], dv[2], file=fff)
    print(pv[0], pv[1], pv[2], file=fff)
    for j in FF:
        F1=j[0]
        F2=j[1]
        F3=j[2]
        print("%12.12e %12.12e %12.12e %12.12e %12.12e %12.12e"%(
            F1.real,F1.imag,F2.real,F2.imag,F3.real,F3.imag),file=fff)
