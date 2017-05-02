from ngsolve import *
from netgen.csg import *
from ngsolve.internal import *
import numpy as np
from time import time

pi=4.*atan(1.)

Ocross=lambda a,b: (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],
                    a[0]*b[1]-b[0]*a[1]) # just gives a tuple
cross=lambda a,b: CoefficientFunction((a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],
                                       a[0]*b[1]-b[0]*a[1]))  # gives a

def SetupFFMesh(ffnetgenMesh):
    FFP=ffnetgenMesh.Points()
    rad=[]
    FFpoints=[]
    for p in FFP:
        xhat=p.p
        nrm=sqrt(xhat[1]**2+xhat[2]**2+xhat[0]**2)
        #print(np,' rad ',nrm)
        rad.append(nrm)
        FFpoints.append(xhat)
    bfaces=[]
    nf=0
    delt=0.0001
    for faces in ffnetgenMesh.Elements2D():
        iv=faces.vertices # Use -1 in next line to allow for
                          #c++/python index mismatch
        if (rad[iv[0].nr-1]>1-delt) & (rad[iv[1].nr-1]>
                                       1-delt) & (rad[iv[2].nr-1]>1-delt):
            bfaces.append([iv[0].nr,iv[1].nr,iv[2].nr])
            nf=nf+1
    print('Found ',nf,' boundary faces')
    return (FFP,bfaces,FFpoints)


def ExportFFMesh(FFP,bfaces, filename):
    """ export FFmesh """
    print ("export FFmesh in neutral format to file = ", filename)
    f = open (filename, 'w')
    print (len(FFP), file=f)
    for p in FFP:
        print (p.p[0], p.p[1], p.p[2], file=f)
    surfels = bfaces;
    print (len(surfels), file=f)
    for el in surfels:
        print (el[0],el[1],el[2], file=f) # no need for +1 since this is a netgen array

def findff6B(ka,FFP,Es,order,mesh,nv,dv,pv,collection,fespace):
    # This version computes the bistatic ff pattern and returns
    # it to the calling program
    #print('Using findff6B....')
    w=fespace.TestFunction()
    FF=[]
    np=0
    vv=GridFunction(fespace)
    zer0=CoefficientFunction((0.,0.,0.))
    vv.Set(zer0,VOL)
    for p in FFP:
        xhat=p.p
        nxhat=sqrt(xhat[0]**2+xhat[1]**2+xhat[2]**2)
        if nxhat>1-.0001:  # the ff mesh contains some interior points
            ff=[0,0,0]
            CEout = exp(1J*ka*(x*xhat[0]+y*xhat[1]+z*xhat[2]))
            # Conjugate to allow for conjugation in linear forms
            for j in range(0,3):
                e=[0.,0.,0.]
                e[j]=1.
                ## Functional for the full FF pattern
                func=cross(Ocross(e,xhat),nv)*CEout*(-1J*ka)
                gf = LinearForm(fespace)
                gf += SymbolicLFI(func*w.Trace(),BND,definedon=
                                  mesh.Boundaries(collection))
                ## Functional for the second part using finite element cutoff
                eT=cross(xhat,Ocross(e,xhat))
                vv.Set(eT*CEout,BND,definedon=mesh.Boundaries(collection))
                gf += SymbolicLFI(curl(vv)*curl(w)-ka*ka*vv*w,VOL,
                                  definedon=mesh.Materials('air'))
                with TaskManager():
                    gf.Assemble()
                # Conjugates are needed since inner product conjugates
                # the second component(?)
                ff1a=gf.vec.InnerProduct(Es.vec)
                ff[j]=ff1a.conjugate()/4./pi
            FF.append(ff)
        else:
            FF.append([-9999,-9999,-9999])
        #print('Done ',np,' of ',len(FFP))
        np=np+1
    return(FF)

def findffnew(ka,FFP,Es,order,mesh,nv,collection,fespace):
    # This version computes the bistatic ff pattern and returns
    # it to the calling program
    print('Using findffnew....')
    start = time()
    FF=[]
    vv=GridFunction(fespace)
    zer0=CoefficientFunction((0.,0.,0.))
    nEs=cross(nv,Es)
    points = []
    for p in FFP:
        xhat = p.p
        nxhat=sqrt(xhat[0]**2+xhat[1]**2+xhat[2]**2)
        if nxhat>1-.0001:  # the ff mesh contains some interior points
            points.append(p.p)
        
    for xhat in points:
        print("p: ",xhat)
        ff0=[0,0,0]
        ff1=[0,0,0]
        Eout = exp(-1J*ka*(x*xhat[0]+y*xhat[1]+z*xhat[2]))
        for j in range(0,3):
            func=Eout * nEs[j]
            with TaskManager():
                ff0[j] = (1j*ka/4./pi)*Integrate(func,mesh,BND,order=order+1,
                                                 definedon=mesh.Boundaries(collection))
        ff0=Ocross(xhat,ff0)
        for j in range(0,3):
            e=[0.,0.,0.]
            e[j]=1.
            ## Second part using finite element cutoff
            eT=cross(xhat,Ocross(e,xhat))
            vv.vec[:] = 0
            #vv.Set(zer0,VOL)
            vv.Set(eT*Eout,BND,definedon=mesh.Boundaries(collection))
            func=CoefficientFunction((curl(vv)*curl(Es)-ka*ka*vv*Es))
            ### Does the above conjugate....?
            with TaskManager():
                ff1[j] = (1./4./pi)*Integrate(func,mesh,VOL,order=order+1,
                                                 definedon=mesh.Materials('air'))
        ff=[ff0[j]+ff1[j] for j in range(3)]
        FF.append(ff)
    #print('Done ',np,' of ',len(FFP))
    print("ffnew needed ", time()-start, " seconds")
    return(FF)

def findffvec(ka,FFP,Es,order,mesh,nv,collection,fespace):
    # This version computes the bistatic ff pattern and returns
    # it to the calling program --- vectorized version
    #print('Using findffnewvec....')
    start = time()
    FF=[]
    funcsBND = []
    funcsVOL = []
    zer0=CoefficientFunction((0.,0.,0.))
    nEs=cross(nv,Es)
    index = 0
    points = []
    for p in FFP:
        xhat = p.p
        nxhat=sqrt(xhat[0]**2+xhat[1]**2+xhat[2]**2)
        if nxhat>1-.0001:  # the ff mesh contains some interior points
            points.append(p.p)
        
    for xhat in points:
        #print("p: ",xhat)
        astart = time()
        ff0=[0,0,0]
        ff1=[0,0,0]
        Eout = exp(-1J*ka*(x*xhat[0]+y*xhat[1]+z*xhat[2]))
        for j in range(0,3):
            func=Eout * nEs[j]
            funcsBND.append(func)
        astart = time()
        for j in range(0,3):
            e=[0.,0.,0.]
            e[j]=1.
            ## Second part using finite element cutoff
            eT=cross(xhat,Ocross(e,xhat))
            #vv.Set(zer0,VOL)
            vv=GridFunction(fespace)
            vv.vec[:] = 0
            vv.Set(eT*Eout,BND,definedon=mesh.Boundaries(collection))
            func=CoefficientFunction(curl(vv)*curl(Es)-ka*ka*vv*Es)
            ### Does the above conjugate....?
            funcsVOL.append(func)
    with TaskManager():
        ff0 = Integrate(CoefficientFunction(tuple(funcsBND)),mesh,order=order+1,definedon=mesh.Boundaries(collection), heapsize=10000000)
        ff1 = Integrate(CoefficientFunction(tuple(funcsVOL)),mesh,order=order+1,definedon=mesh.Materials("air"),heapsize=10000000)
    for i in range(len(points)):
        ff0list = []
        iVOL = []
        for j in range(3):
            ff0list.append((1j*ka/4./pi)*ff0[i*3+j])
            iVOL.append((1./4./pi) * ff1[i*3+j])
        iBND = Ocross(points[i],ff0list)
        ff = [iBND[j] + iVOL[j] for j in range(3)]
        FF.append(ff)
    print("ffvec needed ", time()-start, " seconds")
    return(FF)
    

def findff6B2(ka,FFP,Es0,Es1,order,mesh,nv,collection,fespace):
    # This version computes the bistatic ff pattern and returns
    # it to the calling program, does both polarizations at once
    # print('Using findff6B....')
    start = time()
    w=fespace.TestFunction()
    FF0=[]
    FF1=[]
    np=0
    vv=GridFunction(fespace)
    zer0=CoefficientFunction((0.,0.,0.))
    vv.Set(zer0,VOL)
    points = []
    for p in FFP:
        xhat = p.p
        nxhat=sqrt(xhat[0]**2+xhat[1]**2+xhat[2]**2)
        if nxhat>1-.0001:  # the ff mesh contains some interior points
            points.append(p.p)
        
    for xhat in points:
        print("p: ",xhat)
        ff0=[0,0,0]
        ff1=[0,0,0]
        CEout = exp(1J*ka*(x*xhat[0]+y*xhat[1]+z*xhat[2]))
        # Conjugate to allow for conjugation in linear forms
        for j in range(0,3):
            e=[0.,0.,0.]
            e[j]=1.
            ## Functional for the full FF pattern
            func=cross(Ocross(e,xhat),nv)*CEout*(-1J*ka)
            gf = LinearForm(fespace)
            #gf.vec[:]=0
            gf += SymbolicLFI(func*w.Trace(),BND,definedon=
                              mesh.Boundaries(collection))
            ## Functional for the second part using finite element cutoff
            eT=cross(xhat,Ocross(e,xhat))
            vv.vec[:] = 0
            vv.Set(eT*CEout,BND,definedon=mesh.Boundaries(collection))
            gf += SymbolicLFI(curl(vv)*curl(w)-ka*ka*vv*w,VOL,
                              definedon=mesh.Materials('air'))
            with TaskManager():
                gf.Assemble()
            # Conjugates are needed since inner product conjugates
            # the second component(?)
            ff1a0=gf.vec.InnerProduct(Es0.vec)
            ff0[j]=ff1a0.conjugate()/4./pi
            ff1a1=gf.vec.InnerProduct(Es1.vec)
            ff1[j]=ff1a1.conjugate()/4./pi
        FF0.append(ff0)
        FF1.append(ff1)
    #print('FF Done ',np,' of ',len(FFP)-1)
    print("ff6b2 needed ", time()-start, " seconds")
    return(FF0,FF1)


def findff6(ka,Es,order,mesh,nv,dv,pv,collection,nang,fespace,idir,fff):
    # Use a finite element cutoff function and linear forms
    # This version outputs the ff pattern as a function of azimuthal angle
    # for comparison to MatScat (two calls are needed)
    print('Using findff6....')
    w=fespace.TestFunction()
    FF=[]
    ang=[]
    zer0=CoefficientFunction((0,0,0))
    vv=GridFunction(fespace)
    vv.Set(zer0,VOL)
    if idir==0:
       print('Computing xz plane')
    else:
       print('Computing yz plane')
    for na in range(0,nang):
        ff=[0,0,0]
        theta=pi*na/(nang-1)
        if idir==0:
           xhat=(sin(theta),0,cos(theta))
        else:
           xhat=(0,sin(theta),cos(theta))
        CEout = exp(1J*ka*(x*xhat[0]+y*xhat[1]+z*xhat[2]))
        # Conjugate to to allow for conjugation in linear forms
        for j in range(0,3):
            e=[0.,0.,0.]
            e[j]=1.
            ## Functional for the full FF pattern
            func=cross(Ocross(e,xhat),nv)*CEout*(-1J*ka)
            gf = LinearForm(fespace)
            gf += SymbolicLFI(func*w.Trace(),BND,definedon=
                              mesh.Boundaries(collection))
            ## Functional for the second part using finite element cutoff
            eT=cross(xhat,Ocross(e,xhat))
            vv.Set(eT*CEout,BND,definedon=mesh.Boundaries(collection))
            gf += SymbolicLFI(curl(vv)*curl(w)-ka*ka*vv*w,VOL,
                              definedon=mesh.Materials('air'))
            with TaskManager():
                gf.Assemble()
            # Conjugates are needed since inner product conjugates the
            # second component(?)
            ff1a=gf.vec.InnerProduct(Es.vec)
            ff[j]=ff1a.conjugate()/4./pi
        FF.append(ff)
        ang.append(theta)
        #print('Done ',na,' of ',nang-1)
    print(ka,file=fff)
    print(nang,file=fff)
    print(dv[0], dv[1], dv[2], file=fff)
    print(pv[0], pv[1], pv[2], file=fff)
    print(idir,file=fff)
    for j in range(0,nang):
        F1=FF[j][0]
        F2=FF[j][1]
        F3=FF[j][2]
        theta=ang[j]
        print("%12.12e %12.12e %12.12e %12.12e %12.12e %12.12e %12.12e"%(theta,
                F1.real,F1.imag,F2.real,F2.imag,F3.real,F3.imag),file=fff)
    return(ang,FF)

def findff(FFP,rad,E,Ein,CurlEin,phii):
    print('Using findff...')
    # Specialized code for the impendace problem 
    nEs=cross(nv,E.components[0])
    FF=[];
    np=0
    for p in FFP:
        xhat = p.p
        nxhat=sqrt(xhat[0]**2+xhat[1]**2+xhat[2]**2)
        if nxhat>1-.0001:  # the ff mesh contains some interior points
            points.append(p.p)
        
    for xhat in points:
        print("p: ",xhat)
        Eout = exp(-1J*k*(x*xhat[0]+y*xhat[1]+z*xhat[2]))
        ff0=[0.,0.,0.*1J]
        for j in range(0,3):
            I1 = Integrate(nEs[j]*Eout,mesh,BND,order+1,region_wise=True)
            ff0[j]=I1[0];
        #
        phis=GridFunction(fesH1Bnd)
        phis=E.components[1]
        gphis=cross(cross(nv,grad(phis)),nv)
        nEsn=cross(nEs,nv)
        gphii=cross(cross(nv,grad(phii)),nv)
        nEin=cross(nv,cross(Ein,nv))
        ff1=[0,0,0.*1J]
        G=gphis+gphii
        SE=nEsn+nEin
        for j in range(0,3):
            I1a=Integrate(SE[j]*Eout,mesh,BND,order+1,region_wise=True)
            I1b=Integrate(G[j]*Eout,mesh,BND,order+1,region_wise=True)
            ff1[j]=(l/1J/k)*(I1a[0]+I1b[0])
        ff1=Ocross(ff1,xhat)
        #
        nCurlEin=cross(nv,CurlEin)
        ff2=[0,0,0.*1j]
        for j in range(0,3):
            I1=Integrate(nCurlEin[j]*Eout,mesh,BND,order+1,region_wise=True)
            ff2[j]=-(1/1J/k)*I1[0]
        ff2=Ocross(ff2,xhat)
        #
        pi=4*atan(1.)
        ff=Ocross(xhat,ff0+ff1+ff2)
        ff=[(1J*k/4/pi)*ffi for ffi in ff]
        FF.append(ff)
    return(FF)

def getFF(filename,itype):
    ## Read in the mesh and multistatic far field pattern
    print('Reading data from: ',filename)
    ## First the mesh
    fim=open('ffmesh.msh','r');
    npoints=int(fim.readline())
    print('npoints',npoints)
    Points = np.zeros((3,npoints))
    for j in range(0,npoints):
        line=fim.readline()
        columns=line.split()
        #print(columns)
        Points[:,j]=np.array([float(columns[0]),float(columns[1]),
                              float(columns[2])])
        #print(Points[:,j])
    nb=0
    for j in range(0,npoints):
        if np.sqrt(np.linalg.norm(Points[:,j])>1-1.e-6):
            nb=nb+1
    #print('nb=',nb)
    ntris=int(fim.readline())
    #print('ntris=',ntris)
    tris=np.zeros((3,ntris))
    for j in range(0,ntris):
        line=fim.readline()
        columns=line.split()
        tris[:,j]=np.array([int(columns[0]),int(columns[1]),
                              int(columns[2])])
        #print(tris[:,j])
    fim.close()
    ## Now the far field pattern
    fi=open(filename,'r')
    print('Loading ',filename)
    ka=float(fi.readline())
    print('k=',ka)
    if itype=='L':
        lam=float(fi.readline())
    else:
        lam=np.NaN
    indexes=np.zeros((2,2*nb))
    dv=np.zeros((3,2*nb))
    pv=np.zeros((3,2*nb))
    FFmat=np.zeros((nb,2,3,npoints),dtype='complex')
    for ipoints in range(0,nb):
        for ip in range(0,2):
            line=fi.readline().split()
            indexes[:,2*ipoints+ip]=[int(num) for num in line]
            line=fi.readline().split()
            dv[:,2*ipoints+ip]=[float(num) for num in line]
            line=fi.readline().split()
            pv[:,2*ipoints+ip]=[float(num) for num in line]
            #print('ind:',indexes[:,2*ipoints+ip])
            #print('dv:',dv[:,2*ipoints+ip])
            #print('pv:',pv[:,2*ipoints+ip])
            for j in range(0,npoints):
                line=fi.readline().split()
                row=[float(num) for num in line]
                crow=[row[0]+1J*row[1],row[2]+1J*row[3],row[4]+1J*row[5]]
                FFmat[ipoints,ip,:,j]=crow
    fi.close()
    #print(FFmat)
    return(Points,npoints,tris,ntris,indexes,dv,pv,nb,ka,lam,FFmat)

# function pltsphere(Points,npoints,tris,ntris)
# Lmax=max(max(Points))
# Lmin=min(min(Points))
# axis([-1.5,1.5,-1.5,1.5,-1.5,1.5])
# hold on
# for j=1:ntris
#     %sqrt(Points(1,tris(1,j))^2+Points(2,tris(1,j))^2+Points(3,tris(1,j))^2)-1
#     patch(Points(1,tris(:,j)),Points(2,tris(:,j)),Points(3,tris(:,j)),'y')
# end
# axis('square')
# hold off
# end


# function pltvec(Points,npoints,tris,ntris,FFmat,dvec,pvec)
# ip=1
# dv=0.5*dvec;
# pv=0.5*pvec
# for j=1:npoints
#     if sqrt(Points(:,j)'*Points(:,j))>1-1.e-6
#         pltx(ip)=Points(1,j);
#         plty(ip)=Points(2,j);
#         pltz(ip)=Points(3,j);
#         u(ip)=real(FFmat(1,j));
#         v(ip)=real(FFmat(2,j));
#         w(ip)=real(FFmat(3,j));
#         ip=ip+1;
#     end
# end
# hold on
# quiver3(pltx,plty,pltz,u,v,w,'b')
# h=line([2*dv(1),3*dv(1)],[2*dv(2),3*dv(2)],[2*dv(3),3*dv(3)]);
# set(h,'Color','r','LineWidth',2)
# h=line([2*dv(1),2*dv(1)+pv(1)],[2*dv(2),2*dv(2)+pv(2)],...
#     [2*dv(3),2*dv(3)+pv(3)]);
# set(h,'Color','g','LineWidth',2)
# hold off
# view(dvec)
# drawnow
# end

# function pltnrm(Points,npoints,tris,ntris,FFmat,dvec,pvec)
# FFn=sqrt(abs(FFmat(1,:)).^2+abs(FFmat(2,:)).^2+abs(FFmat(3,:)).^2)
# ii=find(FFn<1.e4);
# FFmin=min(FFn)
# FFmax=max(FFn(ii))
# %FFmax=1
# %axis([-2,2,-2,2,-2,2])
# axis([-1,1,-1,1,-1,1])
# hold on
# for j=1:ntris
#     p1=Points(1,tris(:,j)).*(1+(FFn(tris(:,j))-FFmin)/(FFmax-FFmin));
#     p2=Points(2,tris(:,j)).*(1+(FFn(tris(:,j))-FFmin)/(FFmax-FFmin));
#     p3=Points(3,tris(:,j)).*(1+(FFn(tris(:,j))-FFmin)/(FFmax-FFmin));
#     p1=Points(1,tris(:,j)).*FFn(tris(:,j))/FFmax;
#     p2=Points(2,tris(:,j)).*FFn(tris(:,j))/FFmax;
#     p3=Points(3,tris(:,j)).*FFn(tris(:,j))/FFmax;
#     patch(p1,p2,p3,'y')
# end
# axis('square')
# view(dvec)
# drawnow
# end

#     function []=readff(plt)
