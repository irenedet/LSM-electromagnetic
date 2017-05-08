function [Points,npoints,tris,ntris,indexes,dv,pv,nb,ka,FFmat]=readff(plt)

fidm=fopen('ffmesh.msh');
npoints=fscanf(fidm,'%i\n',1);
Points=fscanf(fidm,'%g %g %g\n',[3, npoints]);
nb=0;
% The ffmesh contains some interior points which need to be avoided
% The boundary points are first in the mesh
for j=1:npoints
    if sqrt(Points(:,j)'*Points(:,j))>1-1.e-6
        nb=nb+1;
        id(nb)=j; % not needed
    end
end
disp(['There are ',num2str(nb),' points on the boundary'])

ntris=fscanf(fidm,'%i\n',1);
tris=fscanf(fidm,'%i %i %i\n',[3, ntris]);
fclose(fidm);

if plt~=1
    figure(10)
    pltsphere(Points,npoints,tris,ntris)
end

fidff=fopen('FFP.txt','r');
FFmat=cell(nb,2);%FF matrices
ka=fscanf(fidff,'%f \n1',1)%wave number
delta=fscanf(fidff,'%f \n1',1)%thickness of thin layer
indexes = zeros(2,nb);

for ipoints=1:nb
    indexes(:,2*ipoints-1)=fscanf(fidff,'%i %i\n',2);
    dv(:,2*ipoints-1)=fscanf(fidff,'%f %f %f\n',3);
    pv(:,2*ipoints-1)=fscanf(fidff,'%f %f %f\n',3);
    FFtmp=fscanf(fidff,'%e %e %e %e %e %e\n',[6,nb]);%FFtmp=fscanf(fidff,'%e %e %e %e %e %e\n',[6,npoints]);
    FFmat{ipoints,1}=[FFtmp(1,:)+1i*FFtmp(2,:);FFtmp(3,:)+1i*FFtmp(4,:);...
        FFtmp(5,:)+1i*FFtmp(6,:)];
    if plt~=1
        figure(11)
        subplot(1,2,1)
        pltnrm(Points,npoints,tris,ntris,FFmat{ipoints,1},dv(:,2*ipoints-1),pv(:,2*ipoints-1))
        axis('square')
    end
    indexes(:,2*ipoints)=fscanf(fidff,'%i %i\n',2);
    dv(:,2*ipoints)=fscanf(fidff,'%f %f %f\n',3);
    pv(:,2*ipoints)=fscanf(fidff,'%f %f %f\n',3);
    FFtmp=fscanf(fidff,'%f %f %f %f %f %f\n',[6,nb]);
    FFmat{ipoints,2}=[FFtmp(1,:)+1i*FFtmp(2,:);FFtmp(3,:)+1i*FFtmp(4,:);...
        FFtmp(5,:)+1i*FFtmp(6,:)];
    if plt~=1
        figure(11)
        subplot(1,2,2)
        %pltsphere(Points,npoints,tris,ntris)
        %pltvec(Points,npoints,tris,ntris,FFmat{ipoints,2},dv(:,2*ipoints),pv(:,2*ipoints))
        pltnrm(Points,npoints,tris,ntris,FFmat{ipoints,2},dv(:,2*ipoints),pv(:,2*ipoints))
        axis('square')
        pause
        clf
    end
end
end

function pltsphere(Points,npoints,tris,ntris)
Lmax=max(max(Points))
Lmin=min(min(Points))
axis([-1.5,1.5,-1.5,1.5,-1.5,1.5])
hold on
for j=1:ntris
    %sqrt(Points(1,tris(1,j))^2+Points(2,tris(1,j))^2+Points(3,tris(1,j))^2)-1
    patch(Points(1,tris(:,j)),Points(2,tris(:,j)),Points(3,tris(:,j)),'y')
end
axis('square')
hold off
end


function pltvec(Points,npoints,tris,ntris,FFmat,dvec,pvec)
ip=1
dv=0.5*dvec;
pv=0.5*pvec
for j=1:npoints
    if sqrt(Points(:,j)'*Points(:,j))>1-1.e-6
        pltx(ip)=Points(1,j);
        plty(ip)=Points(2,j);
        pltz(ip)=Points(3,j);
        u(ip)=real(FFmat(1,j));
        v(ip)=real(FFmat(2,j));
        w(ip)=real(FFmat(3,j));
        ip=ip+1;
    end
end
hold on
quiver3(pltx,plty,pltz,u,v,w,'b')
h=line([2*dv(1),3*dv(1)],[2*dv(2),3*dv(2)],[2*dv(3),3*dv(3)]);
set(h,'Color','r','LineWidth',2)
h=line([2*dv(1),2*dv(1)+pv(1)],[2*dv(2),2*dv(2)+pv(2)],...
    [2*dv(3),2*dv(3)+pv(3)]);
set(h,'Color','g','LineWidth',2)
hold off
view(dvec)
drawnow
end

function pltnrm(Points,npoints,tris,ntris,FFmat,dvec,pvec)
FFn=sqrt(abs(FFmat(1,:)).^2+abs(FFmat(2,:)).^2+abs(FFmat(3,:)).^2);
ii=find(FFn<1.e4);
FFmin=min(FFn);
FFmax=max(FFn(ii));
%FFmax=1
%axis([-2,2,-2,2,-2,2])
axis([-1,1,-1,1,-1,1])
hold on
for j=1:ntris
    p1=Points(1,tris(:,j)).*(1+(FFn(tris(:,j))-FFmin)/(FFmax-FFmin));
    p2=Points(2,tris(:,j)).*(1+(FFn(tris(:,j))-FFmin)/(FFmax-FFmin));
    p3=Points(3,tris(:,j)).*(1+(FFn(tris(:,j))-FFmin)/(FFmax-FFmin));
    p1=Points(1,tris(:,j)).*FFn(tris(:,j))/FFmax;
    p2=Points(2,tris(:,j)).*FFn(tris(:,j))/FFmax;
    p3=Points(3,tris(:,j)).*FFn(tris(:,j))/FFmax;
    patch(p1,p2,p3,'y')
end
axis('square')
view(dvec)
drawnow
end
