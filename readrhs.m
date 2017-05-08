function [SSpoints,Nsample,RHSmat]=readrhs(plt,nb,data_folder)

rhs_file = strcat(data_folder,'/RHS.txt');
Spoints_file = strcat(data_folder,'/Spoints.txt');
fidrhs=fopen(rhs_file,'r');
Nsample=fscanf(fidrhs,'%f \n1',1)%number of samling points
RHSmat=cell(nb,2);%FF matrices
fidsp=fopen(Spoints_file,'r');
SSpoints=zeros(3,Nsample);
for sp=1:Nsample
    SSpoints(:,sp)=fscanf(fidsp,'%f %f %f\n',3);
end
fclose(fidsp);
for ipoints=1:nb
    index=fscanf(fidrhs,'%i %i\n',2);
    dv=fscanf(fidrhs,'%f %f %f\n',3);
    pv=fscanf(fidrhs,'%f %f %f\n',3);
    RHStmp=fscanf(fidrhs,'%e %e %e %e %e %e\n',[6,Nsample]);
    RHSmat{ipoints,1}=[RHStmp(1,:)+1i*RHStmp(2,:);RHStmp(3,:)+1i*RHStmp(4,:);...
        RHStmp(5,:)+1i*RHStmp(6,:)];
    if plt~=1
        figure(11)
        subplot(1,2,1)
        pltnrm(Points,npoints,tris,ntris,RHSmat{ipoints,1},dv(:,2*ipoints-1),pv(:,2*ipoints-1))
        axis('square')
    end
    index=fscanf(fidrhs,'%i %i\n',2);
    dv=fscanf(fidrhs,'%f %f %f\n',3);
    pv=fscanf(fidrhs,'%f %f %f\n',3);
    RHStmp=fscanf(fidrhs,'%f %f %f %f %f %f\n',[6,Nsample]);
    RHSmat{ipoints,2}=[RHStmp(1,:)+1i*RHStmp(2,:);RHStmp(3,:)+1i*RHStmp(4,:);...
        RHStmp(5,:)+1i*RHStmp(6,:)];
    if plt~=1
        figure(11)
        subplot(1,2,2)
        %pltsphere(Points,npoints,tris,ntris)
        %pltvec(Points,npoints,tris,ntris,FFmat{ipoints,2},dv(:,2*ipoints),pv(:,2*ipoints))
        pltnrm(Points,npoints,tris,ntris,RHSmat{ipoints,2},dv(:,2*ipoints),pv(:,2*ipoints))
        axis('square')
        pause
        clf
    end
end
fclose(fidrhs);
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
