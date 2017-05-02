%
%  Simple LSM with fixed regularization parameter
%
clear all
close all
% Read data from FFP.txt
[Points,npoints,tris,ntris,indexes,dv,pv,nb,ka,FFmat]=readff(1);

A=zeros(2*nb,2*nb);
b=zeros(2*nb,1);

% Compute weights w for quadrature on the surface of the sphere
w=zeros(nb,1);
Stot=0;
for j=1:ntris
    AB=Points(:,tris(2,j))-Points(:,tris(1,j));
    AC=Points(:,tris(3,j))-Points(:,tris(1,j));
    S=norm(cross(AB,AC))/2;
    Stot=Stot+S;
    w(tris(1,j))= w(tris(1,j))+S/3;
    w(tris(2,j))= w(tris(2,j))+S/3;
    w(tris(3,j))= w(tris(3,j))+S/3;
end

% Fill matrix
for m=1:nb
    FF1=w(m)*FFmat{m,1};
    FF2=w(m)*FFmat{m,2};
    for l=1:nb
        A(l,m)=pv(:,2*l-1)'*FF1(:,l);%M_thetatheta
        A(nb+l,m)=pv(:,2*l)'*FF1(:,l);%M_phitheta
        A(l,m+nb)=pv(:,2*l-1)'*FF2(:,l);%M_thetaphi
        A(nb+l,nb+m)=pv(:,2*l)'*FF2(:,l);%M_phiphi
    end
end
%%
[SSpoints,Nsample,RHSmat]=readrhs(1,nb);
disp(['Using ',num2str(Nsample),' sampling points in each direction'])
xx=SSpoints(1,:);
yy=SSpoints(2,:);
zz=SSpoints(3,:);
scatter3(xx,yy,zz,'*')
g=zeros(Nsample,1);
% Simple fixed Tikhonov parameter suitable for Stekloff problem
alpha=1.e-8;
disp(['Tikhonov parameter: ',num2str(alpha)])
for j=1:Nsample
    z=SSpoints(:,j);
    for iq=1:2
        q=[0;0;0];
        q(iq)=1;
%         if (z(1)==0.6)||(z(1)==-0.6) 
%             if (iq ==1)
%                 vec=[1,0,0];
%             elseif (iq==2)
%                 vec=[0,1,0];
%             end
%         elseif (z(2)==0.6)||(z(2)==-0.6)
%             if (iq ==1)
%                 vec=[1,0,0];
%             elseif (iq==2)
%                 vec=[0,0,1];
%             end
%         elseif (z(3)==0.6)||(z(3)==-0.6)
%             if (iq==1)
%                 vec=[0,0,1];
%             elseif (iq==2)
%                 vec=[0,1,0];
%             end
%         end
        
    
        for l=1:nb
%             b(l)=q'*RHSmat{l,1}(:,j);% q'*E_b(z_j,-dv_l,pv1)
%             b(l+nb)=q'*RHSmat{l,2}(:,j);
            vec=cross(dv(:,2*l-1),cross(q,dv(:,2*l-1)));%this has to be changed for non-homogeneous background!
            b(l)=1j*ka*(pv(:,2*l-1)'*vec)*exp(-1i*ka*(dv(:,2*l-1)'*z));
            vec=cross(dv(:,2*l),cross(q,dv(:,2*l)));
            b(l+nb)=1j*ka*(pv(:,2*l)'*vec)*exp(-1i*ka*(dv(:,2*l)'*z));
        end
        ball(:,iq)=b;
    end
    g3=(A'*A+alpha*eye(size(A)))\(A'*ball);
    ginv(j)=1/norm(g3(:,1))+1/norm(g3(:,2));%+1/norm(g3(:,3));
        
end
indicator=zeros(3,Nsample);
tol=max(ginv)/3;
for j=1:Nsample
    if ginv(j)>tol
        indicator(:,j)=[1,0,0]';
    end
end
%C=repmat()
S=10*ones(1,Nsample);
scatter3(xx,yy,zz,S,indicator');
%% Plotting
% h=patch(x(:,i),y(:,i),z(:,i),'k')
% s=0.6;
% x=[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]*s;
% y=[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]*s;
% z=[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]*s;
% for i=1:6
%     h=patch(x(:,i),y(:,i),z(:,i),'k');
%     set(h,'edgecolor','w')
% end
l=0.6*ones(1,Nsample);
nosidex1=SSpoints(1,:)-l;nosidex2=SSpoints(1,:)+l;
nosidey1=SSpoints(2,:)-l;nosidey2=SSpoints(2,:)+l;
nosidez1=SSpoints(3,:)-l;nosidez2=SSpoints(3,:)+l;
sidex1=find((nosidex1==0));Sidex1=SSpoints(:,sidex1);
sidex2=find((nosidex2==0));Sidex2=SSpoints(:,sidex2);
sidey1=find((nosidey1==0));Sidey1=SSpoints(:,sidey1);
sidey2=find((nosidey2==0));Sidey2=SSpoints(:,sidey2);
sidez1=find((nosidez1==0));Sidez1=SSpoints(:,sidez1);
sidez2=find((nosidez2==0));Sidez2=SSpoints(:,sidez2);
for i=1:length(sidex1)
     h=patch(Sidex1(1,:),Sidex1(2,:),Sidex1(3,:),'k');
     set(h,'edgecolor','w')
 end
%%
% Lots of test points
nxx=41;
%lsnxx=21
disp(['Using ',num2str(nxx),' sampling points in each direction'])
xx=linspace(-1.5,1.5,nxx);
[zX,zY,zZ]=meshgrid(xx,xx,xx);
zX=reshape(zX,nxx^3,1,1);
zY=reshape(zY,nxx^3,1,1);
zZ=reshape(zZ,nxx^3,1,1);
g=zeros(nxx^3,1);
% Simple fixed Tikhonov parameter suitable for Stekloff problem
alpha=1.e-8;
disp(['Tikhonov parameter: ',num2str(alpha)])

% ball=zeros(2*nb,3);    
ball=zeros(2*nb,2);
for j=1:length(zX)
    if mod(j,nxx*nxx)==0, disp(['Done ',num2str(j), ' of ',num2str(length(zX))]),end
    z=[zX(j);zY(j);zZ(j)];
    for iq=1:3
        q=[0;0;0];
        q(iq)=1;
        for l=1:nb
            vec=cross(dv(:,2*l-1),cross(q,dv(:,2*l-1)));%this has to be changed for non-homogeneous background!
            b(l)=1i*ka*(pv(:,2*l-1)'*vec)*exp(-1i*ka*(dv(:,2*l-1)'*z));
            vec=cross(dv(:,2*l),cross(q,dv(:,2*l)));
            b(l+nb)=1i*ka*(pv(:,2*l)'*vec)*exp(-1i*ka*(dv(:,2*l)'*z));
        end
        ball(:,iq)=b;
    end
    % Solves for three polarizations of artificial source at once
    g3=(A'*A+alpha*eye(size(A)))\(A'*ball);
    ginv(j)=1/norm(g3(:,1))+1/norm(g3(:,2))+1/norm(g3(:,3));
end
ginv=reshape(ginv,nxx,nxx,nxx);
zX=reshape(zX,nxx,nxx,nxx);
zY=reshape(zY,nxx,nxx,nxx);
zZ=reshape(zZ,nxx,nxx,nxx);

level = 0.01;
figure
slice(zX,zY,zZ,ginv,0,0,0)
hold on
isosurface(zX,zY,zZ,ginv,level*(max(max(max(ginv)))-min(min(min(ginv))))+min(min(min(ginv))))
hold off


