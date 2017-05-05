%
% Simple LSM with Ridge regression
%
%clear all
%close all
% Read data from FFP.txt to build prediction matrix
% ||Ax-b||^2 + alpha ||x||^2
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
disp(['The prediction matrix A has been assembled.'])
%% 1.
disp(['cond(A)=',num2str(cond(A))])
%% 2. Noise matrix (to avoid numerical crimes and account for measurement errors)
eps_noise=0.0001;
Noise_mat=eps_noise*(-eye(size(A))+2*rand(size(A)));% uniform random noise
Anoisy=A+Noise_mat;
noise=norm(Noise_mat)/norm(A);
disp(['The level of noise in this experiment is: ',num2str(100*noise),'%'])
A=Anoisy;
%% 3. Assemble right-hand-side 
[SSpoints,Nsample,RHSmat]=readrhs(1,nb,folder_rhs);
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

        for l=1:nb
            b(l)=q'*RHSmat{l,1}(:,j);% q'*E_b(z_j,-dv_l,pv1)
            b(l+nb)=q'*RHSmat{l,2}(:,j);
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
S=10*ones(1,Nsample);
scatter3(xx,yy,zz,S,indicator');
%% Plotting the delamination (defective points)
%h=patch(x(:,i),y(:,i),z(:,i),'k')
figure
s=1.2;
x=[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]*s;
y=[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]*s;
z=[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]*s;
for i=1:6
    h=patch(x(:,i)-0.6,y(:,i)-0.6,z(:,i)-0.6,'k');
    set(h,'edgecolor','b','FaceColor',[.3,.5,.3])
end
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
%figure
hold on
% scatter3(Sidex1(1,:),Sidex1(2,:),Sidex1(3,:),'black');
% scatter3(Sidex2(1,:),Sidex2(2,:),Sidex2(3,:),'black');
% scatter3(Sidey1(1,:),Sidey1(2,:),Sidey1(3,:),'black');
% scatter3(Sidey2(1,:),Sidey2(2,:),Sidey2(3,:),'black');
% scatter3(Sidez1(1,:),Sidez1(2,:),Sidez1(3,:),'black');
% scatter3(Sidez2(1,:),Sidez2(2,:),Sidez2(3,:),'black');
ginvx1=zeros(size(Sidex1));ginvx1(1,:)=50*ginv(sidex1);Ginvx1=Sidex1+ginvx1;
ginvx2=zeros(size(Sidex2));ginvx2(1,:)=50*ginv(sidex2);Ginvx2=Sidex2-ginvx2;
ginvy1=zeros(size(Sidey1));ginvy1(2,:)=50*ginv(sidey1);Ginvy1=Sidey1+ginvy1;
ginvy2=zeros(size(Sidey2));ginvy2(2,:)=50*ginv(sidey2);Ginvy2=Sidey2-ginvy2;
ginvz1=zeros(size(Sidez1));ginvz1(3,:)=50*ginv(sidez1);Ginvz1=Sidez1+ginvz1;
ginvz2=zeros(size(Sidez2));ginvz2(3,:)=50*ginv(sidez2);Ginvz2=Sidez2-ginvz2;
scatter3(Ginvx1(1,:),Ginvx1(2,:),Ginvx1(3,:),'b');
scatter3(Ginvx2(1,:),Ginvx2(2,:),Ginvx2(3,:),'b');
scatter3(Ginvy1(1,:),Ginvy1(2,:),Ginvy1(3,:),'b');
scatter3(Ginvy2(1,:),Ginvy2(2,:),Ginvy2(3,:),'b');
scatter3(Ginvz1(1,:),Ginvz1(2,:),Ginvz1(3,:),'red');
scatter3(Ginvz2(1,:),Ginvz2(2,:),Ginvz2(3,:),'b');
%% !!!!
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
%% 4.
level = 0.4;
figure
slice(zX,zY,zZ,ginv,0,0,0)
%figure
hold on
s=1.2;
x=[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]*s;
y=[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]*s;
z=[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]*s;
for i=1:6
    h=patch(x(:,i)-0.6,y(:,i)-0.6,z(:,i)-0.6,'k');
    set(h,'edgecolor','b','FaceColor',[1,1,1])
end
hold on
isosurface(zX,zY,zZ,ginv,level*(max(max(max(ginv)))-min(min(min(ginv))))+min(min(min(ginv))))
hold off
