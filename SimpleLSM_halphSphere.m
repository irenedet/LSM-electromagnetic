%
% Simple LSM with Ridge regression
%
clear all
close all
% Read data from FFP.txt to build prediction matrix
data_folder = '/home/papalotl/Dropbox/LSM-data/Data_3.0_1';
[Points,npoints,tris,ntris,indexes,dv,pv,nb,ka,FFmat]=...
    readff(1,data_folder);

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
disp(['A has been assembled, and cond(A)=',num2str(cond(A))])
%% 2. Noise matrix (to avoid numerical crimes and account for measurement errors)
eps_noise=0.000001;
Noise_mat=eps_noise*(-eye(size(A))+2*rand(size(A)));% uniform random noise
Anoisy=A+Noise_mat;
noise=norm(Noise_mat)/norm(A);
disp(['The level of noise in this experiment is: ',num2str(100*noise),'%'])
prompt = 'Do you want to continue wih this level of noise? (y/n): ';
str = input(prompt,'s');
if (str == 'y')
    A=Anoisy;
end
%% 3. Assemble right-hand-side
[SSpoints,Nsample,RHSmat]=readrhs(1,nb,data_folder);
disp(['Using ',num2str(Nsample),' sampling points in each direction'])
xx=SSpoints(1,:);
yy=SSpoints(2,:);
zz=SSpoints(3,:);
%scatter3(xx,yy,zz,'*')
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
tol=max(ginv)/2;
for j=1:Nsample
    if ginv(j)>tol
        indicator(:,j)=[1,0,0]';
    end
end
%% Plot interface
figure
hold on
S=5*ones(1,Nsample);
scatter3(xx,yy,zz,S,indicator');
%% Ploting the surface
R=0.6;
[theta,phi] = meshgrid(linspace(0,2*pi,32),linspace(-pi/2,0,32));
x = R.*cos(theta).*cos(phi);
y = R.*sin(theta).*cos(phi);
z = R.*sin(phi);
h1 = surf(x,y,z)
set(h1,'edgecolor','none','FaceColor',[.1,.5,.5])
hold on
[r, theta]=meshgrid(linspace(0,R,32),linspace(0,2*pi,32));
xs = r.*cos(theta);
ys = r.*sin(theta);
zs = r.*sin(0);
h2 = surf(xs,ys,zs)
set(h2,'edgecolor','none','FaceColor',[.1,.5,.5])
%% 4. Plotting the delamination (defective points)
zcoord=SSpoints(3,:);
flatside=find((zcoord==0));Flatside=SSpoints(:,flatside);
xi=50; %amplitude
func=ginv;
ginvflat=zeros(size(Flatside));ginvflat(3,:)=xi*func(flatside);Ginvx1=Flatside+ginvflat;
hold on
scatter3(Ginvx1(1,:),Ginvx1(2,:),Ginvx1(3,:),S(1:length(Ginvx1)),'b');


