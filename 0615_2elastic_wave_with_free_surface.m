%% ------------------- 1. model input -------------------
clear all; close all;

dt = 0.0001;           % 时间步长，单位s
tt = -0.02:dt:0.02;    %主频为60HZ
fm = 60;
A = 1;
% wave=A*(1-2*(pi*fm*tt).^2).*exp(-(pi*fm*tt).^2);
% wave = -A*tt.*exp(-(pi*fm*tt).^2);
wave = -A*exp(-(pi*fm*tt).^2);
wave=[0 0 wave];
T = 0.5;              % 波动传播时间，单位s
wave(round(T/dt))=0;    % 将子波后面部分补零

load '0626_model_multiples.mat'
% Vv=importdata('V1.mat');      % 纵波速度,m/s
% Vs=importdata('V2.mat');      % 横波速度,m/s
V2(:) = 0;
% den=importdata('density.mat');     % 介质密度，kg/m^3
Vv = V1; Vs = V2; den = density;

[Nz,Nx] = size(Vv);
pml = 50;                              % 吸收层的网格数
m = Nx + 2*pml;
n = Nz + 2*pml;
dz = 0.5;
dx = 0.5;
% d=2000*ones(n,m);                    % 介质密度,kg/m^3

V=zeros(n,m);                        % 介质纵波速度,m/s,补边界
V(pml+1:pml+Nz,pml+1:pml+Nx) = Vv;
V(pml+1:pml+Nz,1:pml) = repmat(Vv(:,1),[1 pml]);
V(pml+1:pml+Nz,1+Nx+pml:2*pml+Nx) = repmat(Vv(:,Nx),[1 pml]);
V(1:pml,:) = repmat(V(pml+1,:),[pml 1]);
V(1+pml+Nz:2*pml+Nz,:) = repmat(V(pml+Nz,:),[pml 1]);
clear Vv;
Vmax=max(max(V));         % 纵波最大速度

VS=zeros(n,m);                        % 介质横波速度,m/s,补边界
VS(pml+1:pml+Nz,pml+1:pml+Nx) = Vs;
VS(pml+1:pml+Nz,1:pml) = repmat(Vs(:,1),[1 pml]);
VS(pml+1:pml+Nz,1+Nx+pml:Nx+2*pml) = repmat(Vs(:,Nx),[1 pml]);
VS(1:pml,:) = repmat(VS(pml+1,:),[pml 1]);
VS(1+pml+Nz:2*pml+Nz,:) = repmat(VS(pml+Nz,:),[pml 1]);
clear Vs;

d=zeros(n,m);                        % 介质密度,kg/m^3,补边界
d(pml+1:pml+Nz,pml+1:pml+Nx) = den;
d(pml+1:pml+Nz,1:pml) = repmat(den(:,1),[1 pml]);
d(pml+1:pml+Nz,1+Nx+pml:Nx+2*pml) = repmat(den(:,Nx),[1 pml]);
d(1:pml,:) = repmat(d(pml+1,:),[pml 1]);
d(1+pml+Nz:2*pml+Nz,:) = repmat(d(pml+Nz,:),[pml 1]);
clear den;

R2(1:n,1:m) = d(1:n,1:m).*(V(1:n,1:m).^2)-2*d(1:n,1:m).*(VS(1:n,1:m).^2);
R1(1:n,1:m) = d(1:n,1:m).*(VS(1:n,1:m).^2);
D(1:n,1:m)  = R2(1:n,1:m)+2*R1(1:n,1:m);

%**********************衰减************************
R=1e-6;          % 理论反射系数
ddx=zeros(n,m); 
ddz=zeros(n,m);
plx=pml*dx;
plz=pml*dz;
 for i=1:n
     for k=1:m
         % 区域1
         if i>=1 && i<=pml && k>=1 && k<=pml
             x=pml-k;z=pml-i;
             ddx(i,k)=-log(R)*3*Vmax*x^2/(2*plx^2);
             ddz(i,k)=-log(R)*3*Vmax*z^2/(2*plz^2);
         elseif i>=1 && i<=pml && k>m-pml && k<=m
             x=k-(m-pml);z=pml-i;
             ddx(i,k)=-log(R)*3*Vmax*x^2/(2*plx^2);
             ddz(i,k)=-log(R)*3*Vmax*z^2/(2*plz^2);
         elseif i>n-pml && i<=n && k>=1 && k<=pml
             x=pml-k;z=i-(n-pml);
             ddx(i,k)=-log(R)*3*Vmax*x^2/(2*plx^2);
             ddz(i,k)=-log(R)*3*Vmax*z^2/(2*plz^2);
         elseif i>n-pml && i<=n && k>m-pml && k<=m
             x=k-(m-pml);z=i-(n-pml);
             ddx(i,k)=-log(R)*3*Vmax*x^2/(2*plx^2);
             ddz(i,k)=-log(R)*3*Vmax*z^2/(2*plz^2);
         % 区域2
         elseif i<=pml && k>pml && k<m-pml+1
             x=0;z=pml-i;
             ddx(i,k)=0;ddz(i,k)=-log(R)*3*Vmax*z^2/(2*plz^2);
         elseif  i>n-pml && i<=n && k>pml && k<=m-pml
             x=0;z=i-(n-pml);
             ddx(i,k)=0;ddz(i,k)=-log(R)*3*Vmax*z^2/(2*plz^2);
         % 区域3
         elseif i>pml && i<=n-pml && k<=pml
             x=pml-k;z=0;
             ddx(i,k)=-log(R)*3*Vmax*x^2/(2*plx^2);
             ddz(i,k)=0;
         elseif i>pml && i<=n-pml && k>m-pml && k<=m
             x=k-(m-pml);z=0;
             ddx(i,k)=-log(R)*3*Vmax*x^2/(2*plx^2);
             ddz(i,k)=0;
         end
     end
 end
% ddz(1:pml,1:m)=0; 
clear V;

% 炮检点信息
z0 = pml + 5;                       % 震源z位置
shot0 = pml+(30:3:120);
shotnum = length(shot0);  % 炮数
% shot0=round(m/2);
rz0 = pml + 1;
rx0 = pml+1:pml+Nx;
nt = round(T/dt);
record = zeros(nt,length(rx0),shotnum);
%% --------------------------- END 1 ---------------------------

%% -------------------- 2. Forward modeling --------------------
tic
for shot = 1%:shotnum
    
    x0 = shot0(shot);
    
    % **********************波场模拟********************
    disp(['开始第' num2str(shot) '炮正演计算,共' num2str(shotnum) '炮']);
    
    Pz=zeros(Nz,Nx);
    Px=zeros(Nz,Nx);
    
    pxx=zeros(n,m);
    pzz=zeros(n,m);
    pxz=zeros(n,m);
    Vx=zeros(n,m); 
    Vz=zeros(n,m); 
    pxx1=zeros(n,m);   
    pzz1=zeros(n,m);
    pxz1=zeros(n,m);
    Vx1=zeros(n,m); 
    Vz1=zeros(n,m); 
    pxx2=zeros(n,m);   
    pzz2=zeros(n,m);
    pxz2=zeros(n,m);
    Vx2=zeros(n,m); 
    Vz2=zeros(n,m);
    
    Seis1 = zeros(ceil(T/dt),Nx);
    Seis2 = zeros(ceil(T/dt),Nx);
    
    left = 0;right = 0;down = 0;radius = 0;
    for it=1:nt
        t = it*dt;
        pzz(z0,x0) = pzz(z0,x0) + wave(it);
        pxx(z0,x0) = pxx(z0,x0) + wave(it); 
        
        k = 4:m-4;
        i = 4:n-4;
        Vz1(i,k)=((1-0.5*dt*ddz(i,k)).*Vz1(i,k)+dt*((pzz(i+1,k)-pzz(i,k)).*1.171875-(pzz(i+2,k)-pzz(i-1,k)).*0.065104166666667...
             +(pzz(i+3,k)-pzz(i-2,k)).*0.0046875)./d(i,k)./dz)./(1+0.5*dt*ddz(i,k));
         Vz2(i,k)=((1-0.5*dt*ddx(i,k)).*Vz2(i,k)+dt*((pxz(i,k)-pxz(i,k-1)).*1.171875-(pxz(i,k+1)-pxz(i,k-2)).*0.065104166666667...
             +(pxz(i,k+2)-pxz(i,k-3)).*0.0046875)./d(i,k)./dx)./(1+0.5*dt*ddx(i,k));
         Vz(i,k) = Vz1(i,k)+Vz2(i,k);
         
         Vx1(i,k)=((1-0.5*dt*ddx(i,k)).*Vx1(i,k)+dt*((pxx(i,k+1)-pxx(i,k)).*1.171875-(pxx(i,k+2)-pxx(i,k-1))*0.065104166666667...
             +(pxx(i,k+3)-pxx(i,k-2)).*0.0046875)./d(i,k)./dx)./(1+0.5*dt*ddx(i,k));
         Vx2(i,k)=((1-0.5*dt*ddz(i,k)).*Vx2(i,k)+dt*((pxz(i,k)-pxz(i-1,k)).*1.171875-(pxz(i+1,k)-pxz(i-2,k))*0.065104166666667...
             +(pxz(i+2,k)-pxz(i-3,k)).*0.0046875)./d(i,k)./dz)./(1+0.5*dt*ddz(i,k));
         Vx(i,k) = Vx1(i,k)+Vx2(i,k);
         
%          Vz(z0,x0) = Vz(z0,x0) + wave(it);
         
         i = 4:pml;
         pzz(i,k) = -pzz((pml-i+pml),k);
         i = pml+1:n-4;
         pzz1(i,k)=((1-0.5*dt*ddz(i,k)).*pzz1(i,k)+D(i,k).*dt.*((Vz(i,k)-Vz(i-1,k)).*1.171875-(Vz(i+1,k)-Vz(i-2,k)).*0.065104166666667...
             +(Vz(i+2,k)-Vz(i-3,k)).*0.0046875)./dz)./(1+0.5*dt*ddz(i,k));
         pzz2(i,k)=((1-0.5*dt*ddx(i,k)).*pzz2(i,k)+R2(i,k).*dt.*((Vx(i,k)-Vx(i,k-1)).*1.171875-(Vx(i,k+1)-Vx(i,k-2)).*0.065104166666667...
             +(Vx(i,k+2)-Vx(i,k-3)).*0.0046875)./dx)./(1+0.5*dt*ddx(i,k));
         pzz(i,k) = pzz1(i,k) + pzz2(i,k);
         
         i = 4:n-4;
         pxx1(i,k)=((1-0.5*dt*ddx(i,k)).*pxx1(i,k)+D(i,k).*dt.*((Vx(i,k)-Vx(i,k-1)).*1.171875-(Vx(i,k+1)-Vx(i,k-2)).*0.065104166666667...
             +(Vx(i,k+2)-Vx(i,k-3)).*0.0046875)./dx)./(1+0.5*dt*ddx(i,k)); 
         pxx2(i,k)=((1-0.5*dt*ddz(i,k)).*pxx2(i,k)+R2(i,k).*dt.*((Vz(i,k)-Vz(i-1,k)).*1.171875-(Vz(i+1,k)-Vz(i-2,k)).*0.065104166666667...
             +(Vz(i+2,k)-Vz(i-3,k)).*0.0046875)./dz)./(1+0.5*dt*ddz(i,k)); 
         pxx(i,k) = pxx1(i,k)+pxx2(i,k);
         
         i = 4:pml;
         pxz(i,k) = -pxz((pml-i+pml),k);
         i = pml+1:n-4;
         pxz1(i,k)=((1-0.5*dt*ddx(i,k)).*pxz1(i,k)+R1(i,k).*dt.*((Vz(i,k+1)-Vz(i,k)).*1.171875-(Vz(i,k+2)-Vz(i,k-1)).*0.065104166666667...
             +(Vz(i,k+3)-Vz(i,k-2)).*0.0046875)./dx)./(1+0.5*dt*ddx(i,k));
         pxz2(i,k)=((1-0.5*dt*ddz(i,k)).*pxz2(i,k)+R1(i,k).*dt.*((Vx(i+1,k)-Vx(i,k)).*1.171875-(Vx(i+2,k)-Vx(i-1,k)).*0.065104166666667...
             +(Vx(i+3,k)-Vx(i-2,k)).*0.0046875)./dz)./(1+0.5*dt*ddz(i,k));
         pxz(i,k) = pxz1(i,k)+pxz2(i,k);
         
         Pz(1:n-2*pml,1:m-2*pml) = Vz(1+pml:n-pml,1+pml:m-pml);
         Px(1:n-2*pml,1:m-2*pml) = Vx(1+pml:n-pml,1+pml:m-pml);
         
         % 记录当前波场
         record(it,:,shot) = Vz(rz0,rx0);
         
         %………………波场动画显示程序……………………
%          subplot(1,2,1);
%          imagesc(Pz*10e9);
%          shading interp;
% %        axis square;
%          colormap('gray');
%          xlabel('Z分量');
%          subplot(1,2,2);
%          imagesc(Px*10e9);
%          shading interp;
%          %      axis square;
%          colormap('gray');
%          xlabel('X分量');
%          pause(1e-3);
         
         if rem(it,500)==1
             disp(['t=' num2str(t) ]);
         end
    end
    

end
toc
% save 'record1.mat' record
figure;mwigb(record(:,:,1),3);fig
%%
pic = record(:,:,1);
pic = aec(pic,dt,400*dt,200*dt);
figure;fig;imagesc(pic);colormap(gray)