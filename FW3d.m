%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Seis1, Seis2, Seis3]=FW3d(shotnum,location,Vv,Vs) 


tic
%%******************Ricker**************
dt=0.0005;    % time step,s
tt=-0.02:dt:0.02;                % change
fm=60;                           % change
A=1;
%wave=A*(1-2*(pi*fm*tt).^2).*exp(-(pi*fm*tt).^2);
wave=load('wavelet.mat');wave=wave.w1;
T=0.4;        % total time,s, change
wave(round(T/dt))=0;    % zeros fill
%%***********************************************


[Nz Nx Ny] = size(Vv);
pml=10;                              % pml layers,change
n = Nz + 2*pml;
m = Nx + 2*pml;
l = Ny + 2*pml;

dz=10;                              % grid step,m
dx=10;                              % grid step,m
dy=10;                              % grid step,m

d=2000*ones(n,m,l);                    % density,kg/m^3

%%**********************Vp************************%%
V=zeros(n,m,l);                        % P wave velocity, zeros fill, m/s
V(pml+1:pml+Nz,pml+1:pml+Nx,pml+1:pml+Ny) = Vv;
for i=1:pml
    for j=1:pml
    V(pml+1:pml+Nz,i,j) = Vv(:,1,1);
    V(pml+1:pml+Nz,i+Nx+pml,j) = Vv(:,Nx,1);
    V(pml+1:pml+Nz,i,j+Ny+pml) = Vv(:,1,Ny);
    V(pml+1:pml+Nz,i+Nx+pml,j+Ny+pml) = Vv(:,Nx,Ny);
    end
end
for i=1:Nx
    for j=1:pml
    V(pml+1:pml+Nz,i+pml,j) = Vv(:,i,1);
    V(pml+1:pml+Nz,i+pml,j+Ny+pml) = Vv(:,i,Ny);
    end
end
for j=1:Ny
    for i=1:pml
    V(pml+1:pml+Nz,i,j+pml) = Vv(:,1,j);
    V(pml+1:pml+Nz,i+Nx+pml,j+pml) = Vv(:,Nx,j);
    end
end
for r=1:pml
    for i=1:m
        for j=1:l
            V(r,i,j) = V(pml+1,i,j);
            V(r+Nz+pml,i,j) = V(pml+Nz,i,j);
        end
    end
end
clear Vv;
Vmax=0;         % max velocity of P wave
for r=1:n
    for i=1:m
        for j=1:l
        if V(r,i,j) > Vmax
            Vmax=V(r,i,j);
        end
        end
    end
end
%%**********************Vs************************%%
VS=zeros(n,m,l);                        % S wave velocity, zeros fill, m/s
VS(pml+1:pml+Nz,pml+1:pml+Nx,pml+1:pml+Ny) = Vs;
for i=1:pml
    for j=1:pml
    VS(pml+1:pml+Nz,i,j) = Vs(:,1,1);
    VS(pml+1:pml+Nz,i+Nx+pml,j) = Vs(:,Nx,1);
    VS(pml+1:pml+Nz,i,j+Ny+pml) = Vs(:,1,Ny);
    VS(pml+1:pml+Nz,i+Nx+pml,j+Ny+pml) = Vs(:,Nx,Ny);
    end
end
for i=1:Nx
    for j=1:pml
    VS(pml+1:pml+Nz,i+pml,j) = Vs(:,i,1);
    VS(pml+1:pml+Nz,i+pml,j+Ny+pml) = Vs(:,i,Ny);
    end
end
for j=1:Ny
    for i=1:pml
    VS(pml+1:pml+Nz,i,j+pml) = Vs(:,1,j);
    VS(pml+1:pml+Nz,i+Nx+pml,j+pml) = Vs(:,Nx,j);
    end
end
for r=1:pml
    for i=1:m
        for j=1:l
            VS(r,i,j) = VS(pml+1,i,j);
            VS(r+Nz+pml,i,j) = VS(pml+Nz,i,j);
        end
    end
end
clear Vs;

% lame
for r=1:n
   for i=1:m
       for j=1:l
       R2(r,i,j)=d(r,i,j).*(V(r,i,j).^2)-2*d(r,i,j).*(VS(r,i,j).^2);
       R1(r,i,j)=d(r,i,j).*(VS(r,i,j).^2);
       D(r,i,j)=R2(r,i,j)+2*R1(r,i,j);
       end
   end
end

%%**********************damping************************%%
%% ddx??ddz are damping coefficient of x and z
R=1e-6;          % theoretical reflection coefficient 
ddz=zeros(n,m,l); ddx=zeros(n,m,l); ddy=zeros(n,m,l); 

plz=pml*dz;      % thickness of PML
plx=pml*dx;
ply=pml*dy;

for r=1:n
    for i=1:m
        for j=1:l
        % x zone
        if i>=1 && i<=pml
            x=pml-i;
            ddx(r,i,j)=-log(R)*3*Vmax*x^2/(2*plx^2);
        end
        if i>=m-pml+1 && i<=m
            x=pml-(m-i);
            ddx(r,i,j)=-log(R)*3*Vmax*x^2/(2*plx^2);
        end
        % y zone
        if j>=1 && j<=pml
            y=pml-j;
            ddy(r,i,j)=-log(R)*3*Vmax*y^2/(2*ply^2);
        end
        if j>=l-pml+1 && j<=l
            y=pml-(l-j);
            ddy(r,i,j)=-log(R)*3*Vmax*y^2/(2*ply^2);
        end
        % z zone
        if r>=1 && r<=pml
            z=pml-r;
            ddz(r,i,j)=-log(R)*3*Vmax*z^2/(2*plz^2);
        end
        if r>=n-pml+1 && r<=n
            z=pml-(n-r);
            ddz(r,i,j)=-log(R)*3*Vmax*z^2/(2*plz^2);
        end
        end
    end
end

clear V;
                                                    %change
for shot = 1:shotnum                                  
% z0=pml+2;                       % change
% z0=pml+location(shotnum,1);
% x0=pml+location(shotnum,2);
% y0=pml+location(shotnum,3);


%% **************************************************
%% **********************????????*********************
z0=pml+location(shot,1);
x0=pml+location(shot,2);
y0=pml+location(shot,3);
shotid = ['Source No.' num2str(shot) ',the total number is ' num2str(shotnum) ',z0= ' num2str(z0) ',x0= ' num2str(x0) ',y0= ' num2str(y0) ''];
disp(shotid);

    
Pz=zeros(Nz,Nx,Ny);
Px=zeros(Nz,Nx,Ny);
Py=zeros(Nz,Nx,Ny);
 
pzz=zeros(n,m,l); pzz1=zeros(n,m,l); pzz2=zeros(n,m,l); pzz3=zeros(n,m,l);
pxx=zeros(n,m,l); pxx1=zeros(n,m,l); pxx2=zeros(n,m,l); pxx3=zeros(n,m,l); 
pyy=zeros(n,m,l); pyy1=zeros(n,m,l); pyy2=zeros(n,m,l); pyy3=zeros(n,m,l); 
pyz=zeros(n,m,l); pyz1=zeros(n,m,l); pyz2=zeros(n,m,l); 
pxz=zeros(n,m,l); pxz1=zeros(n,m,l); pxz2=zeros(n,m,l); 
pxy=zeros(n,m,l); pxy1=zeros(n,m,l); pxy2=zeros(n,m,l); 
Vz=zeros(n,m,l); Vz1=zeros(n,m,l); Vz2=zeros(n,m,l); Vz3=zeros(n,m,l); 
Vx=zeros(n,m,l); Vx1=zeros(n,m,l); Vx2=zeros(n,m,l); Vx3=zeros(n,m,l); 
Vy=zeros(n,m,l); Vy1=zeros(n,m,l); Vy2=zeros(n,m,l); Vy3=zeros(n,m,l); 

%Seis1=zeros(ceil(T/dt),Nx,Ny);
%Seis2=zeros(ceil(T/dt),Nx,Ny);
%Seis3=zeros(ceil(T/dt),Nx,Ny);

left = 0;right = 0;down = 0;radius = 0;

ttt = 0;
for t=dt:dt:T
    ttt = ttt + 1;
       pzz(z0,x0,y0)=pzz(z0,x0,y0)+wave(ttt);
       pxx(z0,x0,y0)=pxx(z0,x0,y0)+wave(ttt);
       pyy(z0,x0,y0)=pyy(z0,x0,y0)+wave(ttt);

%%................................ten order difference.............................%%
r=6:n-6;
i=6:m-6;
j=6:l-6;
            Vx1(r,i,j)=((1-0.5.*dt.*ddx(r,i,j)).*Vx1(r,i,j)+dt.*((pxx(r,i,j)-pxx(r,i-1,j)).*19845./16384-(pxx(r,i+1,j)-pxx(r,i-2,j)).*735./8192.0...
                       +(pxx(r,i+2,j)-pxx(r,i-3,j)).*567./40960-(pxx(r,i+3,j)-pxx(r,i-4,j)).*405./229376+(pxx(r,i+4,j)-pxx(r,i-5,j)).*35./294912)./d(r,i,j)./dx)./(1+0.5.*dt.*ddx(r,i,j));
            Vx2(r,i,j)=((1-0.5.*dt.*ddy(r,i,j)).*Vx2(r,i,j)+dt.*((pxy(r,i,j+1)-pxy(r,i,j)).*19845./16384-(pxy(r,i,j+2)-pxy(r,i,j-1)).*735./8192.0...
                       +(pxy(r,i,j+3)-pxy(r,i,j-2)).*567./40960-(pxy(r,i,j+4)-pxy(r,i,j-3)).*405./229376+(pxy(r,i,j+5)-pxy(r,i,j-4)).*35./294912)./d(r,i,j)./dy)./(1+0.5.*dt.*ddy(r,i,j));
            Vx3(r,i,j)=((1-0.5.*dt.*ddz(r,i,j)).*Vx3(r,i,j)+dt.*((pxz(r,i,j)-pxz(r-1,i,j)).*19845./16384-(pxz(r+1,i,j)-pxz(r-2,i,j)).*735./8192.0...
                       +(pxz(r+2,i,j)-pxz(r-3,i,j)).*567./40960-(pxz(r+3,i,j)-pxz(r-4,i,j)).*405./229376+(pxz(r+4,i,j)-pxz(r-5,i,j)).*35./294912)./d(r,i,j)./dz)./(1+0.5.*dt.*ddz(r,i,j));
            Vx(r,i,j)=Vx1(r,i,j)+Vx2(r,i,j)+Vx3(r,i,j);

            Vy1(r,i,j)=((1-0.5.*dt.*ddx(r,i,j)).*Vy1(r,i,j)+dt.*((pxy(r,i+1,j)-pxy(r,i,j)).*19845./16384-(pxy(r,i+2,j)-pxy(r,i-1,j)).*735./8192.0...
                       +(pxy(r,i+3,j)-pxy(r,i-2,j)).*567./40960-(pxy(r,i+4,j)-pxy(r,i-3,j)).*405./229376+(pxy(r,i+5,j)-pxy(r,i-4,j)).*35./294912)./d(r,i,j)./dx)./(1+0.5.*dt.*ddx(r,i,j));
            Vy2(r,i,j)=((1-0.5.*dt.*ddy(r,i,j)).*Vy2(r,i,j)+dt.*((pyy(r,i,j)-pyy(r,i,j-1)).*19845./16384-(pyy(r,i,j+1)-pyy(r,i,j-2)).*735./8192.0...
                       +(pyy(r,i,j+2)-pyy(r,i,j-3)).*567./40960-(pyy(r,i,j+3)-pyy(r,i,j-4)).*405./229376+(pyy(r,i,j+4)-pyy(r,i,j-5)).*35./294912)./d(r,i,j)./dy)./(1+0.5.*dt.*ddy(r,i,j));
            Vy3(r,i,j)=((1-0.5.*dt.*ddz(r,i,j)).*Vy3(r,i,j)+dt.*((pyz(r,i,j)-pyz(r-1,i,j)).*19845./16384-(pyz(r+1,i,j)-pyz(r-2,i,j)).*735./8192.0...
                       +(pyz(r+2,i,j)-pyz(r-3,i,j)).*567./40960-(pyz(r+3,i,j)-pyz(r-4,i,j)).*405./229376+(pyz(r+4,i,j)-pyz(r-5,i,j)).*35./294912)./d(r,i,j)./dz)./(1+0.5.*dt.*ddz(r,i,j));
            Vy(r,i,j)=Vy1(r,i,j)+Vy2(r,i,j)+Vy3(r,i,j);

            Vz1(r,i,j)=((1-0.5.*dt.*ddx(r,i,j)).*Vz1(r,i,j)+dt.*((pxz(r,i+1,j)-pxz(r,i,j)).*19845./16384-(pxz(r,i+2,j)-pxz(r,i-1,j)).*735./8192.0...
                       +(pxz(r,i+3,j)-pxz(r,i-2,j)).*567./40960-(pxz(r,i+4,j)-pxz(r,i-3,j)).*405./229376+(pxz(r,i+5,j)-pxz(r,i-4,j)).*35./294912)./d(r,i,j)./dx)./(1+0.5.*dt.*ddx(r,i,j));
            Vz2(r,i,j)=((1-0.5.*dt.*ddy(r,i,j)).*Vz2(r,i,j)+dt.*((pyz(r,i,j+1)-pyz(r,i,j)).*19845./16384-(pyz(r,i,j+2)-pyz(r,i,j-1)).*735./8192.0...
                       +(pyz(r,i,j+3)-pyz(r,i,j-2)).*567./40960-(pyz(r,i,j+4)-pyz(r,i,j-3)).*405./229376+(pyz(r,i,j+5)-pyz(r,i,j-4)).*35./294912)./d(r,i,j)./dy)./(1+0.5.*dt.*ddy(r,i,j));
            Vz3(r,i,j)=((1-0.5.*dt.*ddz(r,i,j)).*Vz3(r,i,j)+dt.*((pzz(r+1,i,j)-pzz(r,i,j)).*19845./16384-(pzz(r+2,i,j)-pzz(r-1,i,j)).*735./8192.0...
                       +(pzz(r+3,i,j)-pzz(r-2,i,j)).*567./40960-(pzz(r+4,i,j)-pzz(r-3,i,j)).*405./229376+(pzz(r+5,i,j)-pzz(r-4,i,j)).*35./294912)./d(r,i,j)./dz)./(1+0.5.*dt.*ddz(r,i,j));
            Vz(r,i,j)=Vz1(r,i,j)+Vz2(r,i,j)+Vz3(r,i,j);
     
            pxx1(r,i,j)=((1-0.5.*dt.*ddx(r,i,j)).*pxx1(r,i,j)+D(r,i,j).*dt.*((Vx(r,i+1,j)-Vx(r,i,j)).*19845./16384-(Vx(r,i+2,j)-Vx(r,i-1,j)).*735./8192.0...
                       +(Vx(r,i+3,j)-Vx(r,i-2,j)).*567./40960-(Vx(r,i+4,j)-Vx(r,i-3,j)).*405./229376+(Vx(r,i+5,j)-Vx(r,i-4,j)).*35./294912)./dx)./(1+0.5.*dt.*ddx(r,i,j));
            pxx2(r,i,j)=((1-0.5.*dt.*ddy(r,i,j)).*pxx2(r,i,j)+R2(r,i,j).*dt.*((Vy(r,i,j+1)-Vy(r,i,j)).*19845./16384-(Vy(r,i,j+2)-Vy(r,i,j-1)).*735./8192.0...
                       +(Vy(r,i,j+3)-Vy(r,i,j-2)).*567./40960-(Vy(r,i,j+4)-Vy(r,i,j-3)).*405./229376+(Vy(r,i,j+5)-Vy(r,i,j-4)).*35./294912)./dy)./(1+0.5.*dt.*ddy(r,i,j));
            pxx3(r,i,j)=((1-0.5.*dt.*ddz(r,i,j)).*pxx3(r,i,j)+R2(r,i,j).*dt.*((Vz(r,i,j)-Vz(r-1,i,j)).*19845./16384-(Vz(r+1,i,j)-Vz(r-2,i,j)).*735./8192.0...
                       +(Vz(r+2,i,j)-Vz(r-3,i,j)).*567./40960-(Vz(r+3,i,j)-Vz(r-4,i,j)).*405./229376+(Vz(r+4,i,j)-Vz(r-5,i,j)).*35./294912)./dz)./(1+0.5.*dt.*ddz(r,i,j));
            pxx(r,i,j)=pxx1(r,i,j)+pxx2(r,i,j)+pxx3(r,i,j);

            pyy1(r,i,j)=((1-0.5.*dt.*ddx(r,i,j)).*pyy1(r,i,j)+R2(r,i,j).*dt.*((Vx(r,i+1,j)-Vx(r,i,j)).*19845./16384-(Vx(r,i+2,j)-Vx(r,i-1,j)).*735./8192.0...
                       +(Vx(r,i+3,j)-Vx(r,i-2,j)).*567./40960-(Vx(r,i+4,j)-Vx(r,i-3,j)).*405./229376+(Vx(r,i+5,j)-Vx(r,i-4,j)).*35./294912)./dx)./(1+0.5.*dt.*ddx(r,i,j));
            pyy2(r,i,j)=((1-0.5.*dt.*ddy(r,i,j)).*pyy2(r,i,j)+D(r,i,j).*dt.*((Vy(r,i,j+1)-Vy(r,i,j)).*19845./16384-(Vy(r,i,j+2)-Vy(r,i,j-1)).*735./8192.0...
                       +(Vy(r,i,j+3)-Vy(r,i,j-2)).*567./40960-(Vy(r,i,j+4)-Vy(r,i,j-3)).*405./229376+(Vy(r,i,j+5)-Vy(r,i,j-4)).*35./294912)./dy)./(1+0.5.*dt.*ddy(r,i,j));
            pyy3(r,i,j)=((1-0.5.*dt.*ddz(r,i,j)).*pyy3(r,i,j)+R2(r,i,j).*dt.*((Vz(r,i,j)-Vz(r-1,i,j)).*19845./16384-(Vz(r+1,i,j)-Vz(r-2,i,j)).*735./8192.0...
                       +(Vz(r+2,i,j)-Vz(r-3,i,j)).*567./40960-(Vz(r+3,i,j)-Vz(r-4,i,j)).*405./229376+(Vz(r+4,i,j)-Vz(r-5,i,j)).*35./294912)./dz)./(1+0.5.*dt.*ddz(r,i,j));
            pyy(r,i,j)=pyy1(r,i,j)+pyy2(r,i,j)+pyy3(r,i,j);

            pzz1(r,i,j)=((1-0.5.*dt.*ddx(r,i,j)).*pzz1(r,i,j)+R2(r,i,j).*dt.*((Vx(r,i+1,j)-Vx(r,i,j)).*19845./16384-(Vx(r,i+2,j)-Vx(r,i-1,j)).*735./8192.0...
                       +(Vx(r,i+3,j)-Vx(r,i-2,j)).*567./40960-(Vx(r,i+4,j)-Vx(r,i-3,j)).*405./229376+(Vx(r,i+5,j)-Vx(r,i-4,j)).*35./294912)./dx)./(1+0.5.*dt.*ddx(r,i,j));
            pzz2(r,i,j)=((1-0.5.*dt.*ddy(r,i,j)).*pzz2(r,i,j)+R2(r,i,j).*dt.*((Vy(r,i,j+1)-Vy(r,i,j)).*19845./16384-(Vy(r,i,j+2)-Vy(r,i,j-1)).*735./8192.0...
                       +(Vy(r,i,j+3)-Vy(r,i,j-2)).*567./40960-(Vy(r,i,j+4)-Vy(r,i,j-3)).*405./229376+(Vy(r,i,j+5)-Vy(r,i,j-4)).*35./294912)./dy)./(1+0.5.*dt.*ddy(r,i,j));
            pzz3(r,i,j)=((1-0.5.*dt.*ddz(r,i,j)).*pzz3(r,i,j)+D(r,i,j).*dt.*((Vz(r,i,j)-Vz(r-1,i,j)).*19845./16384-(Vz(r+1,i,j)-Vz(r-2,i,j)).*735./8192.0...
                       +(Vz(r+2,i,j)-Vz(r-3,i,j)).*567./40960-(Vz(r+3,i,j)-Vz(r-4,i,j)).*405./229376+(Vz(r+4,i,j)-Vz(r-5,i,j)).*35./294912)./dz)./(1+0.5.*dt.*ddz(r,i,j));
            pzz(r,i,j)=pzz1(r,i,j)+pzz2(r,i,j)+pzz3(r,i,j);

            pyz1(r,i,j)=((1-0.5.*dt.*ddz(r,i,j)).*pyz1(r,i,j)+R1(r,i,j).*dt.*((Vy(r+1,i,j)-Vy(r,i,j)).*19845./16384-(Vy(r+2,i,j)-Vy(r-1,i,j)).*735./8192.0...
                       +(Vy(r+3,i,j)-Vy(r-2,i,j)).*567./40960-(Vy(r+4,i,j)-Vy(r-3,i,j)).*405./229376+(Vy(r+5,i,j)-Vy(r-4,i,j)).*35./294912)./dz)./(1+0.5.*dt.*ddz(r,i,j));
            pyz2(r,i,j)=((1-0.5.*dt.*ddy(r,i,j)).*pyz2(r,i,j)+R1(r,i,j).*dt.*((Vz(r,i,j)-Vz(r,i,j-1)).*19845./16384-(Vz(r,i,j+1)-Vz(r,i,j-2)).*735./8192.0...
                       +(Vz(r,i,j+2)-Vz(r,i,j-3)).*567./40960-(Vz(r,i,j+3)-Vz(r,i,j-4)).*405./229376+(Vz(r,i,j+4)-Vz(r,i,j-5)).*35./294912)./dy)./(1+0.5.*dt.*ddy(r,i,j));
            pyz(r,i,j)=pyz1(r,i,j)+pyz2(r,i,j);

            pxz1(r,i,j)=((1-0.5.*dt.*ddz(r,i,j)).*pxz1(r,i,j)+R1(r,i,j).*dt.*((Vx(r+1,i,j)-Vx(r,i,j)).*19845./16384-(Vx(r+2,i,j)-Vx(r-1,i,j)).*735./8192.0...
                       +(Vx(r+3,i,j)-Vx(r-2,i,j)).*567./40960-(Vx(r+4,i,j)-Vx(r-3,i,j)).*405./229376+(Vx(r+5,i,j)-Vx(r-4,i,j)).*35./294912)./dz)./(1+0.5.*dt.*ddz(r,i,j));
            pxz2(r,i,j)=((1-0.5.*dt.*ddx(r,i,j)).*pxz2(r,i,j)+R1(r,i,j).*dt.*((Vz(r,i,j)-Vz(r,i-1,j)).*19845./16384-(Vz(r,i+1,j)-Vz(r,i-2,j)).*735./8192.0...
                       +(Vz(r,i+2,j)-Vz(r,i-3,j)).*567./40960-(Vz(r,i+3,j)-Vz(r,i-4,j)).*405./229376+(Vz(r,i+4,j)-Vz(r,i-5,j)).*35./294912)./dx)./(1+0.5.*dt.*ddx(r,i,j));
            pxz(r,i,j)=pxz1(r,i,j)+pxz2(r,i,j);

            pxy1(r,i,j)=((1-0.5.*dt.*ddy(r,i,j)).*pxy1(r,i,j)+R1(r,i,j).*dt.*((Vx(r,i,j)-Vx(r,i,j-1)).*19845./16384-(Vx(r,i,j+1)-Vx(r,i,j-2)).*735./8192.0...
                       +(Vx(r,i,j+2)-Vx(r,i,j-3)).*567./40960-(Vx(r,i,j+3)-Vx(r,i,j-4)).*405./229376+(Vx(r,i,j+4)-Vx(r,i,j-5)).*35./294912)./dy)./(1+0.5.*dt.*ddy(r,i,j));
            pxy2(r,i,j)=((1-0.5.*dt.*ddx(r,i,j)).*pxy2(r,i,j)+R1(r,i,j).*dt.*((Vy(r,i,j)-Vy(r,i-1,j)).*19845./16384-(Vy(r,i+1,j)-Vy(r,i-2,j)).*735./8192.0...
                       +(Vy(r,i+2,j)-Vy(r,i-3,j)).*567./40960-(Vy(r,i+3,j)-Vy(r,i-4,j)).*405./229376+(Vy(r,i+4,j)-Vy(r,i-5,j)).*35./294912)./dx)./(1+0.5.*dt.*ddx(r,i,j));
            pxy(r,i,j)=pxy1(r,i,j)+pxy2(r,i,j);

%%...............................de-pml............................%%    
r=1:n-2*pml;
i=1:m-2*pml;
j=1:l-2*pml;
Pz(r,i,j)=Vz(r+pml,i+pml,j+pml);
Px(r,i,j)=Vx(r+pml,i+pml,j+pml);
Py(r,i,j)=Vy(r+pml,i+pml,j+pml);


if mod(ttt,10)==0
   time = ['Current time = ' num2str(t) ' of total ' num2str(T) 's'];disp(time);
end

Seis1(ttt+800*(shot-1),:,:)=Pz(1,1:Nx,1:Ny);
Seis2(ttt+800*(shot-1),:,:)=Px(1,1:Nx,1:Ny);
Seis3(ttt+800*(shot-1),:,:)=Py(1,1:Nx,1:Ny);
end



tips = 'calculation time:';disp(tips);
toc
end