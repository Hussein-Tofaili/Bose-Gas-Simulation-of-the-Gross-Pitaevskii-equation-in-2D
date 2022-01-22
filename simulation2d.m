clearvars; close all; warning('off'); format long;
%% Defining variables

Fs = 10^6;
max=50*10^(-6);
X=-max:1/Fs:max;
Y=-max:1/Fs:max;
[x,y] = meshgrid(X,Y);


psi=y.^0;%initial wave function

N=10^4; %number of interacting atoms
h = 1.054571817*10^(-34); m = 1.44*10^(-25); 
wz=2*pi*4*10^3; w=2*pi*19.3; k = m*w^2; T=2*pi/w;
as=5*10^(-9); lz=(h/(m*wz))^0.5; g=sqrt(8*pi)*as/lz;

%% Potential

%Cube potential well 
V1=y*0; 
for i = 1:size(V1,1)
    if abs(x(1,i))> (27.6/2*10^(-6))
        for j = 1:size(V1,2)
            V1(j,i)=10^-28; %potential ceiling(in Joules)
        end
    end
    if abs(y(i)) > (27.6/2*10^(-6))
        for j = 1:size(V1,2)
            V1(i,j)=10^-28; %potential ceiling(in Joules)
        end
    end
end

V2=0.5*m*w^2*(x.^2+y.^2); %harmonic oscillator potential

%% Determining the ground state
V=V2; %choose between V1 (potential well) and V2 (harmonic potential)
t=0;
dt=10^(-5); %time steps in seconds

for i = 0:10000
    
    L = 4*del2(psi,1/Fs); %Laplacian of psi
    Ek=-h^2/(2*m)*L; %kinetic term
    
    Ein=h^2/m*g*N*psi.^3; %interaction term
    
    dpsi=-1/h*(Ek+Ein+V.*psi);
    psi = psi+dpsi*dt;
    
    
    F = trapz(Y,trapz(X,psi.*psi,2));%Normalization
    Norm = sqrt(1/F);
    psi = psi*Norm;

    t=t+dt;

end


%% Exact solution for harmonic oscillator ground state
sol=(m*w/(pi*h))^(1/2)*exp(-m*w*(x.^2+y.^2)/(2*h));

%% Plotting potential
figure(9)
surf(x,y,V*6.242*10^18) %multiply by 6.242*10^18 to tranfer from J to eV
xlabel('x (m)')
ylabel('y (m)')
zlabel('V (ev)')
title("Trapping Potential")

%% Plotting psi
figure(11)
surf(x,y,psi)
xlabel('x (m)')
ylabel('y (m)')
title("Ground state wave function")

