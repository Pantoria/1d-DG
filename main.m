% Purpose: This program try to solve the one dimensional advection problem using Discontious Galerkin(DG)
% The eqation is Ut+a*Ux=0
% Boundary condition(BC) is U(0,t) = -sin(2*pi*t);
% Initial condition is U(x,0) = sin(x);
% The solve region is T*X = [0,1]*[0,2*pi]

clear all
clc

%=======================
%This block define the constant for calculation 
N=16; %the order of the Legendre polynomial,  using the Legendre polynomial on each element

K=50; %the element number of the space
X=2*pi; %the total length of the space to simulate 
hk = X/K; %the size of each element,the total lengthi is 2*pi

timestep = 0.01;
T = 1;
T_num = T/timestep+1; %The grid number in the time domain 

u_coefficient  = zeros(N,K,T_num);%initialize the coefficient for each element at every timestep

a = 2*pi; %The parameter in the 


x = -1:0.01:1;%This is used for the integral in the interval of -1 and 1
t = 0:timestep:T; %inintialize the time sequence for simulation 

%This block initial the initialize the coefficient before Legendre polynomial for each element
for i =1:K
	for n = 1:N
		y = sin(hk/2*(x+1)+(i-1)*hk).*Legendre(x,n)';  % using the initial condition 
		u_coefficient(n,i,1) = trapz(x,y); %get the coefficient using the orthogonality and normality 
	end
end

%initialize the numerical flux
num_flux = zeros(K+1,T_num); %This value define the numerical flux in each vertex of the element
num_flux(1,:) = -a*sin(a*t); %initialize the numerical flux at point x according to the boudnary ondition

%This block define the left side and right side of different order Legendre polynomial
LegRight=zeros(N,1);
for i = 1:N
	LegRight(i) = Legendre(1,i);
end
LegLeft=zeros(N,1);
for i = 1:N
	LegLeft(i) = Legendre(-1,i);
end

%This calculate the stiff matrix for the calculation 
StiffMatrix =zeros(N,N);
for j = 1:N
	for k = 1:N
		y = GradLegendre(x,j).*Legendre(x,k);
		StiffMatrix(j,k)  = a*trapz(x,y);
		StiffMatrix(j,k) = StiffMatrix(j,k)-a*Legendre(1,j)*Legendre(1,k);
	end
end

%Left is the matrix that used for the calculation after the discretization of time domain 
Left = eye(N,N)-2*timestep/hk*StiffMatrix;
 
 %******************************
 %This is the core algorithm                    
 for j = 1:K %This loop is for the time domain
	for i = 2:T_num %This loop if for the space domain
		u_coefficient(:,j,i) = inv(Left)*(u_coefficient(:,j,i-1)+2*timestep/hk*num_flux(j,i)*LegLeft);
		num_flux(j+1,i) = a*dot(u_coefficient(:,j,i),LegRight); %update the numerical flux in next time step
	end
end


%==================
%This block is used for the plot 
space = 0:X/1000:X-X/1000;
t = 0:timestep:T-timestep;
%calculate the u value at each time point and space point using the coefficient that we get 
for i = 1:100 %represent the time domain index 
	for j = 1:1000 %represent the space domain index 
		value = 0;
		elem_num = ceil(j*X/1000/hk);
		for k = 1:N
			value = value + u_coefficient(k,elem_num,i)*Legendre(2*(j*X/1000-(elem_num-1)*hk)-1,k);
		end
		u(i,j) = value;
		u_real(i,j)=sin(j*X/1000-a*i*timestep); %This is the analytical solution for this equation
		error(i,j) = abs(u(i,j)-u_real(i,j));
	end
end

subplot(1,2,1);
mesh(space,t,u);
ylabel('Time','FontWeight','bold');
xlabel('Displacement','FontWeight','bold');
zlabel('Density','FontWeight','bold');
title('The Discontious Galerkin(DG) solution in time and space domain');

subplot(1,2,2);
mesh(space,t,u_real);
ylabel('Time','FontWeight','bold');
xlabel('Displacement','FontWeight','bold');
zlabel('Density','FontWeight','bold');
title('The analytical solution in time and space domain');

mesh(space,t,error);
ylabel('Time','FontWeight','bold');
xlabel('Displacement','FontWeight','bold');
zlabel('Error','FontWeight','bold');
title('The error between analytical and numerical solution');