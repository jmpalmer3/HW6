%% MECE5397: Homework Assignment #6
% Name: Jordan Palmer
% Last Modified: 11/17/2020

%% Crank-Nicolson Met

% 0 <= t <= T
% 0 <= x <= L

% House keeping commands
clc
clear
close all

%test conditions 1 
L = pi;     %space
T = 10;     %time
n = 10;
dt = T/n;
dx = L/n;
k = 1;
D = 0.1;    %coefficient
%fx = sin(k*x);
g0 = 0;
gL = g0;
U0 = zeros(1, n);
lam = D*dt/(dx)^2;
uexact = @(x) exp(-D*k^2)*sin(k*x);


%boundary conditions
U0(1) = g0;

for x = 2:n-1
    U0(x) = sin(k*x*dx);
end
U0(n) = gL;
Ugraph = U0;
Un(1)=g0; 
Un(n)=gL;

%
%[ a -b  0  0 
% -c  a -b  0 
%  0 -c  a -b
%  0  0 -c  a ] 

%time loop 
for t=2:n
    for i=1:n-2 
        if i == 1
            f(i)=lam*U0(i)+2*(1-lam)*U0(i+1)+lam*U0(i+2)+lam*Un(i);
        elseif i == n-2
            f(i)=lam*U0(i)+2*(1-lam)*U0(i+1)+lam*U0(i+2)+lam*Un(n);
        else
            f(i)=lam*U0(i)+2*(1-lam)*U0(i+1)+lam*U0(i+2);
        end
    end
    
    %set up coefficients
    b = lam*ones(n-3,1);
    c = b;
    a = (2*(1+lam))*ones(n-2,1);
    %put into matrix
    matrix = diag(a)+ diag(-b,1)+ diag(-c,-1);
    
    %divide by function values
    Ufinal = matrix\f';
    %make vector for new row
    Un=[Un(1),Ufinal',Un(n)];
    %Add to graphing value
    Ugraph(t,:)=Un;
    
    %reset U0 with new value
    U0 = Un;
end

mesh(t,x,u)
% axis([0,L,0,T])
% grid on
