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


%boundary conditions
U0(1) = g0;

for x = 1:n-1
    U0(x+1) = sin(k*x*dx);
end
U0(n+1) = gL;
Ugraph = U0;
Un(1)=g0; 
Un(n+1)=gL;

%
%[ a -b  0  0 
% -c  a -b  0 
%  0 -c  a -b
%  0  0 -c  a ] 

%time loop 
for t=2:n+1
    for i=1:n-1 
        if i == 1
            f(i)=lam*U0(i)+2*(1-lam)*U0(i+1)+lam*U0(i+2)+lam*Un(i);
        elseif i == n-1
            f(i)=lam*U0(i)+2*(1-lam)*U0(i+1)+lam*U0(i+2)+lam*Un(n);
        else
            f(i)=lam*U0(i)+2*(1-lam)*U0(i+1)+lam*U0(i+2);
        end
    end
    
    %set up coefficients
    b = lam*ones(n-2,1);
    c = b;
    a = (2*(1+lam))*ones(n-1,1);
    %put into matrix
    matrix = diag(a)+ diag(-b,1)+ diag(-c,-1);
    
    %divide by function values
    Ufinal = matrix\f';
    %make vector for new row
    Un=[Un(1),Ufinal',Un(n+1)];
    %Add to graphing value
    Ugraph(t,:)=Un;
    
    %reset U0 with new value
    U0 = Un;
end



for t = 1:n+1
    for x = 1:n+1
        uexactFunction(t,x) = exp(-D*k^2*(t-1))*sin(k*(x-1)*dx);
    end
end

% uexact = @(x) exp(-D*k^2*1)*sin(k*x);
% x = (0:10/9:10)*dx;
% 
% plot(x,Ugraph(1,:));
% hold on
% fplot(uexact)
% mesh(t,x,Ugraph)
% axis([0,L,0,T])
% grid on


%% Part b

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
fx = 0;
w = 0.1/dt;
%Fxt = (w*cos(w*t)+D*k^2*sin(w*t))*cos(k*x);
lam = D*dt/(dx)^2;


%boundary conditions
U0(1) = sin(w);
U0(n+1) = sin(w)*cos(k);
Ugraph = U0;
Un(1)=U0(1);
Un(n+1)=U0(n+1);

%
%[ a -b  0  0 
% -c  a -b  0 
%  0 -c  a -b
%  0  0 -c  a ] 

%time loop 
for t=2:n+1
    for i=1:n-1 
        if i == 1
            f(i)=lam*U0(i)+2*(1-lam)*U0(i+1)+lam*U0(i+2)+lam*Un(i);
        elseif i == n-1
            f(i)=lam*U0(i)+2*(1-lam)*U0(i+1)+lam*U0(i+2)+lam*Un(n);
        else
            f(i)=lam*U0(i)+2*(1-lam)*U0(i+1)+lam*U0(i+2);
        end
    end
    
    %set up coefficients
    b = lam*ones(n-2,1);
    c = b;
    a = (2*(1+lam))*ones(n-1,1);
    %put into matrix
    matrix = diag(a)+ diag(-b,1)+ diag(-c,-1);
    
    %divide by function values
    Ufinal = matrix\f';
    %make vector for new row
    Un=[Un(1),Ufinal',Un(n+1)];
    %Add to graphing value
    Ugraph(t,:)=Un;
    
    %reset U0 with new value
    U0 = Un;
    U0(1) = sin(w*t);
    U0(n+1) = sin(w*t)*cos(k*t);
    Un = U0;
end

for t = 1:n+1
    for x = 1:n+1
        uexactFunction(t,x) = sin(w*(t-1))*cos(k*(x-1)*dx);
    end
end
