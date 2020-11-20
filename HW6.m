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
    for j=1:n-1 
        if j == 1
            f(j)=lam*U0(j)+2*(1-lam)*U0(j+1)+lam*U0(j+2)+lam*Un(j);
        elseif j == n-1
            f(j)=lam*U0(j)+2*(1-lam)*U0(j+1)+lam*U0(j+2)+lam*Un(n);
        else
            f(j)=lam*U0(j)+2*(1-lam)*U0(j+1)+lam*U0(j+2);
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

%uexact = @(x) exp(-D*k^2*1)*sin(k*x);
x = (0:1:10)*dx;

plot(x,Ugraph(5,:));
hold on
plot(x,uexactFunction(5,:));
axis([0,L,0,1])
grid on
xlabel('Length [x]')
ylabel('Time [t]')
title('Crank-Nicolson at time = 5')

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
w = 1/dt;
%Fxt = @(x) (w*cos(w*t)+D*k^2*sin(w*t))*cos(k*x);
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
    
    %set up coefficients
    b = lam*ones(n-2,1);
    c = b;
    a = (2*(1+lam))*ones(n-1,1);
    %put into matrix
    matrix = diag(a)+ diag(-b,1)+ diag(-c,-1);
    
    for j=1:n-1 
        if j == 1
            f(j)=lam*U0(j)+2*(1-lam)*U0(j+1)+lam*U0(j+2)+lam*Un(j);
        elseif j == n-1
            f(j)=lam*U0(j)+2*(1-lam)*U0(j+1)+lam*U0(j+2)+lam*Un(n);
        else
            f(j)=lam*U0(j)+2*(1-lam)*U0(j+1)+lam*U0(j+2);
        end
    end
    
    %divide by function values
    Ufinal = matrix\f';
    
    for j=1:n-1 
        Ufinal(j) = Ufinal(j) + (w*cos(w*t)+D*k^2*sin(w*t))*cos(k*j*dx);
    end
    %make vector for new row
    Un=[Un(1),Ufinal',Un(n+1)];
    %Add to graphing value
    Ugraph(t,:)=Un;
    
    %reset U0 with new value
    U0 = Un;
    U0(1) = sin(w*t);
    U0(n+1) = sin(w*t)*cos(k*L*dx);
    Un = U0;
end

Ugraph(1,1) = 0;
Ugraph(1,n+1) = 0;

for t = 1:n+1
    for x = 1:n+1
        uexactFunction(t,x) = sin(w*(t-1))*cos(k*(x-1)*dx);
    end
end

e = 0;
for t = 1:n+1
    for x = 1:n+1
        e = (Ugraph(t,x)-uexactFunction(t,x))/uexactFunction(t,x);
    end
end
e = 1/n*e
