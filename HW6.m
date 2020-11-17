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
%uexact(x,t) = exp(-D*k^2*t)*sin(k*x);

%boundary conditions
U0(1) = g0;

for x = 2:n-1
    U0(x) = sin(k*x);
end
    U0(n) = gL;

Un(1)=g0; 
Un(n)=gL;
