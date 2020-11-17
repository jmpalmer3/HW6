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
L = pi;
D = 0.1;
T = 10;
F = 0;
g0 = 0;
gL = g0;
fx = sin(k*x);
uexact(x,t) = exp(-D*k^2*t)*sin(k*x);