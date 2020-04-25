clear all;
close all;
clc;

%% initialize   v*U'' + betta*U' + gamma*U = exp(-2x)
Domain = [0, 1];
Coefficients = [-1 2 1]; % v betta gamma
Number_of_Points = 200;

%% define any function you want

Func = @(y) exp(-2*y);

%% get Orders


Conditions = [];
Conditions(1) = input('Do you want use Dirichlet or Numann for BC0 ?[D / N] \n', 's');
Conditions(2) = input('Enter First Boundry Condition \n');
Conditions(3) = input('Do you want use Dirichlet or Numann for BC1 ?[D / N] \n', 's');
Conditions(4) = input('Enter Second Boundry Condition \n');
Conditions(5) = input('Do you want to use Upwind or Centered?[U / C]\n', 's');



[PDE, F, domain, Method0, Method1, DiscMethod] = ClearOrder(Conditions, Func, Domain, Coefficients, Number_of_Points);

U = PDE\F;


%% Plotting

figure(1);

plot(domain, U)
title(['Solving Eq. BC0 => ', Method0, ' and BC1 => ', Method1, ' by ', DiscMethod,' method']);
xlabel('x axis');
ylabel('y axis');

