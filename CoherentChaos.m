% Code for simulations
% Landau and Sompolinsky, PLOS Comp. Bio. 2018

%clear all
N=2000;         % Size of Network


seed=1;         % Seed for Random Component of Connectivity
seed2 = 1;      % Seed for initial condition of dynamics
rng(seed);      % Initialize seed of random number generator for Connectivity

%% Parameters and Model Specifications

g = 2;          % Std of Random Component
J1=0.5;         % Strength of Structured Component

Row_Balance = false;        % Set 'true' to apply row-balance constraint to Random Component
NonSymmetric_phi = false;  % Set 'true' to define non-symmetric transfer function

if NonSymmetric_phi
    
    beta=4;     % Gain parameter
    p=.5;       % Symmetry parameter
    phi = @(x) (1+exp(-beta*x)).^-p;
else
     phi = @(x) tanh(x);
end

%% Construct Connectivity Matrix

% Random Component
JJ = randn(N)/sqrt(N);

% Input Mode
Xi=ones(N,1);
% Output Mode
Nu = ones(1,N);
Nu(N/2+1:N)=-1;
 
if Row_Balance
    J=g*(JJ-JJ*Xi*Xi'/(Xi'*Xi))+J1*Xi*Nu/sqrt(N);    
else
    J=g*(JJ)+J1*Xi*Nu/sqrt(N);
end

%% Run Simulation


dt = 0.05;      % Time step for output
t_record = 300; % Time to record
t_start= 100;   % Buffer time to run before recording
t_end = t_start+t_record; 

Ts = [0,t_start:dt:t_end];  % Vector of times to output

rng(seed2)      % Initialize seed of random number generator for initial condition
H0 = randn(N,1);% Initial condition

F = @(t,h) (-h+J*phi(h)); % Define dh/dt for differential eqn solver, ode45

tic
[t,Hsol] = ode45(F,Ts,H0);
toc

%% Results

tt=t(2:end)-t_start;    % time steps corresponding to output
H = Hsol(2:end,:)';     % Rearrange output so that each row is a neuron, with recording from t_start until t_end in time steps of dt

h_bar = Xi'*H/N;        % Coherent mode of input current
delta_h = H - Xi*h_bar; % Residuals input current

phis = phi(H);          % Activity
phi_bar = Xi'*phis/N;   % Coherent mode activity

% Calculate Coherence
Coherence = sqrt(mean(h_bar.^2)/mean(mean(H.^2)));

%% Plot


tshow = tt<100; % Time to plot

figure;
plot(tt(tshow),phis(randi(N,10,1),tshow))
hold on
title(sprintf('N = %d, g = %0.1f, Seed = %d, J1 = %0.2f', N, g, seed,J1))

plot(tt(tshow),phi_bar(tshow),'k','LineWidth',2)
box off
ylim([-1,1])
