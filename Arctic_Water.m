clear all
clc 
clf
%% INITIAL CONDITION
% Environment
Zseabed = 3500;
gamma = .0163;

X0 = 0;
Z0 = 600;
C0 = 1450;
theta0 = 4*pi/180;

%% ITERATIVE CALCULATION
dt = 0.01; % timestep size
ts = 3000; % timestep
time = dt *ts;

% Beam Defenition
num_rays = 20;
thetas_beam = linspace(-theta0/2, theta0/2, num_rays);

% Looping for each beam
for j = 1:num_rays

    % Create Containers
    X = zeros(1,ts);
    Z = zeros(1,ts);
    C = zeros(1,ts);
    theta = zeros(1,ts);

    % Initial value at t = 0
    X(1) = X0;
    Z(1) = Z0;
    C(1) = C0 + gamma*Z(1);
    theta(1) = thetas_beam(j);

    % Looping Process
    for i = 1:ts-1
        X(i+1) = X(i) + C(i)*dt*cos(theta(i));
        Z(i+1) = Z(i) + C(i)*dt*sin(theta(i)); 
        C(i+1) = C0 + gamma*Z(i+1); 

        K= cos(theta(i))*C(i+1)/C(i);

        if K>=1 || Z(i+1) < 0 || Z(i+1) >= Zseabed
            theta(i+1) = -theta(i-1);
            Z(i+1) = Z(i);
            C(i+1) = C(i-1);
        else
            theta(i+1) = acos(K) * sign(theta(i));
        end
    end

    % Plot
    plot(X,-Z, 'b')
    xlabel('x')
    ylabel('sea depth')
    title('Accoustic Tracing in Arctic Ocean')
    hold on
    grid on
end
