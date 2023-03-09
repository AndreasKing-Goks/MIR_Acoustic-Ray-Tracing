clear all
clc 
clf
%% INITIAL CONDITION
% Environment
Zseabed = 3500;
Zinv = 700;
gamma_0 = .0163;
gamma_1 = -0.026;

X0 = 0;
Z0 = 600;
C0 = 1450;
theta0 = 4*pi/180;

Cinv = C0 + gamma_1*Zinv;
%% ITERATIVE CALCULATION
dt = 0.01; % timestep size
ts = 1500; % timestep
time = dt *ts;

% Beam Defenition
num_rays = 20;
thetas_beam = linspace(-theta0/2, theta0/2, num_rays);

% Looping for each beam
for j= 1:num_rays

    % Create containers
    X = zeros(1,ts);
    Z = zeros(1,ts);
    C = zeros(1,ts);
    theta = zeros(1,ts);

    % Intial Values at t=0
    X(1) = X0;
    Z(1) = Z0;
    if Z0 >= Zinv
        gamma = gamma_0;
    else
        gamma = gamma_1;
    end
    C(1) = C0 + gamma*Z(1);
    theta(1) = thetas_beam(j);

    % Looping Process
    for i = 1:ts-1
        X(i+1) = X(i) + C(i)*dt*cos(theta(i));
        Z(i+1) = Z(i) + C(i)*dt*sin(theta(i));
        if Z(i+1) <= Zinv
            C(i+1) = C0 + gamma_1*Z(i+1);
        elseif Z(i+1) > Zinv
            C(i+1) = Cinv + gamma_0*(Z(i+1)-Zinv);
        end 
        K= cos(theta(i))*C(i+1)/C(i);
        if K>=1 || Z(i+1) < 0 || Z(i+1) >= Zseabed
            theta(i+1) = -theta(i-1);
            Z(i+1) = Z(i);
            C(i+1) = C(i-1);
        else
            theta(i+1) = acos(K) * sign(theta(i));
        end
        if Z(i+1) >= Zinv
            gamma = gamma_0;
        else
            gamma = gamma_1;
        end
    end

    plot(X,-Z,'b')
    hold on
    grid on
    xlabel('x')
    ylabel('sea depth')
    title('Accoustic Tracing in Mediterranean Sea')
end