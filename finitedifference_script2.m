%{
Eden Zafran
MAE150B HW7
5/25/22
DIY CFD - HW3 Panel Method w/400 Panels Helper Script
NACA 0012 Airfoil AOA = 0deg
%}

clc; clear all; close all;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 1: AOA = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% constants
vinf = 1;   % free stream velocity
c = 1;      % chord length
t = 0.12;   % percentage of thickness to chord of airfoil (for NACA0012, t = 0.12)
aoa = 0;
% fprintf('AOA = %f%c \n', aoa, char(176))
aoa = aoa*(pi/180); %change aoa to rad from deg


n = 400;

%%% generate geometry - XY coordinates
X = zeros(n+1, 1);
Y = zeros(n+1, 1);

% upper surface
for j = 1:(n/2 + 1)
    X(j) = (j-1)*((2*c)/n);
    Y(j) = ((t*c)/0.2)*( 0.2969*sqrt(X(j)/c) - 0.1260*(X(j)/c) - 0.3516*(X(j)/c)^2 + 0.2843*(X(j)/c)^3 - 0.1036*(X(j)/c)^4 );
end

% lower surface
for j = (n/2 + 2): (n+1)
    X(j) = X(n+2-j);
    Y(j) = -1*Y(n+2-j);
end


%%% generate geometry - midpoint, length, phi
l = zeros(n, 1);
x = zeros(n, 1);
y = zeros(n, 1);
phi = zeros(n, 1);


for j = 1 : n

    deltaXj = X(j+1) - X(j);
    deltaYj = Y(j+1) - Y(j);

    %length
    l(j) = sqrt( deltaXj^2 + deltaYj^2 );

    % midpoint
    x(j) = (X(j) + X(j+1))/2;
    y(j) = (Y(j) + Y(j+1))/2;

    % angle ccw from pos x axis
    phi(j) = atan2(deltaYj, deltaXj);

end

%%% compute Iij ("induced velocities")
I = zeros(n,n);
J = zeros(n,n);

for i = 1:n
    for j = 1:n
        if (i == j)
            I(i,j) = pi;
            J(i,j) = 0;
        else
            % 3 step process
            xprime = (x(i) - x(j))*cos(phi(j)) + (y(i) - y(j))*sin(phi(j));
            yprime = -1*(x(i) - x(j))*sin(phi(j)) + (y(i) - y(j))*cos(phi(j));

            h = l(j)/2;
            uprime = (1/2)*log( ((xprime+h)^2 + yprime^2) /((xprime-h)^2 + yprime^2) );
            vprime = atan((xprime+h)/yprime) - atan((xprime-h)/yprime);

            deltaphi = phi(j) - phi(i);

            I(i,j) = uprime*sin(deltaphi) + vprime*cos(deltaphi);
            J(i,j) = uprime*cos(deltaphi) - vprime*sin(deltaphi);
        end
    end
end

%%% solve system
M = I;
b = zeros(n, 1);

for j = 1:n
    b(j) = -1*sin(aoa-phi(j));
end

lambda_prime = linsolve(M,b);
lambda = lambda_prime.*(2*pi*vinf);

%%% calculate Cp, capitalLambda, Cl, Cd, CmLE
% Cp
cp = zeros(n, 1);
for i = 1:n
    sumlambdaiJij = 0;
    for j = 1:n
        sumlambdaiJij = sumlambdaiJij + lambda_prime(j)*J(i,j) ;
    end
    cp(i) = 1 - ( cos(aoa-phi(i)) + sumlambdaiJij )^2;
end

% check total source strength = 0
capitalLambda = 0;
for j = 1:n
    capitalLambda = capitalLambda + lambda(j)*l(j);
end

% Ca/Cn
ca = 0;
cn = 0;
for i = 1:n
    ca = ca + (l(i)/c)*cp(i)*sin(phi(i));
    cn = cn + (l(i)/c)*cp(i)*cos(phi(i));
end
cn = -1*cn;

% Cl
cl = cn*cos(aoa) - ca*sin(aoa);

% Cd
cd = cn*sin(aoa) + ca*cos(aoa);

% CmLE
cmLE = 0;
for i = 1:n
    cmLE = cmLE + cp(i)*l(i)*( y(i)*sin(phi(i)) + x(i)*cos(phi(i)));
end
cmLE = (1/c^2)*cmLE;

% fprintf('Cl (n = %d): %d\n', n, cl);
% fprintf('Cd (n = %d): %d\n', n, cd);
% fprintf('Cm,le (n = %d): %d\n', n, cmLE);
% fprintf('Lambda (sum of source strengths) (n = %d): %d\n\n', n, capitalLambda);


% %%% import dat files
% cp_gregory = importdata("Cpsurf_GregoryOReilly.dat");
% cp_gregory_x = cp_gregory(:,1);
% cp_gregory_cp = cp_gregory(:,2);
% 
% cp_mason = importdata("Cpsurf_Mason.dat");
% cp_mason_x = cp_mason(:,1);
% cp_mason_cp = cp_mason(:,2);


%%% plots
xic = x/c;
xic = xic(1:n);
% 
% % airfoil shape
% figure
% plot(X,Y, 'b.')
% hold on
% plot(X,Y, 'b')
% xlim([0 1])
% ylim([-0.5 0.5])
% title(['1 - Airfoil Geometry (n = ', num2str(n), ')'])
% xlabel('X coordinates')
% ylabel('Y coordinates')
% 
% %lambdai vs xi/c
% figure
% plot(xic, lambda)
% title(['2 - lambda vs x/c (n = ' num2str(n), ')'])
% xlabel('X coordinates/ chord Length')
% ylabel('\lambda')
% 
% % cpi vs xi/c vs mason
% figure
% plot(xic, cp)
% hold on
% scatter(cp_mason_x, cp_mason_cp, 'r', 'x')
% title(['3 - Calculated Cp, Mason vs x/c (n = ', num2str(n), ')'])
% set(gca,'Ydir','reverse')
% xlabel('X coordinates/ chord Length')
% ylabel('Coefficient of Pressure (Cp)')
% legend('Panel Method', 'Mason')
% xlim([0 1])
% 
% % cpi vs xi/c vs gregory
% figure
% plot(xic, cp)
% hold on
% scatter(cp_gregory_x, cp_gregory_cp, 'r', 'x')
% title(['4 - Calculated Cp, Gregory OReilly vs x/c (n = ', num2str(n), ')'])
% set(gca,'Ydir','reverse')
% xlabel('X coordinates/ chord Length')
% ylabel('Coefficient of Pressure (Cp)')
% legend('Panel Method', 'Gregory OReilly')
% xlim([0 1])

cptowrite = [xic, cp];
writematrix( cptowrite, '400panelsmethodcp.txt');