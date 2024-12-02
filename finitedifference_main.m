%{
Eden Zafran
5/25/22
DIY CFD - Main Script
NACA 0012 Airfoil AOA = 0deg
%}

clc; clear all; close all;
%%
% run panel method w/ n = 400
finitedifference_script2; 
clc; clear all; close all;
%% givens
vinf = 1;
c = 1;              % chord length
L = 5*c + c + 5*c;  % length of computational domain
H = 4*c;            % height of computational domain
aoa = 0;            % deg
aoa = aoa*pi/180;   % rad
t = 0.12;           % thicnkess of NACA 0012 airfoil
% grid size
N = 220;
deltaX = L/N;
M = 40;
deltaY = H/M;
epsilon = 1e-6;
error = 1;

%% generate grid
xgrid = zeros(N+1, M+1);
ygrid = zeros(N+1, M+1);

for n = 1:N+1
    for m = 1: M+1
        xgrid(n,m) = (n-1)*deltaX;
        ygrid(n,m) = (m-1)*deltaY;
    end
end

%% generate airfoil coordinates
airfoilcoords = zeros(N+1, 1);
for n = 1:N+1
    if ((xgrid(n,1) >= 5*c) && (xgrid(n,1)<= 6*c))
        xlocation = xgrid(n,1) - 5*c;   % get x location relative to leading edge
        airfoilcoords(n) = ((t*c)/0.2)*( 0.2969*sqrt(xlocation/c) - 0.1260*(xlocation/c) - 0.3516*(xlocation/c)^2 + 0.2843*(xlocation/c)^3 - 0.1036*(xlocation/c)^4 );
    end
end

%% initial guess
psi = zeros(N+1, M+1);
psiOld = zeros(N+1, M+1);

% immediately apply BCs
[psi] = applyBCs(psi, airfoilcoords, H, deltaY, vinf, N, M);

%% apply Jacobi method
while error > epsilon
    for n = 2:N
        for m =2:M
            psi(n,m) = ((deltaY^2/(2*(deltaX^2+deltaY^2)))*(psi(n+1,m) + psi(n-1,m))) + ...
                       ((deltaX^2/(2*(deltaX^2+deltaY^2)))*(psi(n,m+1) + psi(n,m-1)));
        end
    end

    % apply BCs
    [psi] = applyBCs(psi, airfoilcoords, H, deltaY, vinf, N, M);

    % calculate error
    error = max(max(psi-psiOld));

    % update
    psiOld = psi;
end

%% calculate pressure coeff (cp)
cp = zeros(N+1, M+1);
%{
for n = 1:N+1
    for m = 1:M+1
        if (m==1 && n==1)   % bottom left corner
            u = (psi(n,2) - psi(n,1))/deltaY; % forward approx
            v = -(psi(2,m) - psi(1,m))/deltaX;   % forward approx
        elseif (m==1 && n==N+1) % bottom right corner
            u = (psi(n,2) - psi(n,1))/deltaY; % forward approx
            v = -(psi(N+1,m) - psi(N,m))/deltaX;   % backward approx
        elseif (m==M+1 && n==1)   % top left corner
            u = (psi(n,M+1) - psi(n,M))/deltaY; % backward approx
            v = -(psi(2,m) - psi(1,m))/deltaX;   % forward approx
        elseif (m==M+1 && n ==N+1)  % top right corner
            u = (psi(n,M+1) - psi(n,M))/deltaY; % backward approx
            v = -(psi(N+1,m) - psi(N,m))/deltaX;   % backward approx
        elseif (m == 1)     % bottom
            u = (psi(n,2) - psi(n,1))/deltaX; % forward approx
            v = -(psi(n+1,m) - psi(n-1,m))/(2*deltaX); % central approx - unchanged
        elseif (m == M+1)   % top
            u = (psi(n,M+1) - psi(n,M))/deltaX; % backward approx
            v = -(psi(n+1,m) - psi(n-1,m))/(2*deltaX); % central approx - unchanged
        elseif (n == 1)     % right
            u = (psi(n,m+1) - psi(n,m-1))/(2*deltaY);   % central approx - unchanged
            v = -(psi(1,m) - psi(2,m))/(2*deltaX);   % forward approx
        elseif (n == N+1)   % left
            u = (psi(n,m+1) - psi(n,m-1))/(2*deltaY);   % central approx - unchanged
            v = -(psi(1,M+1) - psi(2,M))/(2*deltaX);   % backward approx
        else                % interior points
            u = (psi(n,m+1) - psi(n,m-1))/(2*deltaY);   % central approx
            v = -(psi(n+1,m) - psi(n-1,m))/(2*deltaX);   % central approx
        end

        V = sqrt(u^2 + v^2);
        cp(n,m) = 1 - (V/vinf)^2;
    end
end
%}

for n=1:N
    for m =1:M
        u = (psi(n,m+1) - psi(n,m))/(deltaY);   % forward approx
        v = -(psi(n+1,m) - psi(n,m))/(deltaX);   % forward approx
        V = sqrt(u^2 + v^2);
        cp(n,m) = 1 - (V/vinf)^2;
    end
end

%% plots
% 2.1 - plot mesh
Z = zeros(size(xgrid));
mesh(xgrid,ygrid,Z,'EdgeColor',[0,0,0])
view(2)
xlim([0,L])
ylim([0,H])
title("Grid")
xlabel("x-axis")
ylabel("y-axis")

% 2.2 - distribution of cp vs x along airfoil surface
% get only x coordinates of airfoil
cpxaxis = xgrid(:,1);
cpxaxis = cpxaxis(101:121);
% get cp along airfoil surface
cpalongairfoil = cp(:,1);
cpalongairfoil = cpalongairfoil(101:121);

figure
plot(cpxaxis, cpalongairfoil)
title("Cp as a function of x across airfoil surface")
set(gca,'Ydir','reverse')
xlabel("x - Airfoil Surface")
ylabel("Cp")

%%
%%% import files
cp_gregory = importdata("Cpsurf_GregoryOReilly.dat");
cp_gregory_x = cp_gregory(:,1) + 5;  % add 5 to shift to match LE/TE of this airfoil in middle of domain
cp_gregory_cp = cp_gregory(:,2);

cp_mason = importdata("Cpsurf_Mason.dat");
cp_mason_x = cp_mason(:,1) + 5;     % add 5 to shift to match LE/TE of this airfoil in middle of domain
cp_mason_cp = cp_mason(:,2);

cp_400panels = importdata("400panelsmethodcp.txt");
cp_400panels_x = cp_400panels(:,1) + 5; % add 5 to shift to match LE/TE of this airfoil in middle of domain
cp_400panels_cp = cp_400panels(:,2);

figure
plot(cpxaxis, cpalongairfoil)
hold on
plot(cp_400panels_x, cp_400panels_cp, '-.')
title("Figure 1: CFD vs My Panel Method (n=400)")
set(gca,'Ydir','reverse')
legend("CFD", "My Panel Method")
xlim([5 6])
xlabel("x - Airfoil Surface")
ylabel("Cp")

% experimental results - GregoryOReilly
figure
plot(cpxaxis, cpalongairfoil)
hold on
plot(cp_gregory_x, cp_gregory_cp, '-.')
title("Figure 2: CFD vs Experimental Results")
set(gca,'Ydir','reverse')
legend("CFD", "Experimental Results")
xlim([5 6])
xlabel("x - Airfoil Surface")
ylabel("Cp")

% published panel method results - Mason
figure
plot(cpxaxis, cpalongairfoil)
hold on
plot(cp_mason_x, cp_mason_cp, '-.')
title("Figure 3: CFD vs Published Panel Method Results")
set(gca,'Ydir','reverse')
legend("CFD", "Published Panel Method")
xlim([5 6])
xlabel("x - Airfoil Surface")
ylabel("Cp")

%%
% 2.3 - contour plot for pressure coeff
figure
contourf(xgrid, ygrid, cp, 'LevelStep', 0.05)
title("2D Contours of Cp (Entire Domain)")

% 2.4 contour plot for streamlines
figure
contour(xgrid, ygrid, psi, 'LevelStep', 0.1)
title("Streamlines (Entire Domain)")
xlabel('L')
ylabel('H')

%% function definitions

function [psi] = applyBCs(psi, airfoilcoords, H, deltaY, vinf, N, M)
    % bottom BC
    for n=1:N+1
        psi(n,1) = -psi(n,2)*(airfoilcoords(n)/(deltaY-airfoilcoords(n)));
    end
    
    % top BC
    psi(:,M+1) = vinf*H;
    
    % left BC
    for m=1:M+1
        psi(1,m) = vinf*(m-1)*deltaY;
    end
    
    % right BC
    for m=1:M+1
        psi(N+1,m) = psi(N,m);
    end
end