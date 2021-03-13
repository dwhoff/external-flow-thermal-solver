clear,clc
close all

%%%%%%%%%%%%%%%%%% - External Flow Thermal Solver - %%%%%%%%%%%%%%%%%%%%%%%
% Description:
% This code solves heat transfer problems in external flows. A description 
% of the main settings is given below.
%
% fluid: defines the fluid object for look-up in Cantera
%        options: 'air', 'water', ..., 'user defined'
%        note: if 'user defined', fluid properties must be input manually.
%              This option must be selected if Cantera is not installed.
% BC: defines the type of boundary condition used for the internal fluid
%        options: 'prescribed temp', 'conjugate'
%        note: if 'prescribed temp', the wall temperature must be input
%              manually
% wall_conduction: defines the discretization used in the solid domain
%        options: 'thin', 'thick'
%        note: BC must be set to 'conjugate'
% transient: toggles steady or unsteady analysis
%        options: true, false
%        note: BC must be set to 'conjugate'
% radiation: toggles radiation effects on the outer pipe wall
%        options: true, false
%        note: BC must be set to 'conjugate'
% write: toggles option to write data to text file
%        options: true, false
%
% Based on Stanford's course, ME352C: Convective Heat Transfer
% Written by Davis Hoffman, 2021

% Settings
fluid1 = 'water';
fluid2 = 'water';
BC = 'conjugate';
wall_conduction = 'thick';
transient = false;
radiation = true;
write = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS - Parameters
t = 10e-3; % plate thickness, m
L = 0.2; % plate length, m
P1 = 101e3; % freestream 1 pressure, Pa
T_inf1 = 20; % freestream 1 temperature, deg C
u_inf1 = 1; % freestream 1 velocity, m/s
x_stag1 = L/2; % freestream 1 stagnation line location, m
m1 = 1; % freestream 1 powerlaw acceleration parameter
C1 = u_inf1; % freestream 1 powerlaw constant, u_inf = C*x^m
P2 = 101e3; % freestream 2 pressure, Pa
T_inf2 = 30; % freestream 2 temperature, deg C
u_inf2 = 1; % freestream 2 velocity, m/s
x_stag2 = 0; % freestream 2 stagnation line location, m
m2 = 0; % freestream 2 powerlaw acceleration parameter
C2 = u_inf2; % freestream 2 powerlaw constant, u_inf = C*x^m
rho_p = 5000; % plate wall density, kg/m^3
kp = 30; % plate wall thermal conductivity, W/m.K
cpp = 4000; % plate wall specific heat capacity, J/kg.K
ho1 = 0; % surface 1 imposed convection coefficient, W/m^2.K
ho2 = 0; % surface 2 imposed convection coefficient, W/m^2.K
sigma = 5.67e-8; % Stefan-Boltzmann constant, W/m^2.K^4
eps1 = 0; % surface 1 emissivity
eps2 = 0; % surface 2 emissivity
Trad1 = 0; % surface 1 radiation target temperature, deg C
Trad2 = 0; % surface 2 radiation target temperature, deg C
Re_trans = 2e5; % transition Reynolds number
N = 4; % # of elements in wall-normal direction (if wall-conduction='thick')
M = 32; % # of elements in streamwise direction
To_i = (T_inf1+T_inf2)/2; % initial plate wall temperature (if transient=true), deg C
dt = 1; % time step (if transient=true), s
time_f = 200; % final time (if transient=true), s
tol = 1e-3; % wall temperature tolerance, K
writename = 'output.txt';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert temperatures to Kelvin scale
To_i = To_i+273.15;
T_inf1 = T_inf1+273.15;
T_inf2 = T_inf2+273.15;
Trad1 = Trad1+273.15;
Trad2 = Trad2+273.15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize fluid properties at freestream temperatures
if strcmp(fluid1,'water')
    fl_obj1 = Water;
elseif strcmp(fluid1,'air')
    fl_obj1 = Air;
else
    fl_obj1 = GRI30('Multi'); % open Cantera ideal gas object
end
if strcmp(fluid1,'user defined')
    disp('Enter fluid properties: ')
    rho1 = input(sprintf('\tDensity (kg/m^3): '));
    nu1 = input(sprintf('\tKinematic viscosity (m^2/s): '));
    cp1 = input(sprintf('\tSpecific heat (J/kg.K): '));
    k1 = input(sprintf('\tThermal conductivity (W/m.K): '));
else
    [rho1,nu1,cp1,k1] = get_fluid_props(fl_obj1,fluid1,P1,T_inf1);
end
alpha1 = k1/rho1/cp1; % thermal diffusivity, m^2/s
Pr1 = nu1/alpha1;
Re1 = u_inf1*L/nu1;

if strcmp(fluid2,'water')
    fl_obj2 = Water;
elseif strcmp(fluid2,'air')
    fl_obj2 = Air;
else
    fl_obj2 = GRI30('Multi'); % open Cantera ideal gas object
end
if strcmp(fluid2,'user defined')
    disp('Enter fluid properties: ')
    rho2 = input(sprintf('\tDensity (kg/m^3): '));
    nu2 = input(sprintf('\tKinematic viscosity (m^2/s): '));
    cp2 = input(sprintf('\tSpecific heat (J/kg.K): '));
    k2 = input(sprintf('\tThermal conductivity (W/m.K): '));
else
    [rho2,nu2,cp2,k2] = get_fluid_props(fl_obj2,fluid2,P2,T_inf2);
end
alpha2 = k2/rho2/cp2; % thermal diffusivity, m^2/s
Pr2 = nu2/alpha2;
Re2 = u_inf2*L/nu2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize variables
dx = L/M; % pipe discretization length, m
if strcmp(wall_conduction,'thin')
    N = 1;
end
dy = t/N;
xc = dx/2:dx:L-dx/2; % cell center coordinates
yc = dy/2:dy:t-dy/2;
xe = 0:dx:L; % cell face coordinates
ye = 0:dy:t;
[X,Y] = meshgrid(xc,yc);

% round stagnation locations to nearest cell face
ind_stag1 = round(x_stag1/dx);
x_stag1 = ind_stag1*dx;
ind_stag2 = round(x_stag2/dx);
x_stag2 = ind_stag2*dx;

Acond = dy; % conduction area, m^2
Aconv = dx; % convection area, m^2
mass = rho_p*dx*dy; % element mass, m^2

Bx = kp*Acond/dx; % streamwise conductance parameter, W/K
By = kp*Aconv./dy; % wall-normal conductance parameter, W/K
if transient
    Nt = ceil(time_f/dt); % number of time steps
else
    Nt = 1;
end

% form BC vector (external heat transfer on elements, in Watts)
q_dot = zeros(M,1); % INPUT (when BC='conjugate')
q_dot(end/2:end) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create block diagonal conduction matrix
A_cond = zeros(N*M);
for i = 1:N
    tmpx = Bx.*(diag(2*ones(M,1))-diag(ones(M-1,1),1)-diag(ones(M-1,1),-1));
    tmpx(1,1) = tmpx(1,1)/2;
    tmpx(end,end) = tmpx(end,end)/2;
    if (i == 1 || i==N) && strcmp(wall_conduction,'thick')
        tmpy = By*diag(ones(M,1));
    elseif i>1 && i<N && strcmp(wall_conduction,'thick')
        tmpy = 2*By*diag(ones(M,1));
    else
        tmpy = zeros(M);
    end
    A_cond((i-1)*M+1:i*M,(i-1)*M+1:i*M) = tmpx+tmpy;
end
tmp = repmat(By,M*(N-1),1);
tmp = tmp(:);
tmpy = -diag(tmp,M)-diag(tmp,-M);
A_cond = A_cond+tmpy;

% create wall temperature vector, in K
if strcmp(BC,'prescribed temp')
    T = 10+10*sin(2*pi*xc'/L)+T_in; % INPUT (when BC='prescribed temp')
elseif strcmp(BC,'conjugate')
    T = To_i*ones(N*M,1);
end

clear tmpx tmpy Bx By

[G_star_lam1,b1_1,b2_1] = build_base_DGF_laminar(M,m1);
[G_star_lam2,b1_2,b2_2] = build_base_DGF_laminar(M,m2);
Cw1 = @(k,Pr,nu) 0.561*(2*m1/(m1+1)+0.202)^0.115*sqrt((m1+1)/2)*k*Pr^b1_1* ...
            sqrt(C1*dx^(m1+1)/nu);
Cw2 = @(k,Pr,nu) 0.561*(2*m2/(m2+1)+0.202)^0.115*sqrt((m2+1)/2)*k*Pr^b1_2* ...
    sqrt(C2*dx^(m2+1)/nu);
if m1 == 0
    [G_star_turb1] = build_base_DGF_turbulent(M,m1);
    Ct1 = @(k,Pr,nu) 0.0287*k*Pr^0.6*(u_inf1*dx/nu)^0.8;
else
    G_star_turb1 = G_star_lam1;
    Ct1 = @(k,Pr,nu) 0.561*(2*m1/(m1+1)+0.202)^0.115*sqrt((m1+1)/2)*k*Pr^b1_1* ...
        sqrt(C1*dx^(m1+1)/nu);
end
if m2 == 0
    [G_star_turb2] = build_base_DGF_turbulent(M,m2);
    Ct2 = @(k,Pr,nu) 0.0287*k*Pr^0.6*(u_inf2*dx/nu)^0.8;
else
    G_star_turb2 = G_star_lam2;
    Ct2 = @(k,Pr,nu) 0.561*(2*m2/(m2+1)+0.202)^0.115*sqrt((m2+1)/2)*k*Pr^b1_2* ...
        sqrt(C2*dx^(m2+1)/nu);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iterate through time steps
time = 0;
figure(1); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.125, 0.25, 0.75, 0.5]);
for j = 1:Nt
    fprintf('Timestep %d of %d, t = %.2f s\n',j,Nt,time)
    
    T_prev = T;
    T_old = zeros(N*M,1);
    resid = Inf;
    converged = 0;
    iter = 0;
    % iterate solver until film temperature converges
    while ~converged    
        %%
        % form linear system for the wall temperature
        % Steady: A*T=b
        % Transient: m*cp*dT/dt=-A*T+b

        % compute transition location
        Re1 = u_inf1*L/nu1;
        Re2 = u_inf2*L/nu2;
        x_trans1 = Re_trans*nu1/u_inf1;
        ind_trans1 = round(x_trans1/dx);
        x_trans2 = Re_trans*nu2/u_inf2;
        ind_trans2 = round(x_trans2/dx);
        
        % form DGF matrix
        G_lam1 = Cw1(k1,Pr1,nu1)*G_star_lam1;
        G_lam2 = Cw2(k2,Pr2,nu2)*G_star_lam2;
        G_turb1 = Ct1(k1,Pr1,nu1)*G_star_turb1;
        G_turb2 = Ct2(k2,Pr2,nu2)*G_star_turb2;
        G1 = configure_DGF(G_lam1,G_turb1,ind_stag1,ind_trans1);
        G2 = configure_DGF(G_lam2,G_turb2,ind_stag2,ind_trans2);
        
        G1_tmp = zeros(N*M);
        G1_tmp(1:M,1:M) = G1;
        G2_tmp = zeros(N*M);
        G2_tmp(end-M+1:end,end-M+1:end) = G2;
        G = G1_tmp+G2_tmp;
        
        % create external heat transfer matrix/vector
        A_convo = zeros(N*M);
        A_convo(1:M,1:M) = diag(ho1*Aconv*ones(M,1));
        A_convo(end-M+1:end,end-M+1:end) = diag(ho2*Aconv*ones(M,1));
        b_convo = zeros(N*M,1);
        b_convo(1:M) = ho1*Aconv*T_inf1;
        b_convo(end-M+1:end) = ho2*Aconv*T_inf2;

        % create radiation heat transfer matrix/vector
        if radiation
            Rrad1 = 1./(eps1*sigma*Aconv*(T(1:M)+Trad1).*(T(1:M).^2+Trad1^2));
            Rrad2 = 1./(eps2*sigma*Aconv*(T(end-M+1:end)+Trad2).*(T(end-M+1:end).^2+Trad2^2));
            A_rad = zeros(N*M);
            A_rad(1:M,1:M) = diag(1./Rrad1);
            A_rad(end-M+1:end,end-M+1:end) = diag(1./Rrad2);
            b_rad = zeros(N*M,1);
            b_rad(1:M) = Trad1./Rrad1;
            b_rad(end-M+1:end) = Trad2./Rrad2;
        else
            A_rad = zeros(N*M);
            b_rad = zeros(N*M,1);
        end

        % create vector of external heat transfers
        b_ext = zeros(N*M,1);
        b_ext(1:M) = q_dot;

        % create vector that absorbs the inlet temperature terms
        b_inlet1 = zeros(N*M,1);
        b_inlet1(1:M) = T_inf1*G1(1:M,1:M)*ones(M,1);
        b_inlet2 = zeros(N*M,1);
        b_inlet2(end-M+1:end) = T_inf2*G2(end-M+1:end,end-M+1:end)*ones(M,1);
        b_inlet = b_inlet1+b_inlet2;

        if transient
            f = @(T) 1./(mass*cpp).*(-(A_cond+G+A_convo+A_rad)*T+ ...
                (b_ext+b_convo+b_rad+b_inlet));
            T = advance_RK4(f,T_prev,dt);
        else
            if strcmp(BC,'conjugate')
                T = (A_cond+G+A_convo+A_rad)\(b_ext+b_convo+b_rad+b_inlet);
            end
        end
        
        resid = sqrt(mean((T-T_old).^2));
        if resid<tol
            converged = 1;
        else
            T_old = T;
        end
        
        % get convective heat flux into fluid
        qconv1 = G(1:M,1:M)*(T(1:M)-T_inf1)+ho1*Aconv*(T(1:M)-T_inf1);
        qconv2 = G(end-M+1:end,end-M+1:end)*(T(end-M+1:end)-T_inf2)+ ...
            ho2*Aconv*(T(end-M+1:end)-T_inf2);

        % get fluid properties at the film temperature
        Tfilm1 = (mean(T(1:M))+T_inf1)/2;
        if ~strcmp(fluid1,'user defined')
            [rho1,nu1,cp1,k1] = get_fluid_props(fl_obj1,fluid1,P1,Tfilm1);
            alpha1 = k1/rho1/cp1; % thermal diffusivity, m^2/s
        end
        Tfilm2 = (mean(T(end-M+1:end))+T_inf2)/2;
        if ~strcmp(fluid2,'user defined')
            [rho2,nu2,cp2,k2] = get_fluid_props(fl_obj2,fluid2,P2,Tfilm2);
            alpha2 = k2/rho2/cp2; % thermal diffusivity, m^2/s
        end
        
        iter = iter+1;

        if iter==1
            fprintf('\tIteration: %d\n',iter)
        else
            fprintf('\tIteration: %d, Residual: %.2e\n',iter,resid)
        end
    end
    time = time+dt;
    
    hh = figure(1);
    if  strcmp(wall_conduction,'thin')
        plot(xc*1e3,T,'k-')
        xlabel('$x$ (mm)','interpreter','latex')
        ylabel('$T_o(x)$ ($^o$C)','interpreter','latex')
    elseif  strcmp(wall_conduction,'thick')
        pcolor(X*1e3,Y*1e3,reshape(T,M,N)'); shading interp; colormap hot; colorbar; %caxis([300 375])
        xlabel('$x$ (mm)','interpreter','latex')
        ylabel('$y$ (mm)','interpreter','latex')
    end
    title(sprintf('($t=%d$ s) Wall temperature distribution (K)',time),'interpreter','latex')
    plotfixer
    drawnow
end
To1 = T(1:M);
To2 = T(end-M+1:end);

clear T_old T_prev tmp A_cond A_conv1 A_conv2 A_convo A_rad G G1 G1_tmp ...
    G2 G2_tmp G_lam1 G_lam2 G_star_lam1 G_star_lam2 G_star_turb1 G_star_turb2 ...
    G_turb1 G_turb2 b_convo b_ext b_inlet b_rad q_dot Rrad1 Rrad2  ...
    Ct1 Ct2 Cw1 Cw2 fl_obj1 fl_obj2 j

disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')

%%
% check errors
hm1 = 0.664*k1*Pr1^1/3*sqrt(u_inf1/nu1/L);
hm2 = 0.664*k2*Pr2^1/3*sqrt(u_inf2/nu2/L);
Bi = max([hm1 hm2])*t/kp; % Biot number
if Bi>0.1 && strcmp(wall_conduction,'thin') && strcmp(BC,'conjugate')
    disp('Warning: Bi exceeds 0.1. Consider a thick wall analysis for improved results.')
end
tau_solid = rho_p*cpp*t/(2*max([hm1 hm2])); % solid time response, s
tau_adv = L/min([u_inf1 u_inf2]); % advection (flow-through) time scale
if tau_solid<tau_adv && transient
    disp('Warning: Quasi-steady assumption may be invalid.')
end

%%
% plot
if strcmp(wall_conduction,'thin')
    figure(1)
    plot(xc*1e3,To1-273.15,'k-'); hold on
    xlabel('$x$ (mm)','interpreter','latex')
    ylabel('$T_o(x)$ ($^o$C)','interpreter','latex')
    scale = axis;
    axis([scale(1) 1e3*L scale(3) scale(4)])
elseif strcmp(wall_conduction,'thick')
    figure(1)
    TT = reshape(T,M,N)';
    pcolor(X*1e3,Y*1e3,TT); shading interp; colormap hot; colorbar
    xlabel('$x$ (mm)','interpreter','latex')
    ylabel('$y$ (mm)','interpreter','latex')
    title('Wall temperature distribution (K)')

    figure(2)
    plot(xc*1e3,To1-273.15,'k--'); hold on
    plot(xc*1e3,To2-273.15,'k-.');
    xlabel('$x$ (mm)','interpreter','latex')
    ylabel('$T_o(x)$ ($^o$C)','interpreter','latex')
    legend({'$T_{o,1}(x)$','$T_{o,2}(x)$'},'interpreter','latex')
    legend boxoff; scale = axis;
    axis([scale(1) 1e3*L scale(3) scale(4)])
end

figure(3)
plot(xc*1e3,qconv1/Aconv,'k--'); hold on
plot(xc*1e3,qconv2/Aconv,'k-.');
xlabel('$x$ (mm)','interpreter','latex')
ylabel('$\dot{q}''''_o(x)$ (W/m$^2$)','interpreter','latex')
legend({'$\dot{q}''''_{o,1}(x)$','$\dot{q}''''_{o,2}(x)$'},'interpreter','latex')
legend boxoff; scale = axis;
axis([scale(1) 1e3*L scale(3) scale(4)])

figure(4)
plot(xc*1e3,qconv1/Aconv./(To1-T_inf1),'k--'); hold on
plot(xc*1e3,qconv2/Aconv./(To2-T_inf2),'k-.');
xlabel('$x$ (mm)','interpreter','latex')
ylabel('$h$ (W/m$^2$K)','interpreter','latex')
legend({'$h_1(x)$','$h_2(x)$'},'interpreter','latex')
legend boxoff; scale = axis;
axis([scale(1) 1e3*L scale(3) scale(4)])

plotfixer

%%
% write to file
if write
    To1_tmp = To1;
    To2_tmp = To2;
    To1 = zeros(M+1,1);
    To1(2:end-1) = (To1_tmp(1:end-1)+To1_tmp(2:end))/2;
    To1(1) = 1.5*To1_tmp(1)-0.5*To1_tmp(2);
    To1(end) = 1.5*To1_tmp(end)-0.5*To1_tmp(end-1);
    To2 = zeros(M+1,1);
    To2(2:end-1) = (To2_tmp(1:end-1)+To2_tmp(2:end))/2;
    To2(1) = 1.5*To2_tmp(1)-0.5*To2_tmp(2);
    To2(end) = 1.5*To2_tmp(end)-0.5*To2_tmp(end-1);
    qconv1_tmp = qconv1;
    qconv2_tmp = qconv2;
    qconv1 = zeros(M+1,1);
    qconv2 = zeros(M+1,1);
    qconv1(2:end-1) = (qconv1_tmp(1:end-1)+qconv1_tmp(2:end))/2;
    qconv1(1) = 1.5*qconv1_tmp(1)-0.5*qconv1_tmp(2);
    qconv1(end) = 1.5*qconv1_tmp(end)-0.5*qconv1_tmp(end-1);
    qconv2(2:end-1) = (qconv2_tmp(1:end-1)+qconv2_tmp(2:end))/2;
    qconv2(1) = 1.5*qconv2_tmp(1)-0.5*qconv2_tmp(2);
    qconv2(end) = 1.5*qconv2_tmp(end)-0.5*qconv2_tmp(end-1);
    fileID = fopen(writename,'w');
    fprintf(fileID,'x (m)\t\tTo1 (K)\t\tTo2 (K)\t\tq1 (W/m^2)\tq2 (W/m^2)\n');
    for i = 1:M+1
        fprintf(fileID,'%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n',xe(i),To1(i),To2(i), ...
            qconv1(i)/Aconv,qconv2(i)/Aconv);
    end
    fclose(fileID);
end
clear To_tmp qconv_tmp fileID