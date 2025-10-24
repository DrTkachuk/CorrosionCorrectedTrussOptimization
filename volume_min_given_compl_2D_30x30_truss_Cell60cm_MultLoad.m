% Example 3 for manuscript "A Corrosion Correction for Minimal Weight Design  
% of Truss Structures under Compliance Constraint" by Anton and Mykola
% Tkachuk for journal Technische Mechanik
% This file demonstrates design for truss structures including rust
% All trusses have width w and variable height t=A(i)
% rust properties according to [1]

% Contact: Anton Tkachuk anton.tkachuk@kau.se
% initial code: 25.06.2024
% last update: 17.09.2025

% units kg-cm-ms => derived Pressure = 10^8 Pa, Force = 10 kN

% Note: works only with mosek or sedumi solvers


%
% Literature 
%
% [1] Rodríguez-Yáñez, J. E., Batlle, S. F., Sanabria-Chinchilla, J., & 
% Rojas-Marín, J. F. (2023). Combined effect of the exposure angle and face 
% orientation on the atmospheric corrosion behavior of low carbon steel. 
% Electrochimica Acta, 439, 141567.
% [2] Kanno's textbook on rest

%% 1. Input data
clear; close all;
% 1.1 material properties
% 1.1.1 mechanical properties
E=2.10e3;      % elasticity kg/(cm*ms^2)
rho=7.800e-3;  % density, kg/cm^3
sigma_admissible = 3.50; % admissible stress in tension and compression 10^8 Pa

% 1.1.2 corrosion parameters, consistent with [1]
a=-0.120e-4 % cm/annum/degree
b=20e-4     %cm/annum (average between skyward/groundward)
lifereq=60  % life requirement

% 1.2 optimization requirements
taumax=6;   % compliance in 10 kN*cm

% design ranges
w=2; % cm, width of all trusses
Alow=2;  %cm,  minimum height
Ahigh=8;  % cm, maximum height: base 16, 
Ascaling=1  % scaling of the design parameters

% 1.3 Dimensions of the grid
nn=30;       % number of cells in x- and y-directions
GridStep=60; % grid size in cm

% 1.4 Generation of the grid
nele = 4*nn^2+2*nn
ndof = 2*(nn+1)^2

% 1.4.2 location of nodes, horizontal
x=zeros((nn+1)^2,1);
y=zeros((nn+1)^2,1);
for i=1:(nn+1)
    for j=1:(nn+1)
        x(j+(nn+1)*(i-1)) = (j-1)*GridStep;
        y(j+(nn+1)*(i-1)) = (i-1)*GridStep;
    end
end

% 1.5 load cases
nlc=3; %
F1=zeros(ndof,1);
F1(ndof) = -15e+0; % tuamax ca. 6
F1(ndof-2*nn) = -15e+0;
F2=zeros(ndof,1);
F2(ndof-1) = -7e+0; % tuamax ca. 6
F2(ndof-2*nn-1) = -7e+0;
F3=zeros(ndof,1);
F3((ndof-2*nn-1):2:(ndof-1)) = (12e+0)/nn; % tuamax ca. 6


F=[F1 F2 F3];


% 1.6 connectivities
ID = zeros(4*nn^2+2*nn,2);

for i=1:(nn+1)
    for j=1:nn
        ID(j+nn*(i-1),1) = j+(nn+1)*(i-1);
        ID(j+nn*(i-1),2) = j+(nn+1)*(i-1)+1;
    end
end
for i=1:nn
    for j=1:(nn+1)
        ID(j+(nn+1)*(i-1)+nn^2+nn,1) = j+(nn+1)*(i-1);
        ID(j+(nn+1)*(i-1)+nn^2+nn,2) = j+(nn+1)*(i-1)+nn+1;
    end
end

for i=1:nn
    for j=1:nn
        ID(j+nn*(i-1)+2*nn^2+2*nn,1) = j+(nn+1)*(i-1);
        ID(j+nn*(i-1)+2*nn^2+2*nn,2) = j+(nn+1)*(i-1)+nn+2;
        ID(j+nn*(i-1)+3*nn^2+2*nn,1) = j+(nn+1)*(i-1)+1;
        ID(j+nn*(i-1)+3*nn^2+2*nn,2) = j+(nn+1)*(i-1)+nn+1;
    end
end

% 1.7 displacement boundary conditions
ind_DBC=[1 2:2:(2*nn+2) 5*(2*nn+2)+1];
full_ind = [1:ndof];
free_ind = setdiff(full_ind,ind_DBC);

%% 2. prepare matrices Qred,Jred for solution process

[Qred,Jred,Kmax,nfree,Fred,Lvector,Aredtop,Aredside] = SDTruss2D_AssembleMatrRustMultLoad(E,rho,x,y,ndof,nele,ID,ind_DBC,F,Alow,Ahigh,w,a,b,lifereq);

% 2.1 check spectrum and compliance before optimization
A = Alow*ones(nele,1);
A(1) = Ahigh;

Qred=sparse(Qred);
%compliance
Kred=Qred*diag(A-Aredtop-A.*Aredside)*Qred';
Dred=Kred\Fred;
compliance_base = Fred'*((Qred*diag(A)*Qred')\Fred);
%spectrum
Mred=diag(Jred*A);
omega2 = eig(Kred,Mred);
f = omega2.^(0.5)/(2*pi);

%% 3. vizualization before optimization
Truss_thickness_plot2D_colormap(x,y,ID,A)


%% 4. call solver %%%%
% with corrosion correction
[A,tau,tot_it2,tot_time2]=MinVolumeGivenComplRustMultLoad(nele,nfree,taumax,Ascaling,Alow,Ahigh,Kmax,Qred,nlc,Fred,Lvector,Aredtop,Aredside);

% witout corrosion correction
%[AwoRed,tauwoRed,tot_it3,tot_time3]=MinVolumeGivenComplRustMultLoad(nele,nfree,taumax,Ascaling,Alow,Ahigh,Kmax,Qred,nlc,Fred,Lvector,zeros(nele,1),zeros(nele,1));

%% 5. checking solution for statisfaction of the compliance constraint
compliance_opt = Fred'*((Qred*diag(A-Aredtop-A.*Aredside+Aredtop.*Aredside)*Qred')\Fred);

[lam,zeta]=eig(compliance_opt)

mtot=Lvector'*A*w*rho

%mtotwoRust=Lvector'*AwoRed*w*rho

%% 6. Vizualization
Truss_thickness_plot2D_colormap(x,y,ID,w*A)

%%
%Truss_thickness_plot2D_colormap(x,y,ID,w*(A-AwoRed))

%Truss_thickness_plot2D_colormap(x,y,ID,w*AwoRed)

%(AwoRed-A)

%m_add=Lvector'*(A-AwoRed)*rho*w

%m_add_rel=100*m_add/mtotwoRust

deterHigh=(Ahigh-(b+a*0)*lifereq)*(w-2*(b+a*90)*lifereq)/(Ahigh*w)
deterLow=(Alow-(b+a*0)*lifereq)*(w-2*(b+a*90)*lifereq)/(Alow*w)

deterHigh90=(Ahigh-(b+a*90)*lifereq)*(w-2*(b+a*90)*lifereq)/(Ahigh*w)
deterLow90=(Alow-(b+a*90)*lifereq)*(w-2*(b+a*90)*lifereq)/(Alow*w)

effectiv_area= w*(A-Aredtop-A.*Aredside+Aredtop.*A.*Aredside);

Truss_thickness_plot2D_colormap(x,y,ID,effectiv_area)

%save ("30x30mosek_scaled_objective.dat","AwoRed","mtotwoRust","tot_time3")

%load('30x30mosek_scaled_objective.mat')

%% number of memebers in groups
numEntriesBelowThreshold = sum(A < Alow*1.001)
numEntriesAboveThreshold = sum(A > Ahigh*0.999)
nele - numEntriesBelowThreshold - numEntriesAboveThreshold 

%% 7. buckling estimates and maximum stress estimates
% LC1
[lambda_crit,sigma_min,sigma_max,local_buckling_sf] =SDTruss2D_AssembleBuckl(E,x,y,ndof,nele,ID,ind_DBC,F1,A,w); % stress in 10^8 Pa

% sigma max
stress_flag=(sigma_admissible>sigma_max) && (sigma_admissible<-sigma_min);

local_buckling_sf_max=max(local_buckling_sf)

% LC2
[lambda_crit2,sigma_min2,sigma_max2,local_buckling_sf2] =SDTruss2D_AssembleBuckl(E,x,y,ndof,nele,ID,ind_DBC,F2,A,w); % stress in 10^8 Pa

% sigma max
stress_flag2=(sigma_admissible>sigma_max2) && (sigma_admissible<-sigma_min2);

local_buckling_sf_max2=max(local_buckling_sf2)

% LC3
[lambda_crit3,sigma_min3,sigma_max3,local_buckling_sf3] =SDTruss2D_AssembleBuckl(E,x,y,ndof,nele,ID,ind_DBC,F3,A,w); % stress in 10^8 Pa

% sigma max
stress_flag3=(sigma_admissible>sigma_max3) && (sigma_admissible<-sigma_min3);

local_buckling_sf_max3=max(local_buckling_sf3)

%%
% LC worrst scenario 1
Fref=F1*lam(1,3)+F2*lam(2,3)+F3*lam(3,3);

[lambda_crit4,sigma_min4,sigma_max4,local_buckling_sf4] =SDTruss2D_AssembleBuckl(E,x,y,ndof,nele,ID,ind_DBC,Fref,A,w); % stress in 10^8 Pa

% sigma max
stress_flag4=(sigma_admissible>sigma_max4) && (sigma_admissible<-sigma_min4);

local_buckling_sf_max4=max(local_buckling_sf4)

%%
% LC worrst scenario 2
Fref2=F1*lam(1,2)+F2*lam(2,2)+F3*lam(3,2);

[lambda_crit5,sigma_min5,sigma_max5,local_buckling_sf5] =SDTruss2D_AssembleBuckl(E,x,y,ndof,nele,ID,ind_DBC,Fref2,A,w); % stress in 10^8 Pa

% sigma max
stress_flag4=(sigma_admissible>sigma_max5) && (sigma_admissible<-sigma_min5);

local_buckling_sf_max5=max(local_buckling_sf5)