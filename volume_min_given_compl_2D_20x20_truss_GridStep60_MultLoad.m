% This file demonstrates design for truss structures including rust
% All trusses have width w and variable height t=A(i)
% rust properties according to [1]


%
% Literature 
%
% [1] Rodríguez-Yáñez, J. E., Batlle, S. F., Sanabria-Chinchilla, J., & 
% Rojas-Marín, J. F. (2023). Combined effect of the exposure angle and face 
% orientation on the atmospheric corrosion behavior of low carbon steel. 
% Electrochimica Acta, 439, 141567.
% [2] Kanno's textbook on rest
% 
% Anton Tkachuk anton.tkachuk@kau.se
% 25.06.2024
% units kg-cm-ms

clear; close all;
% material properties
E=2.10e3;      % elasticity kg/(cm*ms^2)
rho=7.800e-3;  % density, kg/cm^3
sigma_admissible = 3.50;% 10^8 Pa
% corrosion parameters, consistent with [1]
a=-0.120e-4 % cm/annum/degree
b=20e-4 %cm/annum (average between skyward/groundward)
lifereq=120 % life requirement

% optimization requirements
Ascaling=1  % scaling of the design parameters
taumax=20/(200/60)  % compliance

CVXSOLVER = 'sdpt'
CVXPRECISION = 'high'

% design ranges
w=2; % cm, width of all trusses
Alow=2;  %cm,  minimum height
Ahigh=8;  % cm, maximum height: base 16, 

% dimensions
nn=20
nele = 4*nn^2+2*nn
ndof = 2*(nn+1)^2

% location of nodes, horizontal

GridStep=60;

x=zeros((nn+1)^2,1)
y=zeros((nn+1)^2,1)
for i=1:(nn+1)
    for j=1:(nn+1)
        x(j+(nn+1)*(i-1)) = (j-1)*GridStep;
        y(j+(nn+1)*(i-1)) = (i-1)*GridStep;
    end
end



%load cases

% nlc=1;
% F=zeros(ndof,1);
% F((ndof-2*nn-1):2:(ndof-1)) = (18e+0)/nn; % tuamax ca. 20
% F=zeros(ndof,1);
% F(ndof) = -15e+0; % tuamax ca. 20
% F(ndof-2*nn) = -15e+0;
% F=zeros(ndof,1);
% F(ndof-1) = -10e+0; % tuamax ca. 20
% F(ndof-2*nn-1) = -10e+0;


% nlc=2;
% F1=zeros(ndof,1);
% F1(ndof) = -15e+0; % tuamax ca. 20
% F1(ndof-2*nn) = -15e+0;
% F2=zeros(ndof,1);
% F2(ndof-1) = -10e+0; % tuamax ca. 20
% F2(ndof-2*nn-1) = -10e+0;
% 
% F=[F1 F2]

 nlc=3;
F1=zeros(ndof,1);
F1(ndof) = -15e+0; % tuamax ca. 20
F1(ndof-2*nn) = -15e+0;
F2=zeros(ndof,1);
F2(ndof-1) = -7e+0; % tuamax ca. 20
F2(ndof-2*nn-1) = -7e+0;
F3=zeros(ndof,1);
F3((ndof-2*nn-1):2:(ndof-1)) = (12e+0)/nn; % tuamax ca. 20


F=[F1 F2 F3];


% number of degrees of freedom and connectivities



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

ind_DBC=[1 2:2:(2*nn+2)];
full_ind = [1:ndof];
free_ind = setdiff(full_ind,ind_DBC);

%% prepare matrices Qred,Jred for solution process

[Qred,Jred,Kmax,nfree,Fred,Lvector,Aredtop,Aredside] = SDTruss2D_AssembleMatrRustMultLoad(E,rho,x,y,ndof,nele,ID,ind_DBC,F,Alow,Ahigh,w,a,b,lifereq);
% and M0
%M0=zeros(nfree);
% M0(15,15) = m0;
% M0(16,16) = m0;

% check spectrum and compliance before optimization
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

%% vizualization before optimization
Truss_thickness_plot2D_colormap(x,y,ID,A)


%% call solver %%%%

[A,tau,tot_it2,tot_time2]=MinVolumeGivenComplRustMultLoad(nele,nfree,taumax,Ascaling,Alow,Ahigh,Kmax,Qred,nlc,Fred,Lvector,Aredtop,Aredside)


[AwoRed,tauwoRed,tot_it3,tot_time3]=MinVolumeGivenComplRustMultLoad(nele,nfree,taumax,Ascaling,Alow,Ahigh,Kmax,Qred,nlc,Fred,Lvector,zeros(nele,1),zeros(nele,1))


compliance_opt = Fred'*((Qred*diag(A-Aredtop-A.*Aredside+Aredtop.*A.*Aredside)*Qred')\Fred)

mtot=Lvector'*A*w*rho

mtotwoRust=Lvector'*AwoRed*w*rho

Truss_thickness_plot2D_colormap(x,y,ID,w*A)

Truss_thickness_plot2D_colormap(x,y,ID,w*(A-AwoRed))

Truss_thickness_plot2D_colormap(x,y,ID,w*AwoRed)

%(AwoRed-A)

m_add=Lvector'*(A-AwoRed)*rho*w

m_add_rel=100*m_add/mtotwoRust

deterHigh=(Ahigh-(b+a*0)*lifereq)*(w-2*(b+a*90)*lifereq)/(Ahigh*w)
deterLow=(Alow-(b+a*0)*lifereq)*(w-2*(b+a*90)*lifereq)/(Alow*w)

deterHigh90=(Ahigh-(b+a*90)*lifereq)*(w-2*(b+a*90)*lifereq)/(Ahigh*w)
deterLow90=(Alow-(b+a*90)*lifereq)*(w-2*(b+a*90)*lifereq)/(Alow*w)

effectiv_area= w*(A-Aredtop-A.*Aredside+Aredtop.*A.*Aredside);

Truss_thickness_plot2D_colormap(x,y,ID,effectiv_area)

%% buckling estimates
%lambda_crit =SDTruss2D_AssembleBuckl(E,x,y,ndof,nele,ID,ind_DBC,F1,A,w);

%% buckling estimates and maximum stress estimates
[lambda_crit,sigma_min,sigma_max,local_buckling_sf] =SDTruss2D_AssembleBuckl(E,x,y,ndof,nele,ID,ind_DBC,F,A,w); % stress in 10^8 Pa

% sigma max
stress_flag=(sigma_admissible>sigma_max) && (sigma_admissible<-sigma_min);

local_buckling_sf_max=max(local_buckling_sf)
