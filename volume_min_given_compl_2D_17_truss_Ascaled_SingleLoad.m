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
% [3] Berke, L., & Khot, N. S. (1987). Structural optimization using optimality criteria.
% In Computer aided optimal design: structural and mechanical systems (pp. 271-311).
% Berlin, Heidelberg: Springer Berlin Heidelberg.
% 
% Anton Tkachuk anton.tkachuk@kau.se
% 25.06.2024
% units kg-cm-ms
% force 10 kN, complainace 100 J

clear; close all;
% material properties
E=2.06842e3;      % elasticity kg/(cm*ms^2)
rho=7.418e-3;  % density, kg/cm^3
sigma_admissible = 3.44738;% 10^8 Pa
% corrosion parameters, consistent with [1]
%a=0
a=0%-0.120e-4 % cm/annum/degree
b=20e-4 %cm/annum (average between skyward/groundward)
lifereq=60 % life requirement

% optimization requirements
Ascaling=1  % scaling of the design parameters
taumax=200  % compliance 
m0=0
CVXSOLVER = 'sdpt'
CVXPRECISION = 'high'

% design ranges
w=5; % cm, width of all trusses
Alow=5;  %cm,  minimum height
Ahigh=20;  % cm, maximum height: base 16, 

% location of nodes, horizontal
 x= 254*[0; 0; 1; 1; 2; 2; 3; 3; 4]
 y= 254*[0; 1; 0; 1; 0; 1; 0; 1; 0]
% location of nodes, vertical


%load cases

nlc=1;
F=[0;0; 0;0; 0;0;  0;0; 0;0; 0;0;   0;0; 0;0;  0.0e+1; -5.0e+1] % unit, 10 kN

% nlc=2;
% F1=[0;0; 0;0; 0;0;  0;0; 0;0; 0;0;   0;0; 0;0;  0.0e+1; 0.8e+1]
% F2=[0;0; 0;0; 0;0;  0;0; 0;0; 0;0;   0;0; 0;0;  1.0e+1; 0.0e+1]
% F=[ F1 F2]

% nlc=3;
% F1=[0;0; 0;0; 0;0;  0;0; 0;0; 0;0;   0;0; 0;0;  0.0e+1; 0.8e+1]
% F2=[0;0; 0;0; 0;0;  0;0; 0;0; 0;0;   0;0; 0;0;  1.0e+1; 0.0e+1]
% F3=[0;0; 0;0; 0;0;  0;0; 0;0; 0;5;   0;0; 0;0;  0.0e+1; 0.0e+1]
% F=[ F1 F2 F3]

% number of degrees of freedom and connectivities
nele = 17
ndof = 18
%     1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17 
ID = [2 4; 2 3; 1 3; 3 4; 4 6; 4 5; 3 5; 5 6; 6 8; 6 7; 5 7; 7 8; 8 9; 7 9; 1 4; 3 6; 5 8]
ind_DBC=[1;2;3;4]

%% prepare matrices Qred,Jred for solution process

[Qred,Jred,Kmax,nfree,Fred,Lvector,Aredtop,Aredside] = SDTruss2D_AssembleMatrRustMultLoad(E,rho,x,y,ndof,nele,ID,ind_DBC,F,Alow,Ahigh,w,a,b,lifereq)
% and M0
M0=zeros(nfree)

% check spectrum and compliance before optimization
A = Alow*ones(nele,1);
A(1) = Ahigh;



%compliance
Kred=Qred*diag(A-Aredtop-A.*Aredside)*Qred'
Dred=Kred\Fred
compliance_base = Fred'*((Qred*diag(A)*Qred')\Fred)
%spectrum
Mred=diag(Jred*A)
omega2 = eig(Kred,Mred);
f = omega2.^(0.5)/(2*pi)

Dred(nfree)

%% vizualization before optimization
Truss_thickness_plot2D_colormap(x,y,ID,A)


%% call solver %%%%

[A,tau,tot_it2,tot_time2]=MinVolumeGivenComplRustMultLoad(nele,nfree,taumax,Ascaling,Alow,Ahigh,Kmax,Qred,nlc,Fred,Lvector,Aredtop,Aredside)


[AwoRed,tauwoRed,tot_it3,tot_time3]=MinVolumeGivenComplRustMultLoad(nele,nfree,taumax,Ascaling,Alow,Ahigh,Kmax,Qred,nlc,Fred,Lvector,zeros(nele,1),zeros(nele,1))


Dredopt=(Qred*diag(A-Aredtop-A.*Aredside+Aredtop.*A.*Aredside)*Qred')\Fred;
Dredopt(14)

compliance_opt = Fred'*((Qred*diag(A-Aredtop-A.*Aredside+Aredtop.*A.*Aredside)*Qred')\Fred)

mtot=Lvector'*A*w*rho

mtotwoRust=Lvector'*AwoRed*w*rho

Truss_thickness_plot2D_colormap(x,y,ID,w*A)

Truss_thickness_plot2D_colormap(x,y,ID,w*(A-AwoRed))

Truss_thickness_plot2D_colormap(x,y,ID,w*AwoRed)

(AwoRed-A)

m_add=Lvector'*(A-AwoRed)*rho*w

m_add_rel=100*m_add/mtotwoRust

deterHigh=(Ahigh-(b+a*0)*lifereq)*(w-2*(b+a*90)*lifereq)/(Ahigh*w)
deterLow=(Alow-(b+a*0)*lifereq)*(w-2*(b+a*90)*lifereq)/(Alow*w)

deterHigh90=(Ahigh-(b+a*90)*lifereq)*(w-2*(b+a*90)*lifereq)/(Ahigh*w)
deterLow90=(Alow-(b+a*90)*lifereq)*(w-2*(b+a*90)*lifereq)/(Alow*w)

effectiv_area= w*(A-Aredtop-A.*Aredside+Aredtop.*A.*Aredside)

Truss_thickness_plot2D_colormap(x,y,ID,effectiv_area)

OutputFoPaper=sprintf('%3.1f &',A*w)


sprintf('%3.1f &',mtot, m_add,m_add_rel)


%% buckling estimates and maximum stress estimates
[lambda_crit,sigma_min,sigma_max,local_buckling_sf] =SDTruss2D_AssembleBuckl(E,x,y,ndof,nele,ID,ind_DBC,F,A,w); % stress in 10^8 Pa

% sigma max
stress_flag=(sigma_admissible>sigma_max) && (sigma_admissible<-sigma_min);