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
% 13.06.2024
% units kg-cm-ms

clear; close all;
E=2.10e3;
rho=7.800e-3;
%rank_tol=1e+3
Ascaling=1
taumax=1.5
m0=40
CVXSOLVER = 'sdpt'
CVXPRECISION = 'high'

%w=2; % wifth of all trusses

% rethink parameters
a=-0.120e-4 % cm/annum/degree
b=20e-4 %cm/annum (average between skyward/groundward)
lifereq=60 % years up to 175 years with total mass increase by 342 kg from 625 kg

niter=7;
alpha=2; %aspect ratio
Alow=4  % minimum height
Ahigh=16  %maximum height: base 16, 


 x= 100*[0; 0; 0; 1; 1; 1; 2; 2; 2; 3; 3; 3]
 y= 100*[2; 1; 0; 2; 1; 0; 2; 1; 0; 2; 1; 0]

% y= 100*[0; 0; 0; 1; 1; 1; 2; 2; 2; 3; 3; 3]
% x= -100*[2; 1; 0; 2; 1; 0; 2; 1; 0; 2; 1; 0]

ifbaseflag=true
omega2base=(2*pi*0.100)^2


nele = 33
ndof = 24

F=[0;0; 0;0; 0;0;  0;0; 0;0; 0;0;   0;0; 0;0; 0;0;   0;0; 0.0e+1; 0.8e+1;     0;0;      ]
%F=[0;0; 0;0; 0;0;  0;0; 0;0; 0;0;   0;0; 0;0; 0;0;   0;0; 1.2e+1;0;     0;0;      ]
%F=[0;0; 0;0; 0;0;  0;0; 0;0; 0;0;   0;0; 0;0; 0;0;   0;0; 0;1.2e+1;     0;0;      ] % rotation by 90 degrees counterclockwise

%ID = [1 4; 4 5; 5 6; 3 6; 2 5; 1 5; 2 4; 1 6; 4 3; 2 6; 3 5;      4 7; 7 8; 8 9; 6 9; 5 8; 4 8; 5 7; 4 9; 7 6; 5 9; 6 8;   7 10; 10 11; 11 12; 9 12; 8 11; 7 11; 8 11; 7 12; 10 9; 8 12; 9 11]
ID = [1 4; 4 5; 5 6; 3 6; 2 5; 1 5; 2 4; 1 6; 4 3; 2 6; 3 5;      4 7; 7 8; 8 9; 6 9; 5 8; 4 8; 5 7; 4 9; 7 6; 5 9; 6 8;   7 10; 10 11; 11 12; 9 12; 8 11; 7 11; 8 10; 7 12; 10 9; 8 12; 9 11]
ind_DBC=[1;2;3;4;5;6]

%% prepare matrices Qred,Jred for solution process

[Qred,Jred,Kmax,nfree,Fred,Lvector,Aredtop,Aredside,Alowvec,Ahighvec] = SDTruss2D_AssembleMatrRust_fixed_aspect_ratio(E,rho,x,y,ndof,nele,ID,ind_DBC,F,Alow,Ahigh,a,b,lifereq,alpha)
% and M0
M0=zeros(nfree)
M0(15,15) = m0
M0(16,16) = m0

% chack spectrum and compliance before optimization
A = Alow*ones(nele,1);
A(1) = Ahigh;



%compliance
Kred=Qred*diag(A)*Qred'
Dred=Kred\Fred
compliance_base = Fred'*((Qred*diag(A)*Qred')\Fred)
%spectrum
Mred=diag(Jred*A)
omega2 = eig(Kred,Mred);
f = omega2.^(0.5)/(2*pi)

%% vizualization before optimization
Truss_thickness_plot2D_colormap(x,y,ID,A)


%% call solver %%%%

[Ar,tau,tot_it2,tot_time2,ArNormEvol2]=MinVolumeGivenComplAndFreq_fixed_aspect_ratio(nele,nfree,taumax,Ascaling,Alowvec,Ahighvec,Kmax,Qred,Fred,Lvector,Aredtop,Aredside,alpha,niter)

[Ar2,tau2,tot_it22,tot_time22,ArNormEvol22]=MinVolumeGivenComplAndFreq_fixed_aspect_ratio_weight2(nele,nfree,taumax,Ascaling,Alowvec,Ahighvec,Kmax,Qred,Fred,Lvector,Aredtop,Aredside,alpha,niter)

[AwoRed,tauwoRed,tot_it3,tot_time3,ArNormEvol3]=MinVolumeGivenComplAndFreq_fixed_aspect_ratio(nele,nfree,taumax,Ascaling,Alow,Ahigh,Kmax,Qred,Fred,Lvector,zeros(nele,1),zeros(nele,1),alpha,1)


omega2 = eigs(Qred*diag(A-Aredtop-A.*Aredside)*Qred',diag(Jred*A)+ M0,10,'smallestabs')
f = omega2.^(0.5)/(2*pi)

%A=Ar+sqrt(Ar).*(Aredside/sqrt(alpha)+Aredtop*sqrt(alpha))+Aredside.*Aredtop;
A=Ar+0.5*(Aredside+Aredtop).*sqrt(4*Ar -4*Aredside.*Aredtop+ (Aredside+Aredtop).^2)-Aredside.*Aredtop + 0.5*(Aredside+Aredtop).^2;

Dred=(Qred*diag(Ar)*Qred')\Fred;
compliance_opt = Fred'*(Dred);

mtot=Lvector'*A*rho

mtotwoRust=Lvector'*AwoRed*rho

Truss_thickness_plot2D_colormap(x,y,ID,A)

Truss_thickness_plot2D_colormap(x,y,ID,(A-AwoRed))

Truss_thickness_plot2D_colormap(x,y,ID,AwoRed)

(AwoRed-A)

m_add=Lvector'*(A-AwoRed)*rho

m_add_rel =100*m_add/mtotwoRust

diff_in_designs=norm(Ar-Ar2,2)/norm(Ar2,2)

figure
hold on
plot(log(ArNormEvol2))
plot(log(ArNormEvol22))

%deterHigh=(Ahigh-(b+a*0)*lifereq)*(w-2*(b+a*90)*lifereq)/(Ahigh)
%deterLow=(Alow-(b+a*0)*lifereq)*(w-2*(b+a*90)*lifereq)/(Alow)

%deterHigh90=(Ahigh-(b+a*90)*lifereq)*(w-2*(b+a*90)*lifereq)/(Ahigh)
%deterLow90=(Alow-(b+a*90)*lifereq)*(w-2*(b+a*90)*lifereq)/(Alow)

Dred(16)