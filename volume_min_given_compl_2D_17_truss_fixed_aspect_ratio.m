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
E=2.06842e3;      % elasticity kg/(cm*ms^2)
rho=7.418e-3;  % density, kg/cm^3
Ascaling=1
taumax=200
CVXSOLVER = 'sdpt'
CVXPRECISION = 'high'

%w=2; % wifth of all trusses

% rethink parameters
%a=0
a=-0.120e-4 % cm/annum/degree
b=20e-4 %cm/annum (average between skyward/groundward)
lifereq=60 % years up to 175 years with total mass increase by 342 kg from 625 kg

niter=7;
alpha=1; %aspect ratio
Alow=25;  % area limit below
Ahigh=100  %area limit above




 x= 254*[0; 0; 1; 1; 2; 2; 3; 3; 4]
 y= 254*[0; 1; 0; 1; 0; 1; 0; 1; 0]

% y= 100*[0; 0; 0; 1; 1; 1; 2; 2; 2; 3; 3; 3]
% x= -100*[2; 1; 0; 2; 1; 0; 2; 1; 0; 2; 1; 0]

ifbaseflag=true
omega2base=(2*pi*0.100)^2


% number of degrees of freedom and connectivities
nele = 17
ndof = 18

%nlc=1;
F=[0;0; 0;0; 0;0;  0;0; 0;0; 0;0;   0;0; 0;0;  0.0e+1; 5.0e+1]
%F=[0;0; 0;0; 0;0;  0;0; 0;0; 0;0;   0;0; 0;0; 0;0;   0;0; 1.2e+1;0;     0;0;      ]
%F=[0;0; 0;0; 0;0;  0;0; 0;0; 0;0;   0;0; 0;0; 0;0;   0;0; 0;1.2e+1;     0;0;      ] % rotation by 90 degrees counterclockwise

%     1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17 
ID = [2 4; 2 3; 1 3; 3 4; 4 6; 4 5; 3 5; 5 6; 6 8; 6 7; 5 7; 7 8; 8 9; 7 9; 1 4; 3 6; 5 8]
ind_DBC=[1;2;3;4]

%% prepare matrices Qred,Jred for solution process

[Qred,Jred,Kmax,nfree,Fred,Lvector,Aredtop,Aredside,Alowvec,Ahighvec] = SDTruss2D_AssembleMatrRust_fixed_aspect_ratio(E,rho,x,y,ndof,nele,ID,ind_DBC,F,Alow,Ahigh,a,b,lifereq,alpha)
% and M0
M0=zeros(nfree)


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

[Ar2,tau2,tot_it22,tot_time22,ArNormEvol22]=MinVolumeGivenComplAndFreq_fixed_aspect_ratio_scaling_off(nele,nfree,taumax,Ascaling,Alowvec,Ahighvec,Kmax,Qred,Fred,Lvector,Aredtop,Aredside,alpha,niter)

%scaling of E modulus
uc=1e4;
[Ar3,tau3,tot_it23,tot_time23,ArNormEvol23]=MinVolumeGivenComplAndFreq_fixed_aspect_ratio_scaling_off(nele,nfree,taumax,Ascaling,Alowvec,Ahighvec,Kmax,Qred*uc,Fred*uc,Lvector,Aredtop,Aredside,alpha,niter)


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

%deterHigh=(Ahigh-(b+a*0)*lifereq)*(w-2*(b+a*90)*lifereq)/(Ahigh)
%deterLow=(Alow-(b+a*0)*lifereq)*(w-2*(b+a*90)*lifereq)/(Alow)

%deterHigh90=(Ahigh-(b+a*90)*lifereq)*(w-2*(b+a*90)*lifereq)/(Ahigh)
%deterLow90=(Alow-(b+a*90)*lifereq)*(w-2*(b+a*90)*lifereq)/(Alow)

OutputFoPaper=sprintf('%3.1f &',A)


sprintf('%3.1f &',mtot, m_add,m_add_rel)

%% effect of scaling (paper has a graph by SDPT3 with medium accuracy)
figure
hold on
plot(log10(ArNormEvol2),'-b')
plot(log10(ArNormEvol22),'-o')
plot(log10(ArNormEvol23),'-*')
xlabel('#iter')
ylabel('log_1_0|s^j^+^1-s^j|')
legend({'scaling on','scaling off', 'scaling off, SI units'},...
    'Location','southwest') 
hold off

%save ("17goodconv.dat","goodConv") 

Dred(14)
