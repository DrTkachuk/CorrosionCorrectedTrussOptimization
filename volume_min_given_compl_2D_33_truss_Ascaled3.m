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
Ascaling=1
taumax=1.5

CVXSOLVER = 'sdpt'
CVXPRECISION = 'high'

w=2; % wifth of all trusses

% rethink parameters
a=-0.120e-4 % cm/annum/degree
b=20e-4 %cm/annum (average between skyward/groundward)
lifereq=250 % years up to 175 years with total mass increase by 342 kg from 625 kg

Alow=2  % minimum height
Ahigh=8  %maximum height: base 16, 
% with 18 mtot reduces to 640/806 kg wo/with rust over 175 years
% with 20 mtot reduces to 632/770 kg wo/with rust over 175 years

%
x= 100*[0; 0; 0; 1; 1; 1; 2; 2; 2; 3; 3; 3]
y= 100*[2; 1; 0; 2; 1; 0; 2; 1; 0; 2; 1; 0]

 % y= 100*[0; 0; 0; 1; 1; 1; 2; 2; 2; 3; 3; 3]
 % x= -100*[2; 1; 0; 2; 1; 0; 2; 1; 0; 2; 1; 0]

% ifbaseflag=true
% omega2base=(2*pi*0.100)^2


nele = 33
ndof = 24
F=[0;0; 0;0; 0;0;  0;0; 0;0; 0;0;   0;0; 0;0; 0;0;   0;0; 0.0e+1; 0.8e+1;     0;0;      ]
%F=[0;0; 0;0; 0;0;  0;0; 0;0; 0;0;   0;0; 0;0; 0;0;   0;0; 0.8e+1; 0.0e+1;     0;0;      ] % rotation by 90 degrees counterclockwise
%F=[0;0; 0;0; 0;0;  0;0; 0;0; 0;0;   0;0; 0;0; 0;0;   0;0; 1.2e+1; 0.6e+1;     0;0;      ]
%F=[0;0; 0;0; 0;0;  0;0; 0;0; 0;0;   0;0; 0;0; 0;0;   0;0; 0;1.2e+1;     0;0;      ] % rotation by 90 degrees counterclockwise

%ID = [1 4; 4 5; 5 6; 3 6; 2 5; 1 5; 2 4; 1 6; 4 3; 2 6; 3 5;      4 7; 7 8; 8 9; 6 9; 5 8; 4 8; 5 7; 4 9; 7 6; 5 9; 6 8;   7 10; 10 11; 11 12; 9 12; 8 11; 7 11; 8 11; 7 12; 10 9; 8 12; 9 11]
ID = [1 4; 4 5; 5 6; 3 6; 2 5; 1 5; 2 4; 1 6; 4 3; 2 6; 3 5;      4 7; 7 8; 8 9; 6 9; 5 8; 4 8; 5 7; 4 9; 7 6; 5 9; 6 8;   7 10; 10 11; 11 12; 9 12; 8 11; 7 11; 8 10; 7 12; 10 9; 8 12; 9 11]
ind_DBC=[1;2;3;4;5;6]

%% prepare matrices Qred,Jred for solution process

[Qred,Jred,Kmax,nfree,Fred,Lvector,Aredtop,Aredside] = SDTruss2D_AssembleMatrRust(E,rho,x,y,ndof,nele,ID,ind_DBC,F,Alow,Ahigh,w,a,b,lifereq)
% and M0
M0=zeros(nfree)

% chack spectrum and compliance before optimization
A = Ahigh*ones(nele,1);




%compliance
Kred=Qred*diag(A-Aredtop-A.*Aredside+Aredtop.*Aredside)*Qred'
Dred=Kred\Fred
compliance_base_withRust = Fred'*Dred
compliance_base_woRust = Fred'*((Qred*diag(A)*Qred')\Fred)
%spectrum
Mred=diag(Jred*A)
omega2 = eig(Kred,Mred);
f = omega2.^(0.5)/(2*pi)
A(1) = Alow;

%% vizualization before optimization
Truss_thickness_plot2D_colormap(x,y,ID,A)


%% call solver %%%%

[A,tau,tot_it2,tot_time2]=MinVolumeGivenComplRust(nele,nfree,taumax,Ascaling,Alow,Ahigh,Kmax,Qred,Fred,Lvector,Aredtop,Aredside)


[AwoRed,tauwoRed,tot_it3,tot_time3]=MinVolumeGivenComplRust(nele,nfree,taumax,Ascaling,Alow,Ahigh,Kmax,Qred,Fred,Lvector,zeros(nele,1),zeros(nele,1))


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

figure
plot ([5,10,20,40,60,120, 180,240],[0.76 , 1.54 , 2.72, 6.51   ,10.16    , 27.1 , 56.03, 114.64])