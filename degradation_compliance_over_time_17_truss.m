% This file demonstrates deterioration of compliance for a single load 
% for truss structures due to rust rust
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
% 09.09.2025
% units kg-cm-ms
% force 10 kN, complainace 100 J

clear; close all;
% material properties
E=2.06842e3;      % elasticity kg/(cm*ms^2)
rho=7.418e-3;  % density, kg/cm^3


lifereq_max=120 % life requirement
stepmax=60; %steps in which degradation is computed

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
F=[0;0; 0;0; 0;0;  0;0; 0;0; 0;0;   0;0; 0;0;  0.0e+1; 5.0e+1] % unit, 10 kN



% number of degrees of freedom and connectivities
nele = 17
ndof = 18
%     1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17 
ID = [2 4; 2 3; 1 3; 3 4; 4 6; 4 5; 3 5; 5 6; 6 8; 6 7; 5 7; 7 8; 8 9; 7 9; 1 4; 3 6; 5 8]
ind_DBC=[1;2;3;4]

% corrosion parameters, consistent with [1]
a=-0.120e-4 % 0 cm/annum/degree
b=20e-4 %cm/annum (average between skyward/groundward)
% cross-sections, for case 1
A =  [ 20.0000
    5.7225
   20.0000
    5.0000
   20.0000
    8.1217
   20.0000
    5.0000
   15.9036
    5.0000
   12.1915
    5.0000
   13.1931
    9.4450
    7.6458
    5.2466
    9.3090]

compliance_base_case1 = zeros(stepmax,1)
for i=1:stepmax
    % prepare matrices Qred,Jred for solution process
    [Qred,Jred,Kmax,nfree,Fred,Lvector,Aredtop,Aredside] = SDTruss2D_AssembleMatrRustMultLoad(E,rho,x,y,ndof,nele,ID,ind_DBC,F,Alow,Ahigh,w,a,b,lifereq_max*i/stepmax)
    
    %compliance
    Kred=Qred*diag(A-Aredtop-A.*Aredside)*Qred'
    Dred=Kred\Fred
    compliance_base_case1(i) = Fred'*((Qred*diag(A-Aredtop-A.*Aredside+Aredtop.*Aredside)*Qred')\Fred)

end

% case 2 solution and parameters
a= 0 % cm/annum/degree
b=20e-4 %cm/annum (average between skyward/groundward)
A=[  20.0000
    6.0522
   20.0000
    5.0000
   20.0000
    8.9153
   20.0000
    5.0000
   17.2008
    5.0000
   12.8326
    5.0000
   14.1716
   10.0912
    8.3595
    5.4964
   10.2946]
compliance_base_case2 = zeros(stepmax,1)
for i=1:stepmax
    % prepare matrices Qred,Jred for solution process
    [Qred,Jred,Kmax,nfree,Fred,Lvector,Aredtop,Aredside] = SDTruss2D_AssembleMatrRustMultLoad(E,rho,x,y,ndof,nele,ID,ind_DBC,F,Alow,Ahigh,w,a,b,lifereq_max*i/stepmax)
    
    %compliance
    Kred=Qred*diag(A-Aredtop-A.*Aredside)*Qred'
    Dred=Kred\Fred
    compliance_base_case2(i) = Fred'*((Qred*diag(A-Aredtop-A.*Aredside+Aredtop.*Aredside)*Qred')\Fred)

end

%% vizualization before optimization
Truss_thickness_plot2D_colormap(x,y,ID,A)

years=1:lifereq_max/stepmax:lifereq_max;
figure
hold on
plot(years,0.1*compliance_base_case1,LineWidth=2)
plot(years,0.1*compliance_base_case2,'-.',LineWidth=2)
plot([0,lifereq_max],0.1*[taumax,taumax],'--',LineWidth=2)
xlabel('time in years')
ylabel('compliance, kN*m')
legend({'case 1', 'case 2','required compliance, max'},'Location','southeast')




