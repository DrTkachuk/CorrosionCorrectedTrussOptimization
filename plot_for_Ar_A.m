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


%w=2; % wifth of all trusses

% rethink parameters
%a=0
a=-0.120e-4 % cm/annum/degree
b=20e-4 %cm/annum (average between skyward/groundward)
lifereq=120 % years up to 175 years with total mass increase by 342 kg from 625 kg

alpha=1; %aspect ratio
Alow=9;  % area limit below
Ahigh=100  %area limit above

w=Alow^0.5;

timeCorr=[0:1:lifereq];
d05=
Ar=@(x) x.^a
figure
plot();