function [lambda_crit,sigma_max,sigma_min,local_buckling_sf] = SDTruss2D_AssembleBuckl(E,x,y,ndof,nele,ID,ind_DBC,F,A,w)
% return estimates for a global buckling, maximum and minimum stress, and
% elementwise estimate for local buckling for 2D truss structure
% Contact: Anton Tkachuk anton.tkachuk@kau.se
% initial code: 15.09.2025
% last update: 17.09.2025
% Input:
% E - elasticity modulus, scalar
% rho - density, scalar
% x,y - nodal location, two vector of length nele
% ndof,nele - number of all degrees of nodes and elements, scalar
% ID - connectivity of nodes in mesh, nele x 2
% ind_BC - a list of notes that is fixed, list
% F - vector with one load cases, ndof
% A - heights of elements
% w - width of elements
% 
% Output:
% lambda_crit - the lowest positive eigenvalue det (Kred - Lambda KGred) = 0
% sigma_max,sigma_min - maximum and minimal (compression) stress in trusses 
% local_buckling_sf - estimate with 2nd Euler case for individual buckling
% of elements (value above 1 means, that the element buckles)

full_ind = [1:ndof];
free_ind = setdiff(full_ind,ind_DBC);

Qfull = zeros(ndof,nele);
Lvector = zeros(nele,1);
Cfull = diag(nele);

% Assembly
for i=1:nele
    x1 = x(ID(i,1));
    x2 = x(ID(i,2));
    y1 = y(ID(i,1));
    y2 = y(ID(i,2));
    L0=sqrt((x2-x1)^2+(y2-y1)^2);
    Lvector(i,1)=L0;
    Biloc = [ -(x2-x1)  -(y2-y1)  (x2-x1)  (y2-y1)]/L0;
    sind = [2*ID(i,1)-1,2*ID(i,1),2*ID(i,2)-1,2*ID(i,2)];
    Qfull(sind,i) = Biloc;
    Cfull(i,i) = E*A(i)*w*L0;
end

% Application of boundary conditions
Qred = Qfull(free_ind,:);
Kred=Qred*Cfull*Qred';
Fred=F(free_ind,1);

% Displacement for a given load vector
DDred= Kred\Fred;

% elemental forces
Nforce=Cfull*Qred'*DDred;

% maximum and minimal stress
sigma_elem=(Nforce./A)/w;
sigma_max = max(sigma_elem);
sigma_min = min(sigma_elem);

%local buckling, safety factor
local_buckling_sf=zeros(nele,1);
for i=1:nele
    if (sigma_elem(i)<0)
        ri=min(w,A(i))/sqrt(12.); %radius of gyration for A
        %ri=sqrt(w*A(i))/sqrt(12.); %radius of gyration for A
        lambda_f=Lvector(i)/ri;
        local_buckling_sf(i) = -sigma_elem(i)/(pi^2*E/lambda_f^2);
    else
        local_buckling_sf(i)=0;
    end
end

%fprintf('%f10.5',min(Nforce))

% global buckling
KG=zeros(ndof,ndof);
kgpattern=[1  0  -1  0; 0  1 0 -1; -1  0  1  0; 0  -1  0  1];

for i=1:nele
   sind = [2*ID(i,1)-1,2*ID(i,1),2*ID(i,2)-1,2*ID(i,2)];
   KG(sind,sind) =KG(sind,sind)+ Nforce(i)/Lvector(i,1)*kgpattern;
end

KGred=KG(free_ind,free_ind);


[v,e] = eig(Kred,-KGred,'vector');

% skipping infinity values and finding the smallest positive
sind = find(isfinite(e) & e>0);
[d,sind] = sort(e(sind));

lambda_crit=d(1);


end