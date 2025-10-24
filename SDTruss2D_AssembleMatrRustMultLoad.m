function [Qred,Jred,Kmax,nfree,Fred,Lvector,Aredtop,Aredside] = SDTruss2D_AssembleMatrRustMultLoad(E,rho,x,y,ndof,nele,ID,ind_DBC,F,Alow,Ahigh,w,a,b,lifereq)
% return assembled matrices for assembly of mass and stiffness matrices
% for 2D truss structure, e.g. Kred = Qred*diag(A)*Qred'; Mred = Jred*diag(A)*Jred'
% Contact: Anton Tkachuk anton.tkachuk@kau.se
% initial code: 17.06.2024
% last update: 17.09.2025
% Input:
% E - elasticity modulus, scalar
% rho - density, scalar
% x,y - nodal location, two vector of length nele
% ndof,nele - number of all degrees of nodes and elements, scalar
% ID - connectivity of nodes in mesh, nele x 2
% ind_BC - a list of notes that is fixed, list
% F - vector with several load cases, ndof x nlc
% Alow,Ahigh - upper and lower limits on height of the elements, scalar
% w - width of elements
% a,b,lifereq - corrosion parameters and required life
% Output:
% Qred - matrix, generating stiffness Kred = Qred*diag(A)*Qred'
% Jred - matrix, generating mass Mred = Jred*diag(A)*Jred'
% Kmax - reference stiffness for scaling
% nfree - number of free degrees of freedom
% Fred - matrix, load vector reduced to free degrees of freedom
% Lvector - vector with element length information, Volume = Lvector'*A*w
% Aredtop,Aredside - reduction of vertical and horizontal dimentions due to
% corrosion
 
full_ind = [1:ndof];
free_ind = setdiff(full_ind,ind_DBC);
nfree=size(free_ind,2);


Qfull = zeros(ndof,nele);
Jfull = zeros(ndof,nele);
Lvector = zeros(nele,1);
Aredtop = zeros(nele,1);
Aredside = zeros(nele,1);

% assebly of matrices
for i=1:nele
    x1 = x(ID(i,1));
    x2 = x(ID(i,2));
    y1 = y(ID(i,1));
    y2 = y(ID(i,2));
    L0=sqrt((x2-x1)^2+(y2-y1)^2);
    Lvector(i,1)=L0;
    Biloc = [ -(x2-x1)  -(y2-y1)  (x2-x1)  (y2-y1)]/L0;
    ind = [2*ID(i,1)-1,2*ID(i,1),2*ID(i,2)-1,2*ID(i,2)];
    Qfull(ind,i) = sqrt(E*w/L0)*Biloc;
    Jfull(ind,i) = 0.5*rho*L0*w* [1  1  1  1];
    phi=asin((y2-y1)/L0); % vertical gives minimum
    Aredtop(i)=2*(rad2deg(abs(phi))*a+b)*lifereq;
    Aredside(i)=2*(b+90*a)*lifereq/w;
end


% Application of b.c.
Qred = Qfull(free_ind,:);
Jred = Jfull(free_ind,:);

Fred=F(free_ind,:);

% Maximum stiffness as parameter for scaling of the problem
A = 0.5*(Alow+Ahigh)*ones(nele,1); % average properties for height
K=Qfull*diag(A)*Qfull';
Kmax=0.1*max(max(K));

end