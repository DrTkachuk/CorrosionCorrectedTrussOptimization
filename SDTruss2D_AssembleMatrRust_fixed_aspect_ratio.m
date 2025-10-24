function [Qred,Jred,Kmax,nfree,Fred,Lvector,Aredtop,Aredside,Alowvec,Ahighvec] = SDTruss2D_AssembleMatrRust_fixed_aspect_ratio(E,rho,x,y,ndof,nele,ID,ind_DBC,F,Alow,Ahigh,a,b,lifereq,alpha)

full_ind = [1:ndof]
free_ind = setdiff(full_ind,ind_DBC)
nfree=size(free_ind,2)


Qfull = zeros(ndof,nele);
Jfull = zeros(ndof,nele);
Lvector = zeros(nele,1);
Aredtop = zeros(nele,1);
Aredside = zeros(nele,1);
Alowvec = zeros(nele,1);
Ahighvec = zeros(nele,1);

for i=1:nele
    x1 = x(ID(i,1));
    x2 = x(ID(i,2));
    y1 = y(ID(i,1));
    y2 = y(ID(i,2));
    L0=sqrt((x2-x1)^2+(y2-y1)^2);
    Lvector(i,1)=L0;
    Biloc = [ -(x2-x1)  -(y2-y1)  (x2-x1)  (y2-y1)]/L0;
    ind = [2*ID(i,1)-1,2*ID(i,1),2*ID(i,2)-1,2*ID(i,2)];
    Qfull(ind,i) = sqrt(E/L0)*Biloc;
    Jfull(ind,i) = 0.5*rho*L0* [1  1  1  1];
    phi=asin((y2-y1)/L0); % vertical gives minimum
    Aredtop(i)=2*(rad2deg(abs(phi))*a+b)*lifereq;
    Aredside(i)=2*(b+90*a)*lifereq;
    thigh=sqrt(Ahigh*alpha);
    whigh=sqrt(Ahigh/alpha);
    Ahighvec(i)=(thigh-Aredtop(i))*(whigh-Aredside(i));
    tlow=sqrt(Alow*alpha);
    wlow=sqrt(Alow/alpha);
    Alowvec(i)=(tlow-Aredtop(i))*(wlow-Aredside(i));
end

A = 0.5*(Alow+Ahigh)*ones(nele,1)


K=Qfull*diag(A)*Qfull'
M=diag(Jfull*A) 

Qred = Qfull(free_ind,:)
Jred = Jfull(free_ind,:)


Fred=F(free_ind)

Kmax=0.1*max(max(K))

end