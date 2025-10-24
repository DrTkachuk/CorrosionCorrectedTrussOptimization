function [Arfinal,tau,tot_it,tot_time,ArNormEvol] = MinVolumeGivenComplAndFreq_fixed_aspect_ratio(nele,nfree,taumax,Ascaling,Alowvec,Ahighvec,Kmax,Qred,Fred,Lvector,Aredtop,Aredside,alpha,niter)
% minimizes volume for given compliance and  frequency constrainta and
% known rust rate

% Anton Tkachuk anton.tkachuk@kau.se
% 02.05.2022
% ToDos
%comment
%add flag and check for base frequency, done
%add M0, done
%add more checks, done
%think about data structures (nele,nfree,Qred,Jred,M0 and Fred are
%structural properties) :: next paper (NP)
%think about structure parameters ::(NP)
%extend to multiple load cases ::(NP)
%exit condition ::(NP)




% setting properties of the cvx solver
%cvx_solver_settings('cvx_precision',CVXPRECISION);
cvx_solver_settings( 'cvx_slvitr', 120);
%cvx_solver_settings('cvx_solver',CVXSOLVER );
 cvx_precision medium
 cvx_solver sdpt3
%cvx_solver sedumi
% cvx_solver_settings( 'cvx_slvitr', 120)
tot_it = 0; % accumulated iteration and cpu time for semi-definite programs
tot_time = 0;

obj_scaling=10*norm(Lvector,1)*min(Alowvec);

Weights=eye(nele);

Arold=zeros(nele,1);
ArNormEvol=zeros(niter,1);

for nn=1:niter % conncave-convex iterations
% beginning of definition of SDP
cvx_begin 
  variable Ar(nele)   %% design variables
   variable tau       %% auxiliary variable for compliance constraint
   minimize( Ar'*Weights*Lvector/obj_scaling ) 
   subject to
   Ar<=Ahighvec*Ascaling  %% upper limit
   Ar>=Alowvec*Ascaling   %% lower limit
   tau <= taumax
   % two versions of compliance constraint with scaling
   %[tau Fred'; Fred (Qred*diag((Ar/Ascaling))*Qred')]/Kmax == semidefinite(nfree+1)
   [tau Fred'/sqrt(Kmax); Fred/sqrt(Kmax) (Qred*diag((Ar/Ascaling))*Qred')/Kmax] == semidefinite(nfree+1)  
cvx_end

 tot_it = tot_it + cvx_slvitr;
 tot_time = tot_time + cvx_cputime;
% update of weights and writing logs
Weights=eye(nele)+0.5*diag((Aredside/sqrt(alpha)+Aredtop*sqrt(alpha))./sqrt(Ar/Ascaling));
ArNormEvol(nn)=norm(Ar/Ascaling-Arold,2)
fprintf("Change of norm %f",ArNormEvol(nn));
Arold=Ar/Ascaling;
end
Arfinal =Ar/Ascaling;
end

%   (Qred*diag(A/Ascaling-Aredtop-A.*Aredside/Ascaling)*Qred'- omega2base*(diag(Jred*A/Ascaling) + M0))/Kmax == semidefinite(nfree)
   %[tau Fred'/sqrt(Kmax); Fred/sqrt(Kmax) (Qred*diag((A/Ascaling-Aredtop-A.*Aredside/Ascaling))*Qred')/Kmax] == semidefinite(nfree+1)
   %   [tau Fred'; Fred (Qred*diag((Ar/Ascaling))*Qred')]/Kmax == semidefinite(nfree+1)