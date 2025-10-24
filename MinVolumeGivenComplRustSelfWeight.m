function [Afinal,tau,tot_it,tot_time] = MinVolumeGivenComplRustSelfWeight(nele,nfree,taumax,Ascaling,Alow,Ahigh,Kmax,Qred,Jred,Fred,Lvector,Aredtop,Aredside)
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
g=9.81e-4; %10 kN/kg

% 
cvx_begin
  variable A(nele) 
   variable tau
   minimize( A'*Lvector ) 
   subject to
   A<=Ahigh*Ascaling
   A>=Alow*Ascaling
   tau <= taumax
   %[tau Fred'; Fred (Qred*diag((A/Ascaling-Aredtop-A.*Aredside/Ascaling+Aredtop.*A.*Aredside/Ascaling))*Qred')]/Kmax == semidefinite(nfree+1)
   [tau (Fred+Jred*A*g)'/sqrt(Kmax); (Fred+Jred*A*g)/sqrt(Kmax) (Qred*diag((A/Ascaling-Aredtop-A.*Aredside/Ascaling+Aredtop.*Aredside))*Qred')/Kmax] == semidefinite(nfree+1)
   %[tau Fred'/sqrt(Kmax); Fred/sqrt(Kmax) (Qred*diag((A/Ascaling-Aredtop-A.*Aredside/Ascaling+Aredtop.*A.*Aredside/Ascaling))*Qred')/Kmax] == semidefinite(nfree+1)
   %[tau Fred'; Fred (Qred*diag((A/Ascaling-Aredtop-A.*Aredside/Ascaling+Aredtop.*A.*Aredside/Ascaling))*Qred')] == semidefinite(nfree+1)
cvx_end

 tot_it = tot_it + cvx_slvitr;
 tot_time = tot_time + cvx_cputime;



Afinal =A/Ascaling;
end
