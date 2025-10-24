function [Afinal,tau,tot_it,tot_time] = MinVolumeGivenComplRustMultLoad(nele,nfree,taumax,Ascaling,Alow,Ahigh,Kmax,Qred,nlc,Fred,Lvector,Aredtop,Aredside)
% minimizes volume for given compliance and  frequency constrainta and
% known rust rate
% Input parameters::
% nele - number of members, each member has one design parameter height A(i)
% nfree - number of free degrees of freedom, consides with row dimension of Qred and Fred
% taumax - maximum admissible compliance, the main constraint of optimization
% Ascaling - scaling of the design variable
% Alow,Ahigh - lower and upper limit for design parameters, height
% Kmax - scaling of the stiffness matrix 
% Qred - matrix containing a scaled B-operators in rows
% nlc - number of load cases
% Fred - matrix of load cases (nfree x nlc)
% Lvector - vector with member length values
% Aredtop - reduction of thickness in height direction
% Aredside - scaled reduction of thicknes in width direction
% Return parameters::
% Afinal - the final distribution of heights
% tau - achieved compliance
% tot_it - total number of required iteration of the selected SDP solver
% tot_time - time of simulation


% Anton Tkachuk anton.tkachuk@kau.se
% 25.06.2024
% ToDos
%comment
%add more checks
%extend to multiple load cases ::(NP)
%exit condition ::(NP)




% setting properties of the cvx solver
%cvx_solver_settings('cvx_precision',CVXPRECISION);
cvx_solver_settings( 'cvx_slvitr', 120);
%cvx_solver_settings('cvx_solver',CVXSOLVER );
%cvx_precision low

%cvx_solver_settings('sdpt3_options', struct('direc', 'NT'))

%cvx_precision high
cvx_precision medium
cvx_solver mosek
 % cvx_solver sdpt3
%cvx_solver sedumi
% cvx_solver_settings( 'cvx_slvitr', 120)
tot_it = 0; % accumulated iteration and cpu time for semi-definite programs
tot_time = 0;

obj_scaling=10*norm(Lvector,1)*Alow;

% 
cvx_begin
  variable A(nele) 
   variable tau
   minimize( A'*Lvector/obj_scaling ) 
   subject to
   A<=Ahigh*Ascaling
   A>=Alow*Ascaling
   tau <= taumax
   %[tau Fred'; Fred (Qred*diag((A/Ascaling-Aredtop-A.*Aredside/Ascaling+Aredtop.*A.*Aredside/Ascaling))*Qred')]/Kmax == semidefinite(nfree+1)
   %[tau*eye(nlc) Fred'/sqrt(Kmax); Fred/sqrt(Kmax) (Qred*diag((A/Ascaling-Aredtop-A.*Aredside/Ascaling+Aredtop.*A.*Aredside/Ascaling))*Qred')/Kmax] == semidefinite(nfree+nlc)
   [tau*eye(nlc) Fred'/sqrt(Kmax); Fred/sqrt(Kmax) (Qred*diag((A/Ascaling-Aredtop-A.*Aredside/Ascaling+Aredtop.*Aredside/Ascaling))*Qred')/Kmax] == semidefinite(nfree+nlc)   
   %[tau Fred'; Fred (Qred*diag((A/Ascaling-Aredtop-A.*Aredside/Ascaling+Aredtop.*A.*Aredside/Ascaling))*Qred')] == semidefinite(nfree+1)
cvx_end

 tot_it = tot_it + cvx_slvitr;
 tot_time = tot_time + cvx_cputime;


% rescaling result
Afinal =A/Ascaling;
end
