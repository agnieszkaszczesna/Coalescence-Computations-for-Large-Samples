% FUNCTION: 
% pdf_t_exp
% probability dernsity function of times in the coalescence tree with
% exponential history of the size change of the population

% PURPOSE: 
% computes probability distribution of time T_k in the coalescence tree with n leaves
% 2 <= k <= n
% for the model with exponential population growth N(t)=N0*exp(-rt)
 
% INPUT: 
% N0 - population size at present 
% r - exponential growth coefficient
% n - coalescence tree size (number of leaves)
% k - index of the coalescence time 
% vec_tim - vector of times


% METHOD: 
% the algorithm and notation are described in the paper 
% "Coalescence computations for large samples drawn from populations of time-varying sizes" 
% Andrzej Polanski, Agnieszka Szczesna, Mateusz Garbulowski, Marek Kimmel
% submitted to PLOS ONE

% OUTPUT:
% vec_p - vector of values of probability density function corresponding to vec_tim

% REMARKS:
% Tested in Matlab 2014 and 2015 

% WARNING
% This function does not have control check of the input data. Providing wrong
% data will result in an unpredictable behavior

% AUTHOR:
% Andrzej Polanski
% email: andrzej.polanski@polsl.pl

function vec_p=pdf_t_exp(N0,r,n,k,vec_tim)

   %%%%%%%%%%%%%%%%% 
   tN0vec=vec_tim/N0;
   rho=N0*r;

   % output buffer
   vec_p=0*vec_tim;

   % scale time
   if r==0
       vec_tau=tN0vec;
   else
       vec_tau=(1/rho)*(exp(rho*tN0vec)-1);
   end
      
   % compute probability densities
   for kk=1:length(vec_tau)
      t=vec_tau(kk);
      q=pdf_t_cc(n,k,t);
      % scale pdf
      vec_p(kk)=max(0,q)*(1+rho*t)/N0;
   end
end


