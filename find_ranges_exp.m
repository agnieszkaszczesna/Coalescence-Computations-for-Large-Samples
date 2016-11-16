% FUNCTION: 
% find_ranges_exp

% PURPOSE: 
% compute range (support) of the probability distribution in the coalescence tree
% for time T_k in the coalescence tree with n leaves
% 2 <= k <= n
% for the model with exponential population growth N(t)=N0*exp(-rt)

% INPUT: 
% N0 - population size at present 
% r - exponential growth coefficient
% n - coalescence tree size (number of leaves)
% k - index of time to coalescence event 

% INTERNAL PARAMETER - threshold for the value of the pdf at the tail and
% the reference value: 
% D=1.0e-5;

% METHOD: 
% heurisitc method - sampling values of the pdf until threshold condition
% is satisfied
% based on the paper 
% "Coalescence computations for large samples drawn from populations of time-varying sizes" 
% Andrzej Polanski, Agnieszka Szczesna, Mateusz Garbulowski, Marek Kimmel
% submitted to PLOS ONE

% OUTPUT:
% [min_t, max_t] - minimal and maximal time  

% WARNING
% This function does not have control check of the input data. Providing wrong
% data will lead to unpredictable behavior

% AUTHOR:
% Andrzej Polanski
% email: andrzej.polanski@polsl.pl

function [min_t, max_t]= find_ranges_exp(N0,r,n,k)
  % product parameter
  rho=N0*r;

  % Expected canonical time
  TRUEEXPT=2*(1/(k-1)-1/n);
  % Variance of the canonical time
  VS1=(k:n).^-2;
  VS2=((k-1):(n-1)).^-2;
  TRUEVART=4*(VS1*VS2');

  % reference value of the pdf
  % 
  t=TRUEEXPT;
  q=pdf_t_cc(n,k,t);
  vpr=max(0,q)*(1+rho*t)/N0;

  % INTERNAL PARAMETER - threshold for the value of the pdf at the tail and
  % the reference value: 
  D=1.0e-5;

  % sample probability distribution left tail
  vt=TRUEEXPT;
  st=sqrt(TRUEVART)/4;
  t=vt;
  while 1
     t=t-st;
     if t <= 0
        t=0;
        break
     end
     % canonical coalescence pdf
     q=pdf_t_cc(n,k,t);
     % scale pdf to exponential scenario
     v_p=max(0,q)*(1+rho*t)/N0;
     if (v_p/vpr) < D
         break
     end
  end
  % scale time
  if r==0
       min_t=N0*t;
  else
       min_t=(1/r)*log(1+rho*t);
  end

  % sample probability distribution right tail
  vt=TRUEEXPT;
  st=sqrt(TRUEVART);
  t=vt; 
  while 1
     t=t+st;
     % canonical coalescence pdf
     q=pdf_t_cc(n,k,t);
     % scale pdf to exponential scenario
     v_p=max(0,q)*(1+rho*t)/N0;
     if (v_p/vpr) < D
         break
     end
  end
  % scale time
  t=1.5*t;
  if r==0
        max_t=N0*t;
  else
        max_t=(1/r)*log(1+rho*t);
  end
end



