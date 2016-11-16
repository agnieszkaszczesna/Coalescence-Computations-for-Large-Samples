SOFTWARE FOR COMPUTING PROBABILITY DISTRIBUTIONS OF TIMES IN THE COALESCENCE 
PROCESS.

Programmed on the basis of algorithms described in the paper  
"Coalescence computations for large samples drawn from populations of time-varying sizes" 
Andrzej Polanski, Agnieszka Szczesna, Mateusz Garbulowski, Marek Kimmel
submitted to PLOS ONE

Contact
andrzej.polanski@polsl.pl
agnieszka.szczesna@polsl.pl

Contain following functions and scripts:

% FUNCTION: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find_ranges_exp
% PURPOSE: 
% compute range (support) of the probability distribution in the coalescence tree
% for time T_k in the coalescence tree with n leaves
% 2 <= k <= n
% for the model with exponential population growth N(t)=N0*exp(-rt)


% FUNCTION: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pdf_t_cc
% PURPOSE: 
% computes probability distribution of time T_k in the canonical coalescence tree with
% constant population size N0=1
% 2 <= k <= n


% FUNCTION: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pdf_t_exp
% PURPOSE: 
% computes probability distribution of time T_k in the coalescence tree with n leaves
% 2 <= k <= n
% for the model with exponential population growth N(t)=N0*exp(-rt)


% SCRIPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example_pdf 
% Example of the usage of functions 
% find_ranges_exp(N0,r,n,k)
% pdf_t_exp(N0,r,n,k,vec_tim)
%
% for computing and plotting probability density function of the time T_k 
% in the coalescence tree with n leaves
% 2 <= k <= n
% for the model with exponential population growth N(t)=N0*exp(-rt)


% SCRIPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example_expectation 
% Example of the usage of functions 
% find_ranges_exp(N0,r,n,k)
% pdf_t_exp(N0,r,n,k,vec_tim)
%
% for computing expectation of the time T_k 
% in the coalescence tree with n leaves
% 2 <= k <= n
% for the model with exponential population growth N(t)=N0*exp(-rt)


% SCRIPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example_variance
% Example of the usage of functions 
% find_ranges_exp(N0,r,n,k)
% pdf_t_exp(N0,r,n,k,vec_tim)
%
% for computing variance and std of the time T_k 
% in the coalescence tree with n leaves
% 2 <= k <= n
% for the model with exponential population growth N(t)=N0*exp(-rt)



Computing probability distributions of times in the coalescence tree for exponential scenario
of population  growth is based on the direct method of numerical integration. 

We have tested these functions for the range of values of the product parameter
$0 <= \rho <= 10^6$ and genealogy sizes $n <= 10^4$.
In the provided programs all parameters are set automatically and the 
relative errors are $10^{-5}$ or better.  






