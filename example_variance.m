% Example of the usage of functions 
%
% find_ranges_exp(N0,r,n,k)
% pdf_t_exp(N0,r,n,k,vec_tim)
%
% for computing variance and std of the time T_k 
% in the coalescence tree with n leaves
% 2 <= k <= n
% for the model with exponential population growth N(t)=N0*exp(-rt)

clear all

% population size at present
N0=2.0e6;

% exponential growth coefficient
r=0.001;

% coalescence tree size (number of leaves)
n=800;

% index of coalescence time 
k=200;

% Find ranges
[min_t, max_t]= find_ranges_exp(N0,r,n,k);

% Define integrand for expecation
timinte=@(t)(pdf_t_exp(N0,r,n,k,t).*t);

% Compute expectation of T_k by numerical integration
[etint,errbnd] = quadgk(timinte,min_t,max_t,'RelTol',1e-8,'AbsTol',1e-11);

% Define integrand for variance
timintv=@(t)(pdf_t_exp(N0,r,n,k,t).*((t-etint).^2));

% Compute variance of T_k by numerical integration
[vtint,errbnd] = quadgk(timintv,min_t,max_t,'RelTol',1e-7,'AbsTol',1e-10);

% Results
disp(['Expectation of time T' num2str(k) ' in the coalescence tree = ' num2str(etint)])
disp(['Variance of time T' num2str(k) ' in the coalescence tree = ' num2str(vtint)])
disp(['Std of time T' num2str(k) ' in the coalescence tree = ' num2str(sqrt(vtint))])

