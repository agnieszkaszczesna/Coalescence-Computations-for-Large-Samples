% Example of the usage of functions 
%
% find_ranges_exp(N0,r,n,k)
% pdf_t_exp(N0,r,n,k,vec_tim)
%
% for computing and plotting probability density function of the time T_k 
% in the coalescence tree with n leaves
% 2 <= k <= n
% for the model with exponential population growth N(t)=N0*exp(-rt)

clear all

% population size at present
N0=2.0e6;

% exponential growth coefficient
r=0.1;

% coalescence tree size (number of leaves)
n=800;

% index of coalescence time 
k=797;

% Find ranges
[min_t, max_t]= find_ranges_exp(N0,r,n,k);

% Set time axis
dt=(max_t-min_t)/100;
vec_tim=min_t:dt:max_t;

% Compute pdf function of T_k
vec_p=pdf_t_exp(N0,r,n,k,vec_tim);

% Plot
plot(vec_tim,vec_p)
grid on
xlabel('Time [number of generations]')
ylabel('pdf')
title(['Pdf of time T' num2str(k) ' in the coalescence tree'])
