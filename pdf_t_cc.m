% FUNCTION: 
% pdf_t_cc
% computing probability density function of times in the canonical coalescence tree with
% constant population size N0=1

% PURPOSE: 
% computes probability distribution of time T_k in the canonical coalescence tree with
% constant population size N0=1
% 2 <= k <= n

% INPUT: 
% n - coalescence tree size (number of leaves)
% k - index of coalescence time  
% tcc - value of the canonical time


% METHOD: 
% the algorithm and notation are described in the paper 
% "Coalescence computations for large samples drawn from populations of time-varying sizes" 
% Andrzej Polanski, Agnieszka Szczesna, Mateusz Garbulowski, Marek Kimmel
% submitted to PLOS ONE

% OUTPUT:
% q - value of probability density function corresponding to tcc

% REMARKS:
% Tested in Matlab 2014 and 2015 

% WARNING
% This function does not have control check of the input data. Providing wrong
% data will result in an unpredictable behavior

% AUTHOR:
% Andrzej Polanski
% email: andrzej.polanski@polsl.pl

function q=pdf_t_cc(n,k,tcc)
   % scale tcc
   tcc=tcc*(k*k);
   % computecoeffficients for transformed pdf
   numvec=zeros(n-k+1,1)';
   denvec=zeros(n-k+1,1)';
   denvec1=ones(n-k+1,1)';
   for kn=k:n
      numvec(kn-k+1)=(nchoosek(kn,2))/(k*k);
      denvec(kn-k+1)=(nchoosek(kn,2))/(k*k);
   end
     
   % compute probability density
   i=sqrt(-1);
   if tcc==0
      if (n-k) <= 3           % compute pdf by partial fraction expansion (speeds up computations)
          q=part_fra_dens(n,k,tcc);
      else                    % 
          q=0;
      end
   else
      if (n-k) <= 3            % compute pdf by partial fraction expansion (speeds up computations)
          q=part_fra_dens(n,k,tcc);
      else                     % compute pdf by integral formula (split integration range into fragments)
          omg_prod=@(omg)((real(prod((ones(length(omg),1)*numvec)./(i*(omg')*denvec1+ones(length(omg),1)*denvec),2)')).*cos(omg*tcc));
          limlo=0;
          limup=pi/(2*tcc);
          q=0;
          while 1
             [qint,errbnd] = quadgk(omg_prod,limlo,limup,'RelTol',1e-8,'AbsTol',1e-12);
              q=q+qint;
              limlo=limup;
              limup=limup+(4*pi)/tcc;
              if abs(qint)<1.0e-11
                 break
              end
          end
          % scale integral
          q=2*q/pi;
      end
   end
   % rescale pdf to tcc
   q=q*(k*k);
end
   
   
% auxiliary FUNCTION: 
% part_fra_dens
% probability distribution of times in the canonical coalescence tree with
% constant population size N0=1
% compute coefficients of tha partial fraction expansion (Polanski,
% Bobrowski, Kimmeel, TPB, 2003)

% WARNING
% This function does not have control check of the input data. Providing wrong
% data will result in an unpredictable behavior

% AUTHOR:
% Andrzej Polanski
% email: andrzej.polanski@polsl.pl
  
function q=part_fra_dens(n,k,tcc)
   
   if k==n
       a_v=1;
   else
      vec_coef=zeros(n-k+1,1);
      for kk=k:n
          vec_coef(kk-k+1)=nchoosek(kk,2)/(k*k);
      end
      a_v=zeros(n-k+1,1);
      onvec=ones(n-k,1);
      for jj=k:n
          prd_vec_num=[vec_coef(1:jj-k); vec_coef(jj-k+2:n-k+1)];
          prd_vec_den=prd_vec_num-onvec*nchoosek(jj,2)/(k*k);
          a_v(jj-k+1)=prod(prd_vec_num./prd_vec_den);
      end
   end
   % compute pdf
   q=0;
   for jj=k:n
       q=q+a_v(jj-k+1)*nchoosek(jj,2)*exp(-tcc*nchoosek(jj,2)/(k*k))/(k*k);
   end
end