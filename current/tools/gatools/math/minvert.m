   function [matinv,minsig,numsig] = minvert(mat,tol,figure_num)
 %
%  SYNTAX:  [matinv,minsig,numsig] = minvert(mat,tol,figure_num)
%
%  PURPOSE:  Perform SVD generalized inverse of non-square or square matrix.
%
%  INPUTS:
%     mat = matrix to invert (can be square or non-square)
%     tol = tolerance for inversion (retains singular values st sv/sv_max>tol)
%		If tol>1, then interprets as # of singular values to keep,,
%		and returns sig(tol)/max(sig)
%     figure_num = specifies figure number of plotting desired (default=no plot)
%
%  OUTPUTS:
%     matinv = SVD generalized inverse of mat.
%     minsig = minimum singular value retained divided by max singular value
%     numsig = number of singular values retained
%
%  RESTRICTIONS:
%     None.
%
%  METHOD: 
%     SVD generalized inverse. Adapted from IDL function minvert.pro.

%  WRITTEN BY:  Dave Humphreys 	ON	4/30/93
%
%  MODIFICATION HISTORY:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)minvert.m	1.2 03/07/12

% Prelims:
   mu0 = 0.4*pi;
   if nargin<2, tol=0; end

% Test special conditions:
   ss=size(mat);
   if(ss(1)==1|ss(2)==1)
      disp('WARNING: Matrix is a vector; returning vector inverse.')
      matinv = mat/norm(mat);
      minsig = 1;
      return
   end

% SVD and truncate:
   [U,sig,V]=svd(mat);
   sigv = diag(sig);
   s=sigv/max(sigv);
   ns = length(s);
   if tol>1
     itmp = (1:min(tol,ns))';   %interpret tol as max idx sv to keep 
   else
     itmp = find(s > tol);  %else interpret tol as tolerance 
   end  %end if tol>1
   sigvi = zeros(ns,1);
   sigvi(itmp) = 1./sigv(itmp);
   minsig = s(max(itmp));
   numsig = max(itmp);

% Make inverse:
   sigi = zeros(size(sig))';
   sigi(1:ns,1:ns) = diag(sigvi);
   matinv = V*sigi*U';

% Plot sv's if figure_num exists:
   if (exist('figure_num','var'))
     figure(figure_num),clf
     plot(s)
     hold on
     plot(s(itmp),'x')
     xlabel('SV idx')
     ylabel('SV/max(SV)')
   end
                                               
