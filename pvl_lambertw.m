function w = pvl_lambertw(z)
%
%
% PVL_LAMBERTW Compute values for the Lambert W Function W(z).
%
% Syntax
%   W = pvl_lambertw(z)
% 
% Description
%   W = LAMBERTW(Z) computes the principal value of the Lambert 
%   W Function, the solution of Z = W*exp(W).  Z may be a 
%   complex scalar or array.  For real Z, the result is real on
%   the principal branch for Z >= -1/e.
%
%   The algorithm uses series approximations as initializations
%   and Halley's method as developed in Corless, Gonnet, Hare,
%   Jeffrey, Knuth, "On the Lambert W Function", Advances in
%   Computational Mathematics, volume 5, 1996, pp. 329-359.

%   Original code by Pascal Getreuer 2005-2006, modified by 
%   Didier Clamond, 2005.  Code downloaded from 
%   http://www.getreuer.info/home/lambertw and modified for inclusion in
%   PVLib.
%
%   Matlab includes a lambertw.m function using a very similar algorithm
%   in the Symbolic Math Toolbox.
%
% Inputs: 
%       Z - scalar or vector of values at which W(Z) will be evaluated.
%
% Outputs:
%       w - scalar or vector of values of W(Z) on the principal branch.
%
% Sources:
%
%   [1] R.M. Corless, G.H. Gonnet, D.E.G. Hare, G.J. Jeffery, and D.E. 
%       Knuth. "On the Lambert W Function." Advances in Computational 
%       Mathematics, vol. 5, 1996

% Use asymptotic expansion w = log(z) - log(log(z)) for most z
tmp = log(z + (z == 0));
w = tmp - log(tmp + (tmp == 0));

% Use a series expansion when close to the branch point -1/e
k = (abs(z + 0.3678794411714423216) <= 1.5);
tmp = sqrt(5.43656365691809047*z + 2) - 1; % [1], Eq. 4.22 and text
w(k) = tmp(k);

for k = 1:36
   % Converge with Halley's method ([1], Eq. 5.9), about 5 iterations satisfies
   % the tolerance for most z
   c1 = exp(w);
   c2 = w.*c1 - z;
   w1 = w + (w ~= -1);
   dw = c2./(c1.*w1 - ((w + 2).*c2./(2*w1)));
   w = w - dw;

   if all(abs(dw) < 0.7e-16*(2+abs(w)))
      break;
   end
end

