function [df, d2f]=numdiff(x,f)

% NUMDIFF : compute first and second order derivative using possibly unequally spaced
% data.
%
% Syntax
%   [df d2f] = numdiff(x, f)
%
% Description
%   numdiff computes first and second order derivatives using a 5th order
%   formula that accounts for possibly unequally spaced data.  Because a
%   5th order centered difference formula is used, numdiff returns NaNs for
%   the first 2 and last 2 points in the input vector for x.
%
% Input :
%   x - a vector of values of x
%   f - a vector of values of the function f for which derivatives are to 
%      be computed.  Must be the same length as x.
%
% Output:
%   df - a vector of length(x) containing the first derivative of f at 
%      each point x except at the first 2 and last 2 points.
%   d2f - a vector of length(x) containing the second derivative of f at 
%      each point x except at the first 2 and last 2 points.
%
% Sources:
% [1] M. K. Bowen, R. Smith, "Derivative formulae and errors for
% non-uniformly spaced points", Proceedings of the Royal Society A, 
% vol. 461 pp 1975-1997, July 2005. DOI: 10.1098/rpsa.2004.1430

% convert to column vectors if needed
f=f(:);
x=x(:);
N=length(f);

df = NaN(size(f));

% first two points are special

df(1:2) = NaN;

% last two points are special

df((N-1):N) = NaN;

% rest of points.  Take reference point to be the middle of each group of 5
% points.  Calculate displacements

F = [f(1:(N-4)) f(2:(N-3)) f(3:(N-2)) f(4:(N-1)) f(5:N)];

A = [x(1:(N-4)) x(2:(N-3)) x(3:(N-2)) x(4:(N-1)) x(5:N)] - repmat(x(3:(N-2)),1,5);

% A(:,3) is always zero but we retain it in case we want to generalize this
% function in the future

U(:,1) = A(:,2).*A(:,3).*A(:,4) + A(:,2).*A(:,3).*A(:,5) + A(:,2).*A(:,4).*A(:,5) + A(:,3).*A(:,4).*A(:,5);
U(:,2) = A(:,1).*A(:,3).*A(:,4) + A(:,1).*A(:,3).*A(:,5) + A(:,1).*A(:,4).*A(:,5) + A(:,3).*A(:,4).*A(:,5);
U(:,3) = A(:,1).*A(:,2).*A(:,4) + A(:,1).*A(:,2).*A(:,5) + A(:,1).*A(:,4).*A(:,5) + A(:,2).*A(:,4).*A(:,5);
U(:,4) = A(:,1).*A(:,2).*A(:,3) + A(:,1).*A(:,2).*A(:,5) + A(:,1).*A(:,3).*A(:,5) + A(:,2).*A(:,3).*A(:,5);
U(:,5) = A(:,1).*A(:,2).*A(:,3) + A(:,1).*A(:,2).*A(:,4) + A(:,1).*A(:,3).*A(:,4) + A(:,2).*A(:,3).*A(:,4);

L(:,1) = (A(:,1) - A(:,2)).*(A(:,1) - A(:,3)).*(A(:,1) - A(:,4)).*(A(:,1) - A(:,5));
L(:,2) = (A(:,2) - A(:,1)).*(A(:,2) - A(:,3)).*(A(:,2) - A(:,4)).*(A(:,2) - A(:,5));
L(:,3) = (A(:,3) - A(:,1)).*(A(:,3) - A(:,2)).*(A(:,3) - A(:,4)).*(A(:,3) - A(:,5));
L(:,4) = (A(:,4) - A(:,1)).*(A(:,4) - A(:,2)).*(A(:,4) - A(:,3)).*(A(:,4) - A(:,5));
L(:,5) = (A(:,5) - A(:,1)).*(A(:,5) - A(:,2)).*(A(:,5) - A(:,3)).*(A(:,5) - A(:,4));

df(3:(N-2)) = sum(-(U./L).*F,2);

% second derivative

U2(:,1) = A(:,2).*A(:,3)+A(:,2).*A(:,4)+A(:,2).*A(:,5)+A(:,3).*A(:,4)+A(:,3).*A(:,5)+A(:,4).*A(:,5);
U2(:,2) = A(:,1).*A(:,3)+A(:,1).*A(:,4)+A(:,1).*A(:,5)+A(:,3).*A(:,4)+A(:,3).*A(:,5)+A(:,4).*A(:,5);
U2(:,3) = A(:,1).*A(:,2)+A(:,1).*A(:,4)+A(:,1).*A(:,5)+A(:,2).*A(:,4)+A(:,2).*A(:,4)+A(:,4).*A(:,5);
U2(:,4) = A(:,1).*A(:,2)+A(:,1).*A(:,3)+A(:,1).*A(:,5)+A(:,2).*A(:,3)+A(:,2).*A(:,5)+A(:,3).*A(:,5);
U2(:,5) = A(:,1).*A(:,2)+A(:,1).*A(:,3)+A(:,1).*A(:,4)+A(:,2).*A(:,3)+A(:,2).*A(:,5)+A(:,3).*A(:,4);

d2f(3:(N-2)) = 2*sum(U2.*F,2);
