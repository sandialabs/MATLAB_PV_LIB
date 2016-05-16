function [outA, outxk, outY, kflag]=Schumaker_QSpline(X,Y)
% Schumaker_QSpline Fit a quadratric spline which preserves
% monotonicity and convexity in the data
%
% Syntax
%   [outA, outxk, outY, kflag] = Schumaker_QSpline(X,Y)
%
% Description
%   Calculates coefficients for C1 quadratic spline interpolating data X,Y,
%   where length(X) = N and length(Y) = N, which preserves monotonicity and
%   convexity in the data.
%
% Input:
%   X, Y - vectors of length N containing (X,Y) points between which the 
%      spline will interpolate.
%
% Output:
%   outA - a Nx3 matrix of coefficients where the ith row defines the 
%      quadratic interpolant between xk_i to xk_(i+1), i.e., 
%       Y=A(i,1)*(X-xk(i))^2 + A(i,2)*(X-xk(i)) + A(i,3)
%   outxk - an ordered vector of knots, i.e., values xk_i where the 
%      spline changes coefficients.  All values in X are used as knots.  
%      However, the algorithm may insert additional knots between data 
%      points in X where changes in convexity are indicated by the 
%      (numerical) derivative. Consequently output outxk has 
%      length >= length(X).
%   outY - y values corresponding to the knots in outxk.  Contains the
%      original data points Y, and also y-values estimated from the spline
%      at the inserted knots.
%   kflag - a vector of length(outxk) of logicals, which are set to true
%      for elements of outxk that are knots inserted by the algorithm.
% Sources:
%   [1] L. L. Schumaker, "On Shape Preserving Quadratic Spline 
%      Interpolation", SIAM Journal on Numerical Analysis 20(4), 
%      August 1983, pp 854-864
%   [2] M. H. Lam, "Monotone and Convex Quadratic Spline Interpolation",
%      Virginia Journal of Science 41(1), Spring 1990 

% ---------  a small number used to decide when a slope is equivalent to zero
eps=1e-6;  

% ---------  make input vectors into columns if necessary.
X=X(:);
Y=Y(:);

N=length(X);

% --------- compute various values used by the algorithm: differences, length of line
% segments between data points, and ratios of differences.
delX = diff(X);  % delX(i) = X(i+1)-X(i)
delY = diff(Y);

L=(delX.^2 + delY.^2).^0.5;
delta=delY./delX;

% -------- calculate first derivative at each x value per [2]

s=zeros(size(X));

left=[0; delta];
right=[delta; 0];

pdelta=left.*right;

u=pdelta>0;

% [2], Eq. 9 for interior points
s(u)=pdelta(u)./(0.5*left(u)+0.5*right(u));  % fix tuning parameters in [2], Eq. 9 at chi = 0.5 and eta = 0.5

% [2], Eq. 7 for left endpoint
if delta(1)*(2*delta(1)-s(2))>0
    s(1)=2*delta(1)-s(2);
end;

% [2] Eq. 8 for right endpoint
if delta(N-1)*(2*delta(N-1)-s(N-1))>0
    s(N)=2*delta(N-1)-s(N-1);
end;

% -------- determine knots.  Start with initial points x
 % [1], Algorithm 4.1 first 'if' condition of Step 5 defines intervals which won't get internal knots
tests = s(1:(N-1))+s(2:N);
u=abs(tests-2*delta(1:(N-1)))<=eps;
% u=true for an interval which will not get an internal knot

P=sum(u);

K=N+sum(~u);  % total number of knots = original data + inserted knots

% -------- set up output arrays
xk=zeros(K,1);  % knot locations, first N-1 and very last (N+K) are original data
yk=zeros(K,1);  % function values at knot locations
flag=false(K,1); % logicals that will indicate where additional knots are inserted
A=zeros(K,3);   % coefficients, e.g., A(1,1)*x^2+A(1,2)*x+A(1,3)

% --------- structures needed to compute coefficients, have to be maintained in association with each knot

tmpX=X(1:(N-1));
tmpY=Y(1:(N-1));
tmpX2=X(2:N);
tmpY2=Y(2:N);
tmps=s(1:(N-1));
tmps2=s(2:N);
diffs=diff(s);

% --------- structure to contain information associated with each knot, used to calculate coefficients
U = zeros(K,6); 

U(1:(N-1),:) = [tmpX tmpX2 tmpY tmps tmps2 delta];    

% --------- [1], Algorithm 4.1 subpart 1 of Step 5
xk(u)=tmpX(u);  % original X values that are left points of intervals without internal knots
yk(u)=tmpY(u);
A(u,3)=tmpY(u);   % constant term for each polynomial for intervals without knots
A(u,2)=s(u);
A(u,1)=0.5*diffs(u)./delX(u);  % leading coefficient

% --------- [1], Algorithm 4.1 subpart 2 of Step 5
xk(~u)=tmpX(~u); % original X values that are left points of intervals with internal knots
yk(~u)=tmpY(~u);

a=s(1:(N-1))-delta(1:(N-1));
b=s(2:N)-delta(1:(N-1));

sbar=zeros(K,1);
eta=sbar;
xi=zeros(K,1);  % will contain mapping from the left points of intervals containing an added knot to each interval's internal knot value

v=(~u)&(a.*b>=0); % first 'else' in Algorithm 4.1 Step 5
Q=sum(v); % number of this type of knot to add
if Q>0
    xk(N:(N+Q-1))=0.5*(tmpX(v)+tmpX2(v));  % knot location
    U(N:(N+Q-1),:)=[tmpX(v) tmpX2(v) tmpY(v) tmps(v) tmps2(v) delta(v)];
    xi(v)=xk(N:(N+Q-1));
end

w=~u&~v&(abs(a)>abs(b)); % second 'else' in Algorithm 4.1 Step 5
R=sum(w);
if R>0
    xk((N+Q):(N+Q+R-1))=tmpX2(w)+a(w).*delX(w)./diffs(w);
    U((N+Q):(N+Q+R-1),:)=[tmpX(w) tmpX2(w) tmpY(w) tmps(w) tmps2(w) delta(w)];
    xi(w)=xk((N+Q):(N+Q+R-1));
end

z=~u&~v&~w; % last 'else' in Algorithm 4.1 Step 5
S=sum(z);
if S>0
    xk((N+Q+R):(N+Q+R+S-1))=tmpX(z)+b(z).*delX(z)./diffs(z);
    U((N+Q+R):(N+Q+R+S-1),:)=[tmpX(z) tmpX2(z) tmpY(z) tmps(z) tmps2(z) delta(z)];
    xi(z)= xk((N+Q+R):(N+Q+R+S-1));
end

% define polynomial coefficients for intervals with added knots
ff=~u;
sbar(ff)=(2*U(ff,6)-U(ff,5))+(U(ff,5)-U(ff,4)).*(xi(ff)-U(ff,1))./(U(ff,2)-U(ff,1));
eta(ff)=(sbar(ff)-U(ff,4))./(xi(ff)-U(ff,1));

ff=(N:(N+Q+R+S-1));
sbar(ff)=(2*U(ff,6)-U(ff,5))+(U(ff,5)-U(ff,4)).*(xk(ff)-U(ff,1))./(U(ff,2)-U(ff,1));
eta(ff)=(sbar(ff)-U(ff,4))./(xk(ff)-U(ff,1));

A(~u,3)=U(~u,3);  % constant term for polynomial for intervals with internal knots
A(~u,2)=U(~u,4);
A(~u,1)=0.5*eta(~u);  % leading coefficient

A((N):(N+Q+R+S-1),3)=U((N):(N+Q+R+S-1),3)+U((N):(N+Q+R+S-1),4).*(xk((N):(N+Q+R+S-1))-U((N):(N+Q+R+S-1),1))+...
    0.5*eta((N):(N+Q+R+S-1)).*(xk((N):(N+Q+R+S-1))-U((N):(N+Q+R+S-1),1)).^2;
A((N):(N+Q+R+S-1),2)=sbar((N):(N+Q+R+S-1));
A((N):(N+Q+R+S-1),1)=0.5*(U((N):(N+Q+R+S-1),5)-sbar((N):(N+Q+R+S-1)))./...
    (U((N):(N+Q+R+S-1),2)-U((N):(N+Q+R+S-1),1));

yk(N:(N+Q+R+S-1))=A((N):(N+Q+R+S-1),3);

xk(N+Q+R+S)=X(N);
yk(N+Q+R+S)=Y(N);
flag(N:(N+Q+R+S-1))=true;  % these are all inserted knots

tmp=[xk A yk flag];
tmp2=sortrows(tmp,1); % sort output in terms of increasing x (original plus added knots)
outxk=tmp2(:,1);
outN=length(outxk);
outA=tmp2(1:(outN-1),2:4);
outY=tmp2(:,5);
kflag=tmp2(:,6);