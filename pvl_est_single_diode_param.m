function [Io, Iph, Rs, Rsh, n] = pvl_est_single_diode_param(I, V, NsVth)
% PVL_EST_SINGLE_DIODE_PARAM regression technique to fit the single diode 
%   equation to data for a single IV curve
%
% Syntax
%   [Io Iph Rs Rsh n] = pvl_est_single_diode_param(I, V, NsVth)
%
% Description
%   pvl_est_single_diode_param uses a regression technique based on [1] to
%   fit the single diode equation to data for a single IV curve.  Although
%   values for each of the five parameters are returned, testing has shown  
%   only Rsh to be stable.  The other parameters, Rs, Io and n may be 
%   negative or imaginary even for IV curve data without obvious flaws.
%   Method coded here uses a principal component transformation of (V,I)
%   prior to regression to attempt to overcome effects of strong
%   colinearity between V and I over much of the I-V curve.
%
% Input:
%   I - a vector of length N of current for the IV curve.  The first value
%   is taken as Isc, the last value must be 0.
%   V - a vector of length N of voltage for the IV curve corresponding to
%   the current values in the input vector I.  The first value must be 0,
%   the last value is taken as Voc.
%   NsVth - the thermal voltage for the module, equal to Ns (number of
%   cells in series) x Vth (thermal voltage per cell).
%
% Output:
%   Io - the dark current value (A) for the IV curve
%   Iph - the light current value (A) for the IV curve
%   Rs - series resistance (ohm) for the IV curve
%   Rsh - shunt resistance (ohm) for the IV curve
%   n - diode (ideality) factor (unitless) for the IV curve
%
% Sources:
%   [1] A. Ortiz-Conde, F. Garci'a Sa'nchez, J. Murci, "New method to
%   extract the model parameters of solar cells from the explict analytic
%   solutions of their illuminated I-V characteristics", Solar Energy 
%   Materials and Solar Cells 90, pp 352-361, 2006. 
%   [2] C. Hansen, Parameter Estimation for Single Diode Models of 
%   Photovoltaic Modules, Sandia National Laboratories Report SAND2015-2065



if length(I)~=length(V)
    msgbox({'Error in pvl_est_single_diode_param'; 'current and voltage vectors are different length'})
    return
end

I=I(:);
V=V(:);
Isc = I(1);
Voc = V(end);

warning('off')

% fit quadratic spline to IV curve in order to compute the co-content 
% (i.e., integral of Isc - I over V) more accurately
[A, xk, xI, kflag]=Schumaker_QSpline(V,I);

% calculate Co-content integral by numerical integration of quadratic spline for (Isc-I) over V 
xN=length(xk);
xk2=xk(2:xN);
xk1=xk(1:(xN-1));
delX=xk2-xk1;
tu=kflag(1:(xN-1))==true;
tmp=[1/3 1/2 1];
ss=repmat(tmp,xN-1,1);
cc=A.*ss;
tmpInt=[0; sum(cc.*[delX.^3 delX.^2 delX],2)];

sCC=zeros(xN,1);

% Use trapezoid rule for the first 5 intervals due to spline being unreliable near
% the left endpoint
%sCC(1:5)=cumtrapz(xk(1:5),Isc-xI(1:5));  % by trapezoid
sCC(1:5)=Isc*xk(1:5) - cumsum(tmpInt(1:5));  % by spline
sCC(6:(xN-5))=Isc*(xk(6:(xN-5))-xk(5)) - cumsum(tmpInt(6:(xN-5))) + sCC(5);

% Use trapezoid rule for the last 5 intervals due to spline being unreliable near
% the right endpoint
sCC((xN-4):xN)=Isc*(xk((xN-4):xN)-xk(xN-5)) - cumsum(tmpInt((xN-4):xN)) ...
    + sCC(xN-5); % by spline

% For estimating diode equation parameters only use original data points,
% not at any knots added by the quadratic spline fit
CC=sCC(~kflag); % co-content integral, i.e., Int_0^Voc (Isc - I) dV


X=[V(:) Isc-I(:) V(:).*(Isc-I(:)) V(:).*V(:) (I(:)-Isc).^2];  % predictor variables for regression of CC

% define principal components transformation to shift, scale and rotate
% V and I before the regression
% These two lines use functions from the Statistics toolbox, and are
% replaced by the basic Matlab code in the following lines (up to setting
% the value of ev1)

% [tR]=princomp(zscore(X(:,1:2)));  % returns eigenvectors (loadings) in columns ordered by decreasing eigenvalue magnitude
% ev1=tR(:,1);  % first component

tmpX = X(:,1:2);
tmpX_length = size(tmpX,1);

tmpX_mean = mean(tmpX);
tmpX_std = std(tmpX);
tmpX_zscore = (tmpX - repmat(tmpX_mean,tmpX_length,1))./repmat(tmpX_std,tmpX_length,1);

[tmpX_V tmpX_D] = eig(cov(tmpX_zscore));
[~, idx] = sort(diag(tmpX_D),'descend');

ev1 = tmpX_V(:,idx(1));

% second component set to be orthogonal and rotated counterclockwise by 90.
% The Matlab function princomp sometimes returns the second eigenvector as 
% a clockwise rotation because of its embedded sign convention
ev2=[0 -1;1 0]*ev1; 
R=[ev1 ev2];  % principal components transformation

%S=zscore(X(:,1:2))*R'; % [V I] shift and scaled by zscore, rotated by R
S = tmpX_zscore*R';

sCC=(CC(:)-mean(CC)); % center co-content values
col1=ones(length(sCC),1);

% predictors. Shifting makes a constant term necessary in the regression model
sX=[S(:,1) S(:,2) S(:,1).*S(:,2) S(:,1).*S(:,1) S(:,2).*S(:,2) col1];


gamma=sX\sCC; % coefficients from regression in rotated coordinates

% matrix which relates principal components transformation R to the mapping
% between [V' I' V'I' V'^2 I'^2] and sX, where prime ' indicates shifted
% and scaled data.  Used to translate from regression coefficients in
% rotated coordinates to coefficients in initial V,I coordinates.
mB = [R(1,1) R(1,2)   0                             0                0;...
      R(2,1) R(2,2)   0                             0                0;...
      0      0        R(1,1)*R(2,2)+R(1,2)*R(2,1)   2*R(1,1)*R(2,1)  2*R(1,2)*R(2,2);...
      0      0        R(1,1)*R(1,2)                 R(1,1)^2         R(1,2)^2;...
      0      0        R(2,1)*R(2,2)                 R(2,1)^2         R(2,2)^2];

% matrix which is used to undo effect of shifting and scaling on regression
% coefficients.
mA = [std(V)    0          std(V)*mean(Isc-I) 2*std(V)*mean(V)  0;...
      0         std(Isc-I) std(Isc-I)*mean(V) 0                 2*std(Isc-I)*mean(Isc-I);...
      0         0          std(V)*std(Isc-I)  0                 0;...
      0         0          0                  std(V)^2          0;...
      0         0          0                  0                 std(Isc-I)^2];

% translate from coefficients in rotated space (gamma) to
% coefficients in original coordinates (beta)
beta=(mB*mA)\gamma(1:5);

%beta = X\(sX*gamma+repmat(mean(CC),size(sX*gamma)));

% Extract five parameter values from coefficients in original coordinates
% Equation 11, [1]
betaGp = beta(4)*2;

% Equation 12, [1]
betaRs = (sqrt(1+16*beta(4)*beta(5)) - 1)/(4*beta(4));

% Equation 13, [1]
betan = (beta(1)*(sqrt(1+16*beta(4)*beta(5)) - 1) + 4*beta(2)*beta(4))/(4*beta(4)*NsVth);

% Single diode equation at Voc, approximating Iph + Io by Isc

betaIo = (Isc - Voc*betaGp)/(exp(Voc/(betan*NsVth)));

% Single diode equation at Isc, using Rsh, Rs, n and Io determined here

betaIph = Isc - betaIo + betaIo*exp(Isc/(betan*NsVth)) ...
            + Isc*betaRs*betaGp;

Iph = betaIph;
Rs  = betaRs;
Rsh = 1/betaGp;
n   = betan;
Io  = betaIo;

