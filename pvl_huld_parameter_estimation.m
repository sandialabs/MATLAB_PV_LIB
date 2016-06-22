function [Model] = pvl_huld_parameter_estimation(Pmp, Ee, Tm)
% PVL_HULD_PARAMETER_ESTIMATION estimates parameters for the Huld module performance model
%
% Syntax
%   [Model] = pvl_huld_parameter_estimation(Pmp, Ee, Tm)
%
% Description
%   pvl_huld_parameter_estimation estimates parameters for the Huld module
%   performance model [1]. The estimation uses robust regression to fit the
%   Huld model, a polynomial in Tm and log(Ee), to Pmp.
%
% Input:
%   Pmp - a N x 1 vector of power (W) at the maximum power point.
%   Ee - a N x 1 vector of effective irradiance (suns).
%   Tm - a N x 1 vector of module (not cell) temperature (C).
%
% Output:
%   Model - a structure containing the model parameters
%     Model.Pmp0 - estimated Pmp at STC.
%     Model.k - a vector of length 6 containing the coefficients k1 through k6.
%
% Sources:
%   [1] A power-rating model for crystalline silicon PV modules, T. Huld,
%   G. Friesen, A. Skoczek, R. Kenny, T. Sample, M. Field, E. Dunlop, Solar
%   Energy Materials and Solar Cells 95(2011), pp 3359-3369.


Y = Pmp./(Ee);

x1 = log(Ee);
x2 = Tm - 25;

beta = pvl_robustfit([x1 x1.^2 x2 x2.*x1 x2.*x1.^2 x2.^2],Y,true);

Model.Pmp0 = beta(1);
Model.k = beta(2:7);


