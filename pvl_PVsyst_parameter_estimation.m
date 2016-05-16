function [PVsyst, oflag] = pvl_PVsyst_parameter_estimation(IVCurves, Specs, Const, maxiter, eps1, graphic)
% PVL_PVSYST_PARAMETER_ESTIMATION estimates parameters for the PVsyst module
% performance model
%
% Syntax
%   [PVsyst oflag] = pvl_PVsyst_parameter_estimation(IVCurves, Specs, Const, maxiter, eps1, graphic)
%
% Description
%   pvl_PVsyst_parameter_estimation estimates parameters for the PVsyst module
%   performance model [1,2,3]. Estimation methods are
%   documented in [4, 5, 6].
%
% Input:
%   IVCurves - a structure containing IV curve data in the following fields
%     IVCurves(i).I = vector of current (A) (same length as V)
%     IVCurves(i).V = vector of voltage (V) (same length as I)
%     IVCurves(i).Ee = effective irradiance (W/m^2), i.e., POA broadband
%       irradiance adjusted by solar spectrum modifier
%     IVCurves(i).Tc = cell temperature (C)
%     IVCurves(i).Isc = short-circut current of IV curve (A)
%     IVCurves(i).Voc = open-curcut voltage of IV curve (V)
%     IVCurves(i).Imp = current at max power point of IV curve (A)
%     IVCurves(i).Vmp = voltage at max power point of IV curve (V)
%
%   Specs - a structure containing module-level values
%     Specs.Ns - number of cells in series
%     Specs.aIsc - temperature coefficient of Isc (A/C)
%
%   Const - a structure containing physical and other constants
%     Const.E0 - effective irradiance at STC, normally 1000 W/m2
%     Const.T0 - cell temperature at STC, normally 25 C
%     Const.k = 1.38066E-23 J/K (Boltzmann's constant)
%     Const.q = 1.60218E-19 Coulomb (elementary charge)
%
% Optional inputs
%   maxiter - an integer setting the maximum number of iterations for the 
%     parameter updating part of the algorithm. Default value is 5
%
%   eps1 - the desired tolerance for convergence for the IV curve fitting.
%     The iterative parameter updating stops when absolute values of the
%     percent change in mean, max and standard deviation of Imp, Vmp and Pmp
%     between iterations are all less than eps1, or when the number of 
%     iterations exceeds maxiter.  Default value of eps1 is 1e-3 (0.0001%). 
%
%   graphic - a boolean, if true then plots are produced during the 
%     parameter estimation process. Default is false
%
%
% Output:
%   PVsyst - a structure containing the model parameters
%     PVsyst.IL_ref - light current (A) at STC
%     PVsyst.Io_ref - dark current (A) at STC
%     PVsyst.eG - effective band gap (eV) at STC
%     PVsyst.Rsh_ref - shunt resistance (ohms) at STC
%     PVsyst.Rsh0 - shunt resistance (ohms) at zero irradiance 
%     PVsyst.Rshexp - exponential factor defining decrease in 
%          Rsh with increasing effective irradiance
%     PVsyst.Rs_ref - series resistance (ohms) at STC
%     PVsyst.gamma_ref - diode (ideality) factor at STC
%     PVsyst.mugamma - temperature coefficient for diode (ideality) factor
%     PVsyst.Iph - vector of values of light current Iph estimated for each IV
%       curve
%     PVsyst.Io - vector of values of dark current Io estimated for each IV
%       curve
%     PVsyst.Rsh - vector of values of shunt resistance Rsh estimated for each IV
%       curve
%     PVsyst.Rs - vector of values of series resistance Rs estimated for each IV
%       curve
%     PVsyst.u - filter indicating IV curves with parameter values deemed 
%       reasonable by the private function filter_params
%     oflag - Boolean indicating success or failure of estimation of the
%    diode (ideality) factor parameter.  If failure, then no parameter values
%    are returned.
%
% Sources:
%
% [1] K. Sauer, T. Roessler, C. W. Hansen, Modeling the Irradiance and 
%     Temperature Dependence of Photovoltaic Modules in PVsyst, 
%     IEEE Journal of Photovoltaics v5(1), January 2015.
%
% [2] A. Mermoud, PV modules modelling, Presentation at the 2nd PV
%     Performance Modeling Workshop, Santa Clara, CA, May 2013
%
% [3] A. Mermoud, T. Lejeune, Performance Assessment of a Simulation Model
%     for PV modules of any available technology, 25th European Photovoltaic
%     Solar Energy Conference, Valencia, Spain, Sept. 2010
%
% [4] C. Hansen, Estimating Parameters for the PVsyst Version 6 Photovoltaic
%     Module Performance Model, Sandia National Laboratories Report
%     SAND2015-8598
%
% [5] C. Hansen, Parameter Estimation for Single Diode Models of 
%     Photovoltaic Modules, Sandia National Laboratories Report
%     SAND2015-2065
%
% [6] C. Hansen, Estimation of Parameters for Single Diode Models using
%     Measured IV Curves, Proc. of the 39th IEEE PVSC, June 2013.


% Set max iterations to timeout if convergence parameters are not met
if isnan(maxiter)
    maxiter = 5; % default value
end

if isnan(eps1)
    eps1 = 1e-3; % default value
end

Ee = [IVCurves.Ee]';
Tc = [IVCurves.Tc]';
TcK = Tc + 273.15;
Isc = [IVCurves.Isc]';
Voc = [IVCurves.Voc]';
Imp = [IVCurves.Imp]';
Vmp = [IVCurves.Vmp]';

% Cell thermal voltage
Vth = Const.k/Const.q*(Tc+273.15);

N = length(IVCurves);

% Initial estimate of Rsh used to obtain the diode factor gamma0 and diode
% temperature coefficient mugamma.  Rsh is estimated using the co-content
% integral method.

pIo = NaN(N,1);
pIph = NaN(N,1);
pRsh = NaN(N,1);
pRs = NaN(N,1);
pn = NaN(N,1);

for i=1:N
    [I, V] = pvl_rectify_IV_curve(IVCurves(i).I, IVCurves(i).V, Voc(i), Isc(i));
    % Initial estimate of Rsh, from integral over voltage and regression
    % [4] Step 3a; [5] Step 3a 
    [pIo(i), pIph(i), pRs(i), pRsh(i), pn(i)] = est_single_diode_param(I, V, Vth(i)*Specs.Ns);
end;

pRsh = pRsh(:);

% Estimate the diode factor gamma from Isc-Voc data. Method incorporates temperature dependence
% by means of the equation for Io
Y = log(Isc - Voc./pRsh) - 3*log(TcK./(Const.T0+273.15));
X1 = Const.q/Const.k*(1./(Const.T0+273.15) - 1./TcK);
X2 = Voc./(Vth*Specs.Ns);
uu = isnan(Y)|isnan(X1)|isnan(X2);

X = [ones(size(X1(~uu))) X1(~uu) -X1(~uu).*(TcK(~uu) - (Const.T0+273.15)) X2(~uu) -X2(~uu).*(TcK(~uu)-(Const.T0+273.15))];
alpha = X\Y(~uu);

gamma_ref = 1/alpha(4);
mugamma = alpha(5)/alpha(4)^2;

if isnan(gamma_ref) || isnan(mugamma) || imag(gamma_ref)~=0 || imag(mugamma)~=0
    badgamma = true;
else
    badgamma = false;
end

if ~badgamma
    
    gamma = gamma_ref + mugamma*(Tc - Const.T0);

    PVsyst.gamma = gamma;

    if graphic
        figure
        plot(X2,Y,'b+')
        hold on
        plot(X2, X*alpha, 'r.')
        xlabel('X = Voc / Ns \times Vth')
        ylabel('Y = log(Isc - Voc/Rsh)')
        legend('I-V Data', 'Regression model', 'location', 'NorthWest')
        box on
        text(min(X2)+0.85*(max(X2)-min(X2)),1.05*max(Y),['\gamma_0 = ' num2str(gamma_ref)]);
        text(min(X2)+0.85*(max(X2)-min(X2)),0.98*max(Y),['\mu_\gamma = ' num2str(mugamma)]);
    end

    nNsVth = gamma.*(Vth*Specs.Ns);

    % display progress bar, which shows fraction of iterations complete
    hw=waitbar(0,'Initial values');


    %% For each IV curve, sequentially determine initial values for Io, Rs, and Iph
    % [4] Step 3a; [5] Step 3

    Io = NaN(N,1);
    Iph = NaN(N,1);
    Rs = NaN(N,1);
    Rsh = pRsh;

    for i=1:N
        [I, V] = pvl_rectify_IV_curve(IVCurves(i).I, IVCurves(i).V, Voc(i), Isc(i));

        if Rsh(i)>0
            % Initial estimate of Io, evaluate the single diode model at Voc
            % and approximate Iph + Io = Isc
            % [4] Step 3a; [5] Step 3b 
            Io(i) = (Isc(i) - Voc(i)/Rsh(i))*exp(-Voc(i)/nNsVth(i));

            % Initial estimate of Rs from dI/dV near Voc
            % [4] Step 3a; [5] Step 3c
            dIdV = numdiff(V,I);
            u = V>0.5*Voc(i) & V<0.9*Voc(i);
            tmp = -Rsh(i)*dIdV-1;
            v = u & (tmp>0);
            if sum(v)>0
                vtRs = nNsVth(i)/Isc(i)* ...
                    (log(tmp(v)*nNsVth(i)/(Rsh(i)*Io(i))) - ...
                    V(v)/nNsVth(i));
                Rs(i) = mean(vtRs(vtRs>0));
            else
                Rs(i) = 0;
            end

            % Initial estimate of Iph, evaluate the single diode model at Isc
            % [4] Step 3a; [5] Step 3d
            Iph(i) = Isc(i) - Io(i) + Io(i)*exp(Isc(i)/nNsVth(i)) ...
                + Isc(i)*Rs(i)/Rsh(i);
        else % Rsh came back negative
            Io(i) = NaN;
            Rs(i) = NaN;
            Iph(i) = NaN;
        end
    end


    % Filter IV curves for good initial values
    % [4] Step 3b
    u = filter_params(Io, Rsh, Rs, Ee, Isc);

    % Refine Io to match Voc
    % [4] Step 3c
    tmpIph = Iph;
    tmpIo = update_Io_known_n(Rsh(u), Rs(u), nNsVth(u), Io(u), tmpIph(u), Voc(u));
    Io(u) = tmpIo;

    % Calculate Iph to be consistent with Isc and current values of other parameters
    % [5], Step 3c
    Iph = Isc - Io + Io.*exp(Rs.*Isc./nNsVth) + Isc.*Rs./Rsh;

    %% Refine Rsh, Rs, Io and Iph in that order.
    i = 1;  % counter variable for parameter updating while loop, counts iterations
    PrevConvergeParams = struct('State',0,'VmpErrMeanChange',Inf); % Initialize a struct for PrevConvergeParams, required for first run through of check_converge
    if graphic
        h = figure();
    end
    if graphic
        ConvergeParamsFig = figure();   % Create a new handle for the Convergence Parameter Figure
    end

    waitbar(0,hw,'Updating parameters');
    drawnow; pause(0.05);

    while (((PrevConvergeParams.VmpErrMeanChange >= eps1) || ...
            (PrevConvergeParams.ImpErrMeanChange >= eps1) || ...
            (PrevConvergeParams.PmpErrMeanChange >= eps1) || ...
            (PrevConvergeParams.VmpErrStdChange >= eps1) || ...
            (PrevConvergeParams.ImpErrStdChange >= eps1) || ...
            (PrevConvergeParams.PmpErrStdChange >= eps1) || ...
            (PrevConvergeParams.VmpErrAbsMaxChange >= eps1) || ...
            (PrevConvergeParams.ImpErrAbsMaxChange >= eps1) || ...
            (PrevConvergeParams.PmpErrAbsMaxChange >= eps1)) && (i <= maxiter)) 
        % update waitbar to show number of iterations complete
        waitbar(i/maxiter,hw);

        % Update Rsh to match max power point using a fixed point method.
        [tmpRsh] = update_Rsh_fixed_pt(Rsh(u), Rs(u), Io(u), Iph(u), ...
            nNsVth(u), Imp(u), Vmp(u)); 
        if graphic
            figure(h)
            scatter(i, mean(abs(tmpRsh - Rsh(u))), 5, 'k', 'filled');
            hold on;
            title('update Rsh')
            ylabel('mean(abs(tmpRsh(u) - Rsh(u)))')
            xlabel('Iteration')
        end
        Rsh(u) = tmpRsh;

        % Calculate Rs to be consistent with Rsh and maximum point point
        [~, phi] = calc_theta_phi_exact(Imp(u), Iph(u), Vmp(u), Io(u), ...
            nNsVth(u), Rs(u), Rsh(u));
        Rs(u) = (Iph(u)+Io(u)-Imp(u)).*Rsh(u)./Imp(u) - ...
            nNsVth(u).*phi./Imp(u) - Vmp(u)./Imp(u);

        % Update filter for good parameters
        u = filter_params(Io, Rsh, Rs, Ee, Isc);

        % Update value for Io to match Voc
        [tmpIo] = update_Io_known_n(Rsh(u), Rs(u), nNsVth(u), Io(u), Iph(u), Voc(u));
        Io(u) = tmpIo;

        % Calculate Iph to be consistent with Isc and other parameters
        Iph = Isc - Io + Io.*exp(Rs.*Isc./nNsVth) + Isc.*Rs./Rsh;

        % Update filter for good parameters
        u = filter_params(Io, Rsh, Rs, Ee, Isc);

        % compute the IV curve from the current parameter values
        Results = pvl_singlediode(Iph(u), Io(u), Rs(u), Rsh(u), nNsVth(u));

        % Check convergence criteria
        % [4] Step 3d
        if graphic
            ConvergeParams = check_converge(PrevConvergeParams, Results, Vmp(u), Imp(u), graphic, ConvergeParamsFig, i);
        else
            ConvergeParams = check_converge(PrevConvergeParams, Results, Vmp(u), Imp(u), graphic, 0, i);
        end        
        PrevConvergeParams = ConvergeParams;
        i = i+1;
    end

    if i==maxiter
        waitbar(1,hw)
    end

    %% Extract coefficients for auxillary equations

    % Estimate Io0 and eG
    TcK = Tc + 273.15; % Convert Tc to K
    T0K = Const.T0 + 273.15; % convert T0 to K
    X = Const.q/Const.k*(1/T0K - 1./TcK(u))./gamma(u);
    Y = log(Io(u))-3*log(TcK(u)/T0K);
    beta = robustfit(X,Y);
    Io0 = exp(beta(1));
    eG = beta(2);


    if graphic
        % Predict Io and Eg
        pIo = Io0*((Tc(u)+273.15)/(Const.T0+273.15)).^3.*...
            exp((Const.q/Const.k)*(eG./gamma(u)).*(1./(Const.T0+273.15)-1./(Tc(u)+273.15)));

        figure
        subplot(311)
        plot(Tc(u),Y,'r+')
        hold all
        plot(Tc(u),beta(1) + X*beta(2),'b.')
        xlabel('Cell temp. (C)')
        ylabel('log(Io)-3log(T_C/T_0)')
        legend('Data','Model','Location','NorthWest')

        subplot(312)
        plot(Tc(u),Io(u),'r+')
        hold all
        plot(Tc(u),pIo,'.')
        xlabel('Cell temp. (C)')
        ylabel('I_O (A)')
        legend('Extracted','Predicted','Location','NorthWest')

        subplot(313)
        plot(Tc(u),(pIo-Io(u))./Io(u)*100,'x')
        xlabel('Cell temp. (C)')
        ylabel('Percent Deviation in I_O')
        refline(0,0)

        figure('Position',[1 1 600 300])
        plot(Tc(u),Y  + 3*(Tc(u)/Const.T0),'k.')
        hold all
        plot(Tc(u),beta(1) + X*beta(2) + 3*(Tc(u)/Const.T0),'g.')
        xlabel('Cell temp. (C)')
        ylabel('log(Io)-3log(T_C/T_0)')
        xlabel('Cell temp. (C)', 'FontSize',15,'FontWeight','bold')
        ylabel('ln(Io)-3ln(T_C/T_0)', 'FontSize',15,'FontWeight','bold')
        legend('Data','Regression Model','Location','NorthWest')

        figure
        plot(Tc(u),Io(u),'b+')
        hold all
        plot(Tc(u),pIo,'r.')
        xlabel('Cell temp. \circC')
        ylabel('I_o (A)')
        legend('Extracted from IV Curves','Pred. by Eq. 3','Location','NorthWest')
        text(min(Tc(u))+(max(Tc(u))-min(Tc(u)))*0.0,min(Io(u))+(max(Io(u))-min(Io(u)))*0.83,['I_{O0} = ' num2str(Io0)]);
        text(min(Tc(u))+(max(Tc(u))-min(Tc(u)))*0.0,min(Io(u))+(max(Io(u))-min(Io(u)))*0.75,['\epsilon_G = ' num2str(eG)]);
    end

    % Estimate Iph0
    X = (Tc(u)-Const.T0);
    Y = Iph(u).*(Const.E0./Ee(u));
    Iph0 = nanmean(Y-Specs.aIsc*X);

    if graphic
        % predict Iph
        pIph = (Ee(u)/Const.E0).*(Iph0+Specs.aIsc*(Tc(u)-Const.T0));
        figure
        subplot(311)
        plot(Ee(u), Y,'r+')
        hold all
        line([0 max(Ee(u))],[Iph0 Iph0])
        xlabel('Irradiance (W/m^2)')
        ylabel('I_L')
        legend('Data','I_L at STC','Location','SouthEast')

        subplot(312)
        plot(Ee(u),Iph(u),'r+')
        hold all
        %(E(u)/Const.E0).*(Iphi0+mIsc*(Tc(u)-Const.T0))
        plot(Ee(u),pIph,'.');
        xlabel('Irradiance (W/m^2)')
        ylabel('I_L (A)')
        legend('Extracted','Predicted','Location','NorthWest')

        subplot(313)
        plot(Ee(u),(pIph-Iph(u))./Iph(u)*100,'x')
        xlabel('Irradiance (W/m^2)')
        ylabel('Percent Deviation from I_{ L}')
        refline(0,0)
        
        figure
        plot(Tc(u),Iph(u),'b+');
        hold all
        plot(Tc(u),pIph,'r.');
        line([0 80],[Iph0 Iph0]);
        text(1.1*min(Tc(u)),1.05*Iph0,['I_{L0} = ' num2str(Iph0)]);
        xlabel('Cell temp. \circC');
        ylabel('I_L (W/m^2')
        legend('Extracted from IV Curves','Pred. by Eq. 2','I_L at STC','Location','NorthWest')
        
    end

    % Additional filter for Rsh and Rs; restrict effective irradiance to be
    % greater than 400 W/m2

    v = Ee>400;

    % Estimate Rsh0, Rsh_ref and Rshexp

    % Initial guesses.  Rsh0 is value at Ee=0.
    if any(Ee<400)
        gRsh0 = nanmean(Rsh(Ee<400));
    else
        gRsh0 = max(Rsh);
    end
    % Rsh_ref is value at Ee=1000
    if any(Ee>400)
        gRshref = nanmean(Rsh(Ee>400));
    else
        gRshref = min(Rsh);
    end
    % PVsyst default for Rshexp is 5.5
    Rshexp = 5.5;

    % Here we use a nonlinear least squares technique.  lsqnonlin minimizes the
    % sum of squares of the objective function (here, tf).
    x0 = [gRsh0 gRshref];
    tf = @(x) log10(estRsh(x,Rshexp,Ee(u),Const.E0)) - log10(Rsh(u));
    options = optimset('Display','off');  % Notify only on non-convergence
    beta = lsqnonlin(tf,x0,[1 1],[1e7 1e6],options);

    % Extract PVsyst parameter values
    Rsh0 = beta(1);
    Rshref = beta(2);

    if graphic
        % Predict Rsh
        pRsh = estRsh(beta,Rshexp,Ee,Const.E0);
        figure
        subplot(211)
        plot(Ee(u),log10(Rsh(u)),'r.')
        hold all
        plot(Ee(u),log10(pRsh(u)),'b.')
        %ylim([2 6])
        xlabel('Irradiance (W/m^2)')
        ylabel('log_{10}(R_{sh})')
        legend('Extracted','Predicted','Location','NorthWest')
%        ylim([2 4.5])

        subplot(212)
        plot(Ee(u),(log10(pRsh(u)) - log10(Rsh(u)))./log10(Rsh(u))*100,'x')
        xlabel('Irradiance (W/m^2)')
        ylabel('Percent Deviation in log_{10}(R_{sh})')
        refline(0,0)
%        ylim([-35 15])

        figure
        plot(Ee(u),log10(Rsh(u)),'b.')
        hold all
        plot(Ee(u),log10(pRsh(u)),'r.')
        %ylim([2 6])
        xlabel('Irradiance (W/m^2)')
        ylabel('log_{10}(R_{sh})')
        legend('Extracted from IV Curves','Pred. by Eq. 5','Location','SouthWest')
        text(150,3.65,['R_{SH0} = ' num2str(Rsh0)]);
        text(150,3.5,['R_{SH,ref} = ' num2str(Rshref)]);
        text(150,3.35,['R_{SHexp} = ' num2str(Rshexp)]);
        
    end

    % Estimate Rs0
    Rs0 = mean(Rs(u&v));

    if graphic
        figure
        subplot(211)
        plot(Ee(u&v),Rs(u&v),'r.')
        hold all
        plot(Ee(u&v),Rs0*ones(size(Ee(u&v))),'b.')
        xlabel('Irradiance (W/m^2)')
        ylabel('R_S')
        ylim([0 1]);
        xlim([0 1200])
        legend('R_S values','Model')

        subplot(212)
        plot(Ee(u),(Rs0-Rs(u))./Rs(u)*100,'x')
        xlabel('Irradiance (W/m^2)')
        ylabel('Percent Deviation in R_S')
        refline(0,0)
        
        figure
        plot(Ee(u&v),Rs(u&v),'b.')
        hold all
        line([0 max(Ee(u))],[Rs0 Rs0],'Color','r');
        xlabel('Irradiance (W/m^2)')
        ylabel('R_S')
%         ylim([0 1]);
%         xlim([0 1200])
        legend('Extracted from IV Curves','Pred. by Eq. 7','Location','SouthWest')
        text(800,1.2*Rs0,['R_{S0} = ' num2str(Rs0)]);
    end


    %% Save parameter estimates in output structure
    PVsyst.IL_ref = Iph0;
    PVsyst.I0_ref = Io0;
    PVsyst.eG = eG;
    PVsyst.Rs_ref = Rs0;
    PVsyst.gamma_ref = gamma_ref;
    PVsyst.mugamma = mugamma;
    PVsyst.Iph = Iph;
    PVsyst.I0 = Io;
    PVsyst.Rsh0 = Rsh0;
    PVsyst.Rsh_ref = Rshref;
    PVsyst.Rshexp = Rshexp;
    PVsyst.Rs = Rs;
    PVsyst.Rsh = Rsh;
    PVsyst.Ns = Specs.Ns;
    PVsyst.u = u;
    
    oflag = true;

    close(hw)
    drawnow; pause(0.05);

else
    
    oflag = false;
    
    % Save parameter estimates in output structure
    PVsyst.IL_ref = NaN;
    PVsyst.I0_ref = NaN;
    PVsyst.eG = NaN;
    PVsyst.Rs_ref = NaN;
    PVsyst.gamma_ref = NaN;
    PVsyst.mugamma = NaN;
    PVsyst.Iph = NaN;
    PVsyst.I0 = NaN;
    PVsyst.Rsh0 = NaN;
    PVsyst.Rsh_ref = NaN;
    PVsyst.Rshexp = NaN;
    PVsyst.Rs = NaN;
    PVsyst.Rsh = NaN;
    PVsyst.Ns = Specs.Ns;
    PVsyst.u = zeros(size(IVCurves));
end

end

function pRsh = estRsh(x,Rshexp,G,G0)

% Computes Rsh for PVsyst model where the parameters are in vector x:
% x(1) = Rsh0;
% x(2) = Rshref;
% x(3) = Rshexp;

Rsh0 = x(1);
Rshref = x(2);

Rshb = max((Rshref - Rsh0*exp(-Rshexp))/(1-exp(-Rshexp)),0);
pRsh = Rshb + (Rsh0 - Rshb)*exp(-Rshexp*(G/G0));

pRsh = pRsh(:);

end

function u = filter_params(Io, Rsh, Rs, Ee, Isc)

% Function filter_params identifies bad parameters sets. A bad set contains
% NaN, non-positive or imaginary values for parameters; Rs > Rsh; or data
% where effective irradiance Ee differs by more than 5% from a linear fit
% to Isc vs. Ee.

badRsh = Rsh<0 | isnan(Rsh);
negRs = ~(Rs>0);
badRs = Rs>Rsh | isnan(Rs);
imagRs = imag(Rs)~=0;
badIo = imag(Io)~=0 | Io<=0;
goodR = ~badRsh & ~imagRs & ~negRs & ~badRs & ~badIo;

eff = Ee/1000\Isc;
pIsc = eff.*Ee/1000;
pIsc_error = abs(pIsc-Isc)./Isc;
badIph = pIsc_error > 0.05;  % check for departure from linear relation between Isc and Ee

u = goodR & ~badIph;

end

function [ConvergeParam] = check_converge(PrevParams, Results, Vmp, Imp, graphic, ConvergeParamsFig,i)
% Function check_converge computes convergence metrics for all IV curves.
%
% Inputs
%       PrevParams: Convergence Parameters from the Previous Iteration (used to determine Percent Change in values between iterations)
%       Results:    Performance parameters of the (predicted) single diode fitting, which includes Voc, Vmp, Imp, Pmp, Isc
%       Vmp, Imp:   Measured values for each IV curve
%       graphic:    Argument to determine whether to display Figures
%       ConvergeParamsFig: Handle to the ConvergeParam Plot
%       i:          Index of current iteration in cec_parameter_estimation
%  
% Outputs
%       ConvergeParam - a structure containing the following for Imp, Vmp
%       and Pmp:
%         - maximum percent difference between measured and modeled values
%         - minimum percent difference between measured and modeled values
%         - maximum absolute percent difference between measured and
%           modeled values
%         - mean percent difference between measured and modeled values
%         - standard deviation of percent difference between measured and
%           modeled values
%
%         - absolute difference for previous and current values of
%            maximum absolute percent difference (measured vs. modeled) 
%         - absolute difference for previous and current values of
%            mean percent difference (measured vs. modeled) 
%         - absolute difference for previous and current values of
%            standard deviation of percent difference (measured vs. modeled) 

ConvergeParam.ImpErrMax = max((Results.Imp-Imp)./Imp*100); % max of the error in Imp
ConvergeParam.ImpErrMin = min((Results.Imp-Imp)./Imp*100); % min of the error in Imp
ConvergeParam.ImpErrAbsMax = max(abs((Results.Imp-Imp)./Imp*100)); % max of the error in Imp
ConvergeParam.ImpErrMean = mean((Results.Imp-Imp)./Imp*100); % mean of the error in Imp
ConvergeParam.ImpErrStd = std((Results.Imp-Imp)./Imp*100); % std of the error in Imp


ConvergeParam.VmpErrMax = max((Results.Vmp-Vmp)./Vmp*100); % max of the error in Vmp
ConvergeParam.VmpErrMin = min((Results.Vmp-Vmp)./Vmp*100); % min of the error in Vmp
ConvergeParam.VmpErrAbsMax = max(abs((Results.Vmp-Vmp)./Vmp*100)); % max of the error in Vmp
ConvergeParam.VmpErrMean = mean((Results.Vmp-Vmp)./Vmp*100); % mean of the error in Vmp
ConvergeParam.VmpErrStd = std((Results.Vmp-Vmp)./Vmp*100); % std of the error in Vmp

ConvergeParam.PmpErrMax = max((Results.Pmp-(Imp.*Vmp))./(Imp.*Vmp)*100); % max of the error in Pmp % CKC Added 2012-07-24 to compute std of Pmp
ConvergeParam.PmpErrMin = min((Results.Pmp-(Imp.*Vmp))./(Imp.*Vmp)*100); % min of the error in Pmp % CKC Added 2012-07-24 to compute std of Pmp
ConvergeParam.PmpErrAbsMax = max(abs((Results.Pmp-(Imp.*Vmp))./(Imp.*Vmp)*100)); % max of the error in Pmp % CKC Added 2012-07-24 to compute std of Pmp
ConvergeParam.PmpErrMean = mean((Results.Pmp-(Imp.*Vmp))./(Imp.*Vmp)*100); % mean of the error in Pmp % CKC Added 2012-07-24 to compute std of Pmp
ConvergeParam.PmpErrStd = std((Results.Pmp-(Imp.*Vmp))./(Imp.*Vmp)*100); % std of the error in Pmp % CKC Added 2012-07-24 to compute std of 


if (PrevParams.State ~= 0)
    ConvergeParam.ImpErrStdChange = abs((ConvergeParam.ImpErrStd - PrevParams.ImpErrStd)/PrevParams.ImpErrStd);
    ConvergeParam.VmpErrStdChange = abs((ConvergeParam.VmpErrStd - PrevParams.VmpErrStd)/PrevParams.VmpErrStd);
    ConvergeParam.PmpErrStdChange = abs((ConvergeParam.PmpErrStd - PrevParams.PmpErrStd)/PrevParams.PmpErrStd);
    ConvergeParam.ImpErrMeanChange = abs((ConvergeParam.ImpErrMean - PrevParams.ImpErrMean)/PrevParams.ImpErrMean);
    ConvergeParam.VmpErrMeanChange = abs((ConvergeParam.VmpErrMean - PrevParams.VmpErrMean)/PrevParams.VmpErrMean);
    ConvergeParam.PmpErrMeanChange = abs((ConvergeParam.PmpErrMean - PrevParams.PmpErrMean)/PrevParams.PmpErrMean);
    ConvergeParam.ImpErrAbsMaxChange = abs((ConvergeParam.ImpErrAbsMax - PrevParams.ImpErrAbsMax)/PrevParams.ImpErrAbsMax);
    ConvergeParam.VmpErrAbsMaxChange = abs((ConvergeParam.VmpErrAbsMax - PrevParams.VmpErrAbsMax)/PrevParams.VmpErrAbsMax);
    ConvergeParam.PmpErrAbsMaxChange = abs((ConvergeParam.PmpErrAbsMax - PrevParams.PmpErrAbsMax)/PrevParams.PmpErrAbsMax);
    ConvergeParam.State = 1;
else
   ConvergeParam.ImpErrStdChange = Inf;
   ConvergeParam.VmpErrStdChange = Inf;
   ConvergeParam.PmpErrStdChange = Inf;
   ConvergeParam.ImpErrMeanChange = Inf;
   ConvergeParam.VmpErrMeanChange = Inf;
   ConvergeParam.PmpErrMeanChange = Inf;
   ConvergeParam.ImpErrAbsMaxChange = Inf;
   ConvergeParam.VmpErrAbsMaxChange = Inf;
   ConvergeParam.PmpErrAbsMaxChange = Inf;
   ConvergeParam.State = 1;
end

if graphic

    figure(ConvergeParamsFig)
    subplot(3,3,1)
    plot(i,ConvergeParam.PmpErrMean,'x-')
    hold on;
    title('Mean of Err in Pmp')
    ylabel('mean((pPmp-Pmp)/Pmp*100)')
    xlabel('Iteration') 

    subplot(3,3,2)
    plot(i,ConvergeParam.VmpErrMean,'x-')
    hold on;
    title('Mean of Err in Vmp')
    ylabel('mean((pVmp-Vmp)/Vmp*100)')
    xlabel('Iteration')

    subplot(3,3,3)
    plot(i,ConvergeParam.ImpErrMean,'x-')
    hold on;
    title('Mean of Err in Imp')
    ylabel('mean((pImp-Imp)/Imp*100)')
    xlabel('Iteration')

    subplot(3,3,4)
    plot(i,ConvergeParam.PmpErrStd,'x-')
    hold on;
    title('Std of Err in Pmp')
    ylabel('std((pPmp-Pmp)/Pmp*100)')
    xlabel('Iteration') 

    subplot(3,3,5)
    plot(i,ConvergeParam.VmpErrStd,'x-')
    hold on;
    title('Std of Err in Vmp')
    ylabel('std((pVmp-Vmp)/Vmp*100)')
    xlabel('Iteration')

    subplot(3,3,6)
    plot(i,ConvergeParam.ImpErrStd,'x-')
    hold on;
    title('Std of Err in Imp')
    ylabel('std((pImp-Imp)/Imp*100)')
    xlabel('Iteration')
        
    subplot(3,3,7)
    plot(i,ConvergeParam.PmpErrAbsMax,'x-')
    hold on;
    title('AbsMax of Err in Pmp')
    ylabel('max(abs((pPmp-Pmp)/Pmp*100))')
    xlabel('Iteration') 

    subplot(3,3,8)
    plot(i,ConvergeParam.VmpErrAbsMax,'x-')
    hold on;
    title('AbsMax of Err in Vmp')
    ylabel('max(abs((pVmp-Vmp)/Vmp*100))')
    xlabel('Iteration')

    subplot(3,3,9)
    plot(i,ConvergeParam.ImpErrAbsMax,'x-')
    hold on;
    title('AbsMax of Err in Imp')
    ylabel('max(abs((pImp-Imp)/Imp*100))')
    xlabel('Iteration')
        
end 

end

