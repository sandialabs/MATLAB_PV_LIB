function beta = pvl_robustfit(X, Y, intercept)

% PVL_ROBUSTFIT Regress Y onto X using iteratively reweighted least squares
%
% Syntax
%   beta = pvl_robustfit(X, Y, intercept)
%
% Description
%   PVL_ROBUSTFIT uses the Matlab function robustfit if the Statistics 
%   toolbox is available. Otherwise it uses the RLM function in the Python
%   statsmodels package, which requires a Python installation.
%
% Inputs
%   X - must be a N x M matrix, where M is the number of predictors and N is
%      the number of observations
%   Y - must be N x 1 vector
%   intercept - an integer, 0 indicates no intercept, any other value
%     indicates that the regression has an intercept term
%
% Output
%   beta - a vector of coefficients, M x 1 if no intercept is specified, 
%     (M+1) x 1 if an intercept is specified. The coefficient for the 
%     intercept (if present) is beta(1).

% See if robustfit is available (requires statistics toolbox)
tmp = which('robustfit');

if isempty(tmp)
    % need python code
    msg = 'Statistics toolbox not found';
    disp(msg);
    
    % check for python installation
    tmp2 = dos('python --version');
    
    if tmp2==0
        % found python
        msg = 'Found Python installation; using statsmodels package';
        disp(msg);
    
        loc = which('py_rlm.py');  % find python script

        tmpdir = getenv('TEMP');

        f1 = fopen([tmpdir '\temp_rlm_config.txt'],'w');
        fprintf(f1,'%i\n',intercept);
        fclose(f1);

        csvwrite([tmpdir '\temp_rlm_input.csv'],[X Y])

        dos(['python ' loc]);

        beta = csvread([tmpdir '\temp_rlm_output.csv']);

        % clean up temporary files
        delete([tmpdir '\temp_rlm_config.txt'], [tmpdir '\temp_rlm_input.csv'], [tmpdir '\temp_rlm_output.csv']);

    else
        msg = 'Python installation not found';
        disp(msg);
        beta = NaN; 
    end
        
else
    if intercept
        beta = robustfit(X,Y);
    else
        beta = robustfit(X,Y,'bisquare',4.685,'off'); % const is 'off' to omit the intercept.
    end
end
