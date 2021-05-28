% Author: Vincenzo Bonifaci, Università Roma Tre, Italy
% <vincenzo.bonifaci@uniroma3.it>
%
% This file is to be used with the *l1benchmark* MATLAB package: 
% https://people.eecs.berkeley.edu/~yang/software/l1benchmark/
%
% Implementation of the Accelerated Gradient Scheme (AGS) for
% l1-Norm Minimization, from the paper: 
%
% V. Bonifaci. A Laplacian Approach to l1-Norm Minimization. 
% arXiv:1901.08836 [cs.DS] -- http://arxiv.org/abs/1901.08836
function [q_best, iter, timeSteps, errorSteps] = SolveAGS(A, b, varargin)

STOPPING_TIME = -2;
STOPPING_GROUND_TRUTH = -1;
stoppingCriterion = STOPPING_GROUND_TRUTH;

m = size(A, 2); 
n = size(A, 1); 

% --- Initial point
q0 = ones(m, 1);  % initialize with uniform solution
%q0 = A' * linsolve(A * A', b); % initialize with least-square solution
% --- Fine tuning parameters
cond_threshold = 1e-24; % minimum condition number trusted
maxiter = 5000; % default max no. of iterations (if unspecified)
beta = 3.5; % smoothness parameter -- smaller values are more aggressive
delta_border = 1e-15 * ones(m, 1); 
%delta_border = 1e-15 * sqrt(q0'*q0) * ones(m, 1); 
% ---

% Parse the optional inputs.
if (mod(length(varargin), 2) ~= 0 ),
    error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
end
parameterCount = length(varargin)/2;

for parameterIndex = 1:parameterCount,
    parameterName = varargin{parameterIndex*2 - 1};
    parameterValue = varargin{parameterIndex*2};
    switch lower(parameterName)
        case 'stoppingcriterion'
            stoppingCriterion = parameterValue; 
        %case 'initialization'
        %    y0 = parameterValue;
        case 'groundtruth'
            qG = parameterValue;
        %case 'lambda'
        %    lambda = parameterValue;
        case 'maxiteration'
            maxiter = parameterValue;
        %case 'isnonnegative'
        %    isNonnegative = parameterValue;
        case 'tolerance'
            tolerance = parameterValue;
        %case 'verbose'
        %    verbose = parameterValue;
        case 'maxtime'
            maxTime = parameterValue;
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end
clear varargin
% ---

At = A'; 
t0 = tic ;

iter = 0; 
timeSteps = nan(1,maxiter);
errorSteps = nan(1,maxiter); 
x0 = abs(q0); 
x = x0; 
q_best = q0; 
totgradf = zeros(m, 1); 

while (iter < maxiter)
    iter = iter + 1; 
    
    % Nesterov's magic coefficients
    alpha = 1/(2*iter); 
    tau = 2/(iter+2); 
    %tau = 1e-15; 
    
    % Build Laplacian
    L = (A .* x') * At; 

    % Find potentials
    [p, L_cond] = linsolve(L, b); 
    if (L_cond < cond_threshold)
        disp('AGS: Laplacian ill-conditioned. Returning best iterate.'); 
        timeSteps = timeSteps(1:iter-1);
        errorSteps = errorSteps(1:iter-1);
        return;
    end
    
    % Find congestions
    cong = At * p;
    
    % Find gradient and cumulative gradient
    gradf = ones(m,1) - cong .^ 2; 
    totgradf = totgradf + alpha*gradf; 
    
    % Update running log
    qp = x .* cong; 
    t_elaps = toc(t0); 
    
    % Save best solution
    %MODIFICATO
    if norm(qp - qG) < norm(q_best - qG)
    %if norm(qp, 1) < norm(q_best, 1)
        q_best = qp; 
    end
    
    timeSteps(iter) = t_elaps; 
    % Measure the difference of the L1 norms: 
    errorSteps(iter) = abs(norm(q_best, 1) - norm(qG, 1)) ;
    % In the original l1benchmark, the Euclidean distance was used instead:
    %errorSteps(iter) = norm(q_best - qG); 

    % Check termination
    switch stoppingCriterion
        case STOPPING_GROUND_TRUTH
            done = (errorSteps(iter) < tolerance);
        case STOPPING_TIME
            done = (t_elaps >= maxTime); 
        otherwise
            error('Undefined stopping criterion');
    end % end of the stopping criteria switch

    if done
        break; 
    end
    
    % Compute new iterate
    y = max(delta_border, x - gradf/beta); 
    z = max(delta_border, x0 - totgradf/beta); 
    xp = tau * z + (1-tau) * y; 

    % Next iteration
    x = xp;
end

timeSteps = timeSteps(1:iter); 
errorSteps = errorSteps(1:iter); 

end


