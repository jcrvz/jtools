%% Cuckoo Search via L'evy Flight Algorithm method [1]
%
% CS - performs a metaheuristic procedure to solve a simple-constrained
%        optimisation problem, such as:
%
%                min FOBJ(x),    where X = [X_1, X_2, ..., X_Nd]',
%                 X
%                 subject to     X_lb <= X <= X_ub,
%
% XOPT = CS(FOBJ,BOUNDARIES) searchs a simple-constrained minimal XOPT of
% the fitness function (or objective function) FOBJ. The feasible region is
% drawn by the simple constraints given by BOUNDARIES. BOUNDARIES is a
% matrix of Nd-by-2, where Nd is the number of dimensions or variables.
% Elements of the first column of BOUNDARIES are lower boundaries of each
% design variable, and the elements of the second one are upper boundaries
% of each design variable.
%
% XOPT = CS(FOBJ,BOUNDARIES,PARAMETERS) searchs a simple-constrained
% minimal XOPT of the FOBJ function using a different set of tune parameters
% of the method. PARAMETERS contains tune parameters for CS and it is a
% struct defined as follows:
%
%                                      (default values)
%     PARAMETERS = struct(...
%                           'NAG'       ,   25   , ...  % Number of nests
%                           'BETA'      ,   1.5  , ...  % l'evy distr.
%                           'ALPHA'     ,   0.1  , ...  % step size scale factor
%                           'PA'        ,   0.25  , ...  % probability of discover eggs
%
%                           'EPS1'      ,   0.5  , ...  % historical eps.
%                           'EPS2'      ,   1e-3 , ...  % population eps.
%                           'MSAT'      ,   10   , ...  % max. saturation
%                           'MITE'      ,   1e12 , ...  % max. iteration
%                           'UNCONST'   ,   false  ...  % unconstrained
%                         );
%
% [XOPT,FOPT,DETAILS] = CS(FOBJ,BOUNDARIES,PARAMETERS) performs the above
% process but in this case the final value of fitness function (FOPT) and
% some additional information (DETAILS) are returned. DETAILS is also a
% struct which contains:
%
%       DETAILS.time    - time required (in seconds) to find XOPT
%              .fevs    - evaluations of fitness function done
%              .steps   - number of steps or iterations performed
%              .outmsg  - flag of convergence (1 is convergence)
%
% --------------------------------------------------------------------------
% Reference:
% [1] Yang, X.-S., & Suash Deb. (2009). Cuckoo Search via Levy flights.
%     In 2009 World Congress on Nature & Biologically Inspired Computing
%     (NaBIC) (pp. 210�214). IEEE. http://doi.org/10.1109/NABIC.2009.5393690
% --------------------------------------------------------------------------
%
% {Copyright (c) Jorge M. Cruz-Duarte. All rights reserved.}
%
% Departamento de Ingenier�a El�ctrica,
% Divisi�n de Ingenier�as, Campus Irapuato-Salamanca,
% Universidad de Guanajuato, Salamanca, Guanajuato, M�xico.
%
% Grupo de Investigaci�n CEMOS,
% Escuela de Ingenier�as El�ctrica, Electr�nica y de Telecomunicaciones,
% Universidad Industrial de Santander, Bucaramanga, Santander, Colombia.
%
% Modifications:
%                   2015-jul ver :: v0
%                   2016-apr ver :: v1
%
% Contact us: jorge.cruz@ugto.mx xxxxx ones(nd,1)*details.Constraints,ni,na

function [BestNest,fBest,details] = CS(fObj,bnd)

% Read parameters
% if nargin < 3,
    Na      = 40;                   % Number of nests
    pa      = 0.75;                 % Probability of discover eggs
    beta    = 1.5;                     % Beta parameter to calculate sigma
    alpha   = 1;                 % Step size scale factor

%     eps1    = 1e-3;
%     eps2    = 1e-6;
     M       = 5e2;
    %msat    = 50;

    unconst = false;
% else
%     Na      = parameters.NAG;
%     pa      = parameters.PA;
%     beta    = parameters.BETA;
%     alpha   = parameters.ALPHA;
% 
%     eps1    = parameters.EPS1;
%     eps2    = parameters.EPS2;
%     M       = parameters.MITE;
%     msat    = parameters.MSAT;
% 
%     unconst = parameters.UNCONST;
% end

% Read problem's dimensions
Nd      = size(bnd,1);

% Define a quick function to get Pi values
get_Pi      = @(condition,current_position,best_current_position) ...
    repmat(condition,1,Nd).*current_position + ...
    repmat(~condition,1,Nd).*best_current_position;

% Define a quick function to obtain f(X)
    function the_results = evaluate_function(the_function,the_positions)
        the_results   = nan(Na,1);
        for s = 1 : Na
            the_results(s) = the_function(the_positions(s,:));
        end
    end

% Initial nest positions (random solutions)
bnd     = [min(bnd,[],2) max(bnd,[],2)];
bnd_1       = repmat(bnd(:,1)',Na,1);
bnd_2       = repmat(bnd(:,2)',Na,1);

Nest           = bnd_1 + rand(Na,Nd).*(bnd_2 - bnd_1);

% Initial solutions
Fitness      = evaluate_function(fObj,Nest);

% Find the Nest best position (initial step)
[fBest,g] = min(Fitness); BestNest = Nest(g,:);

% Calculate std dev of u (Mantegna's algorithm)
sigma       = (gamma(1 + beta)*sin(pi*beta/2)/(gamma((1 + beta)/2)*beta*...
    2^((beta - 1)/2)))^(1/beta);

% Set auxiliar variables
% steps    = 0;
% msatc   = 0;
% sumAVG  = 0;
% sumSD   = 0;
%f_hist = [fBest,nan(1,M)];

%% Main process
tic,
for steps = 1 : M %+ 1 %while steps <= M && msatc < msat,
    
    Nest_ = Nest; fNest_ = Fitness;
    
    % Get a cuckoo randomly by L'evy flight
    nu = (randn(Na,Nd)*sigma)./(abs(randn(Na,Nd)).^(1/beta));
    Nest = Nest + alpha*randn(Na,Nd).*nu.*(Nest - repmat(BestNest,Na,1));

    % Check if the particle is in search space
    if unconst == false
        check = Nest < bnd_1; Nest = ~check.*Nest + check.*bnd_1;
        check = Nest > bnd_2; Nest = ~check.*Nest + check.*bnd_2;
     end

    % Evaluate objective function in the new position for each particle
    %fNest      = evaluate_function(fObj,Nest);

    % Find the best position for each particle
    %Nest      = get_Pi(fNest_ < Fitness,Nest,Nest);
    %Fitness   = min(fNest_,Fitness);

    % A fraction (pa) of worse nests are abandoned and new ones are built
    Nest = Nest + rand(Na,Nd).*(Nest(randperm(Na),:) - Nest(randperm(Na),:)).*...
        double(rand(Na,Nd) < pa);

     % Check if the particle is in search space
     if unconst == false
        check = Nest < bnd_1; Nest = ~check.*Nest + check.*bnd_1;
        check = Nest > bnd_2; Nest = ~check.*Nest + check.*bnd_2;
     end

    % Evaluate objective function in the new position for each particle
    Fitness      = evaluate_function(fObj,Nest);

    % Found the best position for each particle
    Nest      = get_Pi(Fitness < fNest_,Nest,Nest_);
    Fitness   = min(fNest_,Fitness);

    % Find the Nest best position (initial step)
    %f_hist(steps)    = fBest;
    [fBest,g] = min(Fitness); BestNest = Nest(g,:);

    % Statistical
    %sumAVG      = sumAVG + fBest;
    %sumSD       = sumSD + fBest^2;
    %currAVG     = sumAVG/steps;
    %currSD      = eps1*sqrt(sumSD/steps - currAVG^2);


    %sortedFitness = sort(Fitness);
    %currAVG     = mean(sortedFitness(1:10));
    %currSD      = eps1*std(sortedFitness(1:10));

    % Statistics related to dispersion of the particles
    %radii      = sqrt(sum((repmat(BestNest,Na,1) - Nest).^2,2));
    %radious     = max(radii);

    % Stop criteria
    %if fBest_ == fBest,% && (abs(fBest) < abs(currAVG + currSD)), %|| ... %&& radious < eps1
            %fBest_ = fBest,
    %    msatc = msatc + 1;
    %else
    %    msatc = 0;
    %end
    %if mod(steps,10) == 0,
    %        fprintf('%d, %.6g :: %s, [msat = %d]\n',steps,fBest,...
    %            sprintf('%.6g ',BestNest),msatc);
    % end
    %fprintf('fBest = %.4g at t = %d, x = (%s)\n',fBest,steps,sprintf('%.4g ',BestNest));
    %fprintf('fBest:%.6g, steps:%d, avg+sd:%.6g, sat:%d \n',fBest, steps, currAVG + currSD, msatc)

    % Plot objective function value evolution                    [toDelete]
 %        topl        = [topl; [steps,fBest,currAVG,currSD]];
 %        topl1       = [topl1; [steps,fBest,radious,0]];

    % Plot objective function value evolution
    % plot_things(topl,topl1,Nest,BestNest,[1 steps+1;1 steps+1;bnd(1,:)],[nan nan;0 eps2*1.5;bnd(2,:)],eps2,fBest,steps);

end
t       = toc;

%if steps >= M,   outmsg = 0; else outmsg = 1; end

details = struct('time',t,'fevs',(Na + 1)*(steps-1),'steps',steps-1);%,...
%     'outmsg',outmsg);%,'favg',currAVG,'fstd',currSD);
end
