function [Xg,fg,details] = SSOA(fObj,bnd)

% Load Statistics' Package
%pkg load statistics

% Read parameters
Na      = 40;
theta   = pi/8;
rl      = 0.75;
M       = 1e12;

Tol     = 1e-12;
msat    = 1000;

% Dimensions
Nd      = size(bnd,1);

% Pre-allocate some variables
R       = ones(Nd);
I       = eye(Nd);

% Define a quick function to obtain f(X)
    function the_results = evaluate_function(the_function,the_positions)
        the_results   = nan(Na,1);
        for s = 1 : Na,
            the_results(s) = the_function(the_positions(s,:));
        end
    end

% Make the rotation matrix
cmbn    = combnk(1:Nd,2);   % Possible combinations by 2D plane
for ij = 1 : size(cmbn,1),
    i           = min(cmbn(ij,:));
    j           = max(cmbn(ij,:));

    R_aux       = eye(Nd);
    R_aux(i,i)  = cos(theta);
    R_aux(j,j)  = cos(theta);
    R_aux(i,j)  = -sin(theta);
    R_aux(j,i)  = sin(theta);

    if ij == 1, R = R.*R_aux;
    else        R = R*R_aux; end
end
S = R;

% Define the boundaries for each dimension
bnd         = [min(bnd,[],2) max(bnd,[],2)];
bnd_1       = repmat(bnd(:,1),1,Na);
bnd_2       = repmat(bnd(:,2),1,Na);

% Calculate the initial positions
X           = bnd_1 + rand(Nd,Na).*(bnd_2 - bnd_1);

% Fitness function
getFitness  = @(x) evaluate_function(fObj,x);

% Calculate the initial fitness values
fX          = getFitness(X');

% Found the best position to be rotation centre
[fg,ig] = min(fX); Xg = X(:,ig);

% Find initial value for max radii
maxradii    = max(sum((repmat(Xg,1,Na) - X).^2,2));

% Set auxiliar variables
steps    = 1;
msatc   = 0;

% Statistical variables
sumAVG  = 0;
sumSD   = 0;
%fv      = nan(M+1,1);
fv(1)   = fg;

%% Main process
tic,
while steps <= M && msatc < msat && maxradii > Tol
    
    for i = 1 : Na,
        % Update the position for each point
        r       = rl + (1 - rl)*rand;
        X(:,i)  = r*S*X(:,i) - rand(Nd,1).*((r*S - I)*Xg);%
    end

    % Check if the particle is inside the search space
    check = X < bnd_1; X = ~check.*X + check.*bnd_1;
    check = X > bnd_2; X = ~check.*X + check.*bnd_2;

    % Evaluate objective function in new positions
    fX          = getFitness(X');

    fg_ = fg;
    % Found the best position to be rotation center
    [fg,ig]     = min(fX);
    if fg < fg_, Xg = X(:,ig);
    else fg = fg_; end

    % Statistical block
    sumAVG      = sumAVG + fg;
    sumSD       = sumSD + fg^2;
    currAVG     = sumAVG/steps;
    currSD      = sqrt(sumSD/steps - currAVG^2);

    % Update step
    steps    = steps + 1;

    fv(steps) = fg;
    
    % Stop criteria:
    % 1. Stagnation
    msatc       = (msatc + 1)*double(fg == fg_); % abs(fBest - currAVG) < currSD && 
    % 2. Max radii population
    maxradii    = max(sum((repmat(Xg,1,Na) - X).^2,2));

end
t = toc; Xg = Xg';

% Final things
if steps >= M,   outmsg = 0; else outmsg = 1; end

details = struct('time',t,'fevs',steps*(Na + 1),'steps',steps,...
    'outmsg',outmsg,'favg',currAVG,'fstd',currSD,'historical',fv);
end
