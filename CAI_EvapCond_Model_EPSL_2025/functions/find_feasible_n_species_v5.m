function [ n_species_feasible, rel_violation ] = ...
    find_feasible_n_species_v4_1e10multiplier_20251018( formula_matrix, b_elements )
% FIND_FEASIBLE_N_SPECIES Given a formula matrix and element abundance
% vector (b_elements), this function finds a feasible species abundance
% vector (n_species_feasible). 
%
% The output 'n_species_feasible' is a species abundance vector that meets
% two criteria: (1) all species abundances are greater than zero, and
% (2) mass balance constraints are satisfied
%
% The function accomplishes this by using the MATLAB function linprog to
% solve a simple optimization problem - a feasibility problem 
% Relevant notes - see section 5.12 in these notes: 
% " https://www.stat.cmu.edu/~ryantibs/convexopt-F15/scribes/15-barr-method-scribed.pdf "
% (but, to be clear, this function mostly just relies on MATLAB's 'linprog',
% it doesn't follow the notes to do so)


%%% Scaling

multiplier = 1/min(b_elements);
% multiplier = 1/median(b_elements);

b_scaled = multiplier * b_elements;


%%% Create sizing variables

M = size(formula_matrix,1); % # of elements
N = size(formula_matrix,2); % # of species

sz_m = [M,1];   % size of vector of length M (# of elements)
sz_n = [N,1];   % size of vector of length N (# of species)


%%% Construct coefficient vector for objective function
f = [ zeros(sz_n); 1];


%%% Construct equality and inequality constraint matrices and vectors

% Construct equality constraint matrix
Aeq = [ formula_matrix, zeros(sz_m) ];

% Construct equality constraint vector
beq = b_scaled;

% Construct inequality constraint matrix
Ain = [ -1 * [ eye(N), ones(sz_n) ]; transpose(zeros(sz_n)), 1];

% Construct inequality contraint vector
bin = zeros(N+1,1);


%%% Employ MATLAB function 'linprog' to solve for feasible point

my_options = optimoptions(@linprog,'Display','off');


% Solve for set of values that corresponds to a feasible point

x = linprog(f, Ain, bin, Aeq, beq, [], [], my_options);


%%% Create vector with feasible values for output

% main output - feasible species abundance vector
n_species_feasible = x(1:N);

% Restore to original scale
n_species_feasible = n_species_feasible ./ multiplier;
b_elements = b_scaled ./ multiplier;

% Relative violation of mass balance constraints
rel_violation = ( formula_matrix * n_species_feasible - b_elements ) ...
    ./ b_elements;

end

