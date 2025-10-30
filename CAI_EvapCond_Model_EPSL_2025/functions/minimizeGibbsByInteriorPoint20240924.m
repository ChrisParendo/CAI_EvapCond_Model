function [n_species, gammas_sacm, mu_species, log_act_species, zeta_vector] = ...
    minimizeGibbsByInteriorPoint20240924( ...
    n_species_0, formula_matrix, b_elements, ...
    mu_0_species, R, T, P, i_gas, i_sacm, i_sol, ...
    W_sacm, Q_m_sacm, tau_sequence, is_sacmSpeciesIncluded, ...
    phaseIdentifier, gamma_species, is_traceSpecies)

%   minimizeGibbsByInteriorPoint20240924
%
%   This function is adapted from min_gibbs_nonstoichiometric_direct_2022x04x25 
%   (which retains obsolete lines of code showing different criteria for
%   convergence that were previously used)
%
%   This function minimizes the Gibbs energy of the system.  It employs 
%   a Newton-Raphson interior-point method (e.g., Leal et al., 2017)
%   and does not presume in advance which phases are stable / unstable. 
%
%   The primary output is the species abundance vector (n_species)
%   corresponding to the system with Gibbs energy minimized — i.e., the
%   system at equilibrium

%% Constants and Setup

% Factor to reduce step size if conditions are not met
step_size_reducer = 1/2;
% M -> # of elements, N -> # of species
[M,N] = size(formula_matrix);
% Size of species abundance vector
sz = size(n_species_0);

%% Initialize values — values correspond to zeroth iteration
n_species = n_species_0;    % inital species abundance vector
n_species_total = sum(n_species);   % sum of species abundances

% Calculate initial chemical potential vector (mu_species) and Jacobian of
% chemical potential vector (J_mu_species).  [Aside: Jacobian of chemical
% potential vector is equivalent to Hessian of Gibbs of system.]

[mu_species, J_mu_species, gammas_sacm, ~, ~, ~, log_act_species] = ...
    constructChemPotentialVectorAndJacobian20240924(sz, n_species, ...
    i_gas, i_sacm, i_sol, mu_0_species, R, T, P, W_sacm, Q_m_sacm, N, ...
    is_sacmSpeciesIncluded, phaseIdentifier, gamma_species, is_traceSpecies);

total_iteration_counter = 0; % initialize counter

%% Main Loop - Iterate Over Tau Sequence
for outer_iteration = 1:length(tau_sequence)

    tau = tau_sequence(outer_iteration);  % Tau for current iteration

    iteration_counter = 0;  % initialize counter for inner loop

    % Calculate KKT multiplier (zeta) vector
    zeta_vector = tau ./ n_species;   

    % Calculate Lagrange multiplier (lambda) vector
    lambda_vector = ...
        transpose(formula_matrix) \ ( zeta_vector - mu_species );       
    % [Aside: above expression shown on p.8 of my notes titled "Problem: 
    % min g(n) given constraints"]

    % Construct (reduced size) KKT matrix - corresp. to eqns for n and lambda
    kkt_matrix = ...
        [ J_mu_species + diag(n_species)^(-1) * diag(zeta_vector), ...
        transpose(formula_matrix); ...
        formula_matrix, zeros(M) ];

    % Construct gradient of Lagrangian function
    kkt_vector = ...
        [ mu_species + transpose(formula_matrix) * lambda_vector ...
        - tau * diag(n_species)^(-1) * ones(length(n_species),1) ; ...
        formula_matrix * n_species - b_elements ];    

    % Calculate delta vector(s):
    % n(iteration+1) - n(iteration) for all species, and 
    % lambda(iteration+1) - lambda(iteration) for all lagrange multipliers
    delta_n_and_lambda = kkt_matrix \ ( - kkt_vector );
    % [Aside: delta_n = delta_n_and_lambda(1:N);]

    % delta lambda
    delta_lambda = delta_n_and_lambda(N+1:end);

    % delta zeta
    delta_zeta = ( - zeta_vector ) + diag(n_species)^(-1) * ...
        ( tau - diag(zeta_vector) * delta_n_and_lambda(1:N) ) ;

    % Calculate gibbs energy of system (extensive)
    gibbs  = transpose(n_species) * mu_species;

    % Obsolete: molar_gibbs = gibbs / n_species_total; % Calculate molar gibbs of system


    %% Inner Loop: Iterate until convergence:
    max_iterations = 50; % Maximum allowed iterations for each Tau value

    while max( abs(delta_n_and_lambda(1:N)) ./ n_species ) > 1e-5 && ...
            iteration_counter < max_iterations

        total_iteration_counter = total_iteration_counter + 1;

        % count iterations
        iteration_counter = iteration_counter + 1;
        % inner_iteration = iteration_counter;

        % Provisional step-size and species abundance vector
        w = 1;
        n_species_prov = n_species + w * delta_n_and_lambda(1:N);

        % Condition #1
        % Check for negative abundances and reduce step-size if necessary
        while any( n_species_prov < 0 )
            w = step_size_reducer * w;
            n_species_prov = ...
                n_species + w * delta_n_and_lambda(1:N);
        end
                                
        % Similar but more strict than above condition (Added May 12, 2023)
        while any( n_species_prov.*(zeta_vector+w*delta_zeta) < tau/50)
            w = step_size_reducer * w;
            n_species_prov = ...
                n_species + w * delta_n_and_lambda(1:N); 
        end
    
        % Provisional total moles
        n_species_total_prov = sum(n_species_prov);

        % Provisional chemical potential vector (mu_species) and Jacobian
        [mu_species_prov, J_mu_species_prov, gammas_sacm, ~, ~, ~, log_act_species ] = ...
            constructChemPotentialVectorAndJacobian20240924( ...
            sz, n_species_prov, i_gas, i_sacm, i_sol, ...
            mu_0_species, R, T, P, W_sacm, Q_m_sacm, N, ...
            is_sacmSpeciesIncluded, phaseIdentifier, gamma_species, is_traceSpecies );

        % Provisional gibbs of system
        gibbs_prov = transpose(n_species_prov) * mu_species_prov;

        % Provisional Lagrange multiplier (lambda) vector
        lambda_vector_prov = lambda_vector + w * delta_lambda;

        % Condition #2: check that gibbs has decreased —
        % Check that Gibbs of system has decreased and reduce step-size if
        % necessary

        % merit function from prior values (reference)
        merit_ref = gibbs + (transpose(formula_matrix * n_species - b_elements) ...
            * lambda_vector) - (transpose(tau * ones(N,1)) * log(n_species));

        % merit function from provisional values
        merit_prov = gibbs_prov + (transpose(formula_matrix * n_species_prov - b_elements) ...
            * lambda_vector_prov) - (transpose(tau * ones(N,1)) * log(n_species_prov));

        % [Aside: objective_importance = abs(gibbs) / ( abs(gibbs) + abs(
        % (transpose(formula_matrix * n_species - b_elements) * lambda_vector) - 
        % (transpose(tau * ones(N,1)) * log(n_species)) ) );]
        % [Aside: objective_importance = gibbs / merit_ref]

        while merit_prov > merit_ref  && ...
                abs((merit_prov - merit_ref)/merit_ref) > 1e14 % added latter part 
                % so that floating point errors don't cause this loop to be engaged

        w = step_size_reducer * w;

        n_species_prov = n_species + w * delta_n_and_lambda(1:N);

        n_species_total_prov = sum(n_species_prov);

        [mu_species_prov, J_mu_species_prov, gammas_sacm, ~, ~, ~, log_act_species ] = ...
            constructChemPotentialVectorAndJacobian20240924( ...
            sz, n_species_prov, i_gas, i_sacm, i_sol, ...
            mu_0_species, R, T, P, W_sacm, Q_m_sacm, N, ...
            is_sacmSpeciesIncluded, phaseIdentifier, gamma_species, is_traceSpecies );

        gibbs_prov = transpose(n_species_prov) * mu_species_prov;

        lambda_vector_prov = lambda_vector + w * delta_lambda;
        merit_prov = gibbs_prov + (transpose(formula_matrix * n_species_prov - b_elements) ...
            * lambda_vector_prov) - (transpose(tau * ones(N,1)) * log(n_species_prov));

        end

        % Update all values -> no longer provisional

        n_species = n_species_prov;

        lambda_vector = lambda_vector + w * delta_n_and_lambda(N+1:end);
        % Alternative?: lambda_vector = transpose(formula_matrix) \ ( - mu_species );  

        zeta_vector = zeta_vector + w * delta_zeta;

        n_species_total = n_species_total_prov;

        mu_species = mu_species_prov;

        J_mu_species = J_mu_species_prov;

        kkt_matrix(1:N,1:N) = J_mu_species + diag(n_species)^(-1) * diag(zeta_vector);

        kkt_vector = ...
            [ mu_species + transpose(formula_matrix) * lambda_vector ...
            - tau * diag(n_species)^(-1) * ones(length(n_species),1) ; ...
            formula_matrix * n_species - b_elements ];

        delta_n_and_lambda = kkt_matrix \ ( - kkt_vector );

        delta_lambda = delta_n_and_lambda(N+1:end);

        delta_zeta = ( - zeta_vector ) + diag(n_species)^(-1) * ...
            ( tau - diag(zeta_vector) * delta_n_and_lambda(1:N) ) ;
    
        gibbs = gibbs_prov;
    
        % Obsolete: molar_gibbs = gibbs / n_species_total;

    end  

end

total_iteration_counter; % may sometimes want to display total iterations

end

