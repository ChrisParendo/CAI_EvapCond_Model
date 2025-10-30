function [ mu_species, J_mu_species, gamma_sacm, xxxx_sacm, ...
    d_dn_x_sacm, d_dn_xxxx_sacm, log_act_species ] ...
    = constructChemPotentialVectorAndJacobian20240928( equilibrium_params, sz, N)


% This function (adapted from calc_gamma_mu_Jmu_2023x09x26) calculates a
% vector with the chemical potentials (mu's) of all species and the
% Jacobian of the chemical potentials for the system

% 2024/09/28: Adds lines for treatment of K2O(melt)
% 2024/09/28: Structure as input

%% Unpack parameters from structure                                        % 2024/09/28

% Unpack parameters for easier access
n_species = equilibrium_params.n_species;                                 
gamma_species = equilibrium_params.gamma_species;

% Unpack indices
i_gas = equilibrium_params.indices.i_gas;
i_sacm = equilibrium_params.indices.i_sacm;
i_sol = equilibrium_params.indices.i_sol;

% Unpack logical flags
is_sacmSpeciesIncluded = equilibrium_params.logicals.is_sacmSpeciesIncluded;
is_traceSpecies = equilibrium_params.logicals.is_traceSpecies;
is_K2O_melt = equilibrium_params.logicals.is_K2O_melt;

% Unpack thermodynamic parameters
R = equilibrium_params.thermodynamics.R;
T = equilibrium_params.thermodynamics.T;
P = equilibrium_params.thermodynamics.P;
mu_0_species = equilibrium_params.thermodynamics.mu_0_species;

% Unpack identifiers
phase_identifier = equilibrium_params.identifiers.phase_identifier;

% Unpack melt parameters
W_sacm = equilibrium_params.cmas_melt.W_sacm;
Q_m_sacm = equilibrium_params.cmas_melt.Q_m_sacm;

%% Prepare/Preallocate

x_in_phase = nan(sz); % Mole fraction of species in its respective phase
act_species = nan(sz);  % Activity of species
log_act_species = nan(sz);  % Natural log of activity of species
mu_species = nan(sz);   % Chemical potential of species

J_mu_species = zeros(N); % Jacobian of chemical potentials
                         % safer to preallocate NaN and replace w/ zeros?


%% Calculate values for gaseous species: ideal gas model

% Calculate mole fractions, activities, and chemical potentials

x_in_phase(i_gas) = n_species(i_gas) ./ sum(n_species(i_gas));
act_species(i_gas) = x_in_phase(i_gas) * P;
log_act_species(i_gas) = log(act_species(i_gas));

mu_species(i_gas) = mu_0_species(i_gas) + R*T*log_act_species(i_gas);

% Calculate the Jacobian of mu w.r.t. n

num_gas = length(i_gas);

J_mu_species(i_gas,i_gas) = R*T*(repmat(-1/sum(n_species(i_gas)), ...
    num_gas, num_gas) + diag(1./n_species(i_gas)));


%% Calculate values for melt species: Berman1983 CMAS melt model

% if melt species are permitted in system, calculate gammas, mu, and J_mu

if any(i_sacm)

x_in_phase(i_sacm) = n_species(i_sacm)./sum(n_species(i_sacm));

X_sacm_input = nan(4,1);
X_sacm_input(~is_sacmSpeciesIncluded) = 0;
X_sacm_input(is_sacmSpeciesIncluded) = x_in_phase(i_sacm);

threshold = 1e-9; % prevents mole fractions from being too close to zero
isLow = X_sacm_input < threshold;
X_sacm_input(isLow) = threshold;
X_sacm_input = X_sacm_input / sum(X_sacm_input); % renormalize to unity

% run functions that together determine activity coefficients (gammas) and 
% J_mu for CMAS melt:

[d_dn_x_sacm, xxxx_sacm, d_dn_xxxx_sacm] = ...
    calc_sacm_JinvX_XXXX_JXXXX( sum(n_species(i_sacm)), ...
        X_sacm_input(1),X_sacm_input(2),X_sacm_input(3),X_sacm_input(4));

[gamma_sacm, ~, d_dn_mu_sacm] = ...
    calc_sacm_gamma_vector_d_dn_mu_matrix_sacm_FIXING( ...
        W_sacm, Q_m_sacm, X_sacm_input, d_dn_x_sacm, xxxx_sacm, ...
        d_dn_xxxx_sacm, R,T, sum(n_species(i_sacm)) );
    
gamma_species(i_sacm) = gamma_sacm(is_sacmSpeciesIncluded);

act_species(i_sacm) = gamma_species(i_sacm) .* x_in_phase(i_sacm);
log_act_species(i_sacm) = log(act_species(i_sacm));

mu_species(i_sacm) = mu_0_species(i_sacm) + R*T*log_act_species(i_sacm);

J_mu_species(i_sacm,i_sacm) = ...
    d_dn_mu_sacm(is_sacmSpeciesIncluded,is_sacmSpeciesIncluded);

else % if no melt species are permitted/included in system

% assign NaN to outputs of function that pertain to melt
gamma_sacm = NaN;
xxxx_sacm = NaN;
d_dn_x_sacm = NaN;
d_dn_xxxx_sacm = NaN;

end


%% Calculate values for melt trace species: Constant activity coefficient model
% (can clean this up later - maybe double-check too, was done quickly)

% Create "is Melt?" logical
falseVector = false(N,1);
isGas = falseVector;
isGas(i_gas) = true;
isSol = falseVector;
isSol(i_sol) = true;
isMelt = ~ ( isGas | isSol );

isTraceInMelt = isMelt & is_traceSpecies;

% Calculate total moles in melt
nt_inMelt = sum(n_species(isMelt));

% Calculate mole fractions of *trace* species in melt
x_in_phase(isTraceInMelt) = n_species(isTraceInMelt) ./ nt_inMelt;

%%% Calculate activity coefficient for K2O in melt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (Added Sept. 28, 2024)

    if any(is_K2O_melt)
        sticking_coefficient = 0.1; % Temporarily hardcoded
        mole_fraction_K2O_melt = x_in_phase(is_K2O_melt);
        [activity_coefficient_K2O_melt, ~] = ... 
            calculateActivityCoefficientK2Omelt20240928(T, sticking_coefficient, ...
            mole_fraction_K2O_melt);
        gamma_species(is_K2O_melt) = activity_coefficient_K2O_melt;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate activities of *trace* species in melt
act_species(isTraceInMelt) = ...
    gamma_species(isTraceInMelt) .* x_in_phase(isTraceInMelt);
 
log_act_species(isTraceInMelt) = log(act_species(isTraceInMelt));

% Calculate chemical potential (mu)
mu_species(isTraceInMelt) = mu_0_species(isTraceInMelt) + R*T*log_act_species(isTraceInMelt);

% Calculate and assign entries to J_mu
num_traceSpeciesInMelt = sum(isTraceInMelt);

J_mu_species(isTraceInMelt,isTraceInMelt) = ...
    R*T*(repmat(-1/nt_inMelt,num_traceSpeciesInMelt,num_traceSpeciesInMelt) ...
    + diag(1./n_species(isTraceInMelt)));
J_mu_species(isTraceInMelt,i_sacm) = -1/nt_inMelt;
J_mu_species(i_sacm,isTraceInMelt) = -1/nt_inMelt;  % Think this is right, but should double-check.


%% Calculate/Assign Values for Solid Species

% for each solid phase, calculate mu and the entry/entries of J_mu for all 
% species in that solid phase

uniqueSolidPhases = unique(phase_identifier(i_sol));

numSolidPhases = length(uniqueSolidPhases);

for i = 1:numSolidPhases

    idx = uniqueSolidPhases(i);

    % logical indicating which species are part of this phase
    isPhase = idx == phase_identifier;

    % total moles in this phase
    nt_in_phase = sum(n_species(isPhase));

    % mole fractions of species that comprise this phase
    x_in_phase(isPhase) = n_species(isPhase) ./ sum(n_species(isPhase));

    % calculate activities ...
    % ... for trace species
    act_species(isPhase & is_traceSpecies) = ...
        gamma_species(isPhase & is_traceSpecies) ...
        .* x_in_phase(isPhase & is_traceSpecies);
    % ... for major species, treating them as pure phases
    act_species(isPhase & ~is_traceSpecies) = 1;

    log_act_species(isPhase) = log(act_species(isPhase));

    % calculate mu
    mu_species(isPhase) = ...
        mu_0_species(isPhase) + R*T*log_act_species(isPhase);

    % calculate and assign entries to J_mu

    num_speciesInPhase = sum(isPhase);

    J_mu_species(isPhase,isPhase) = ...
    R*T*(repmat(-1/nt_in_phase,num_speciesInPhase,num_speciesInPhase) ...
    + diag(1./n_species(isPhase)));

end


end % End of function