function [ n_sat_vector, mu_sat_vector, V_output, rxnRates, ...
    evaporationFluxes_all, productionFluxes] = ...
    calculate_vapor_pressures_fluxes_prescribed_condensate_20240915( ...
    simulation_state )
% CALCULATE_VAPOR_PRESSURES_FLUXES_PRESCRIBED_CONDENSATE
%
% This function is adapted from calcVaporPressuresAndFluxes5
%
%   This function calculates equilibrium (vapor) pressures and evaporative 
%   fluxes.  The calculation consists of solving a system of equations 
%   comprised of (1) a set of equations for the free energies of reactions
%   for the relevant chemical reactions *one equation for each reaction*
%   and (2) a set of equations for the constraint that, at steady state,
%   the rate of production of gas species by chemical reactions must equal 
%   the rate at which same gas species are transported to the ambient gas
%   reservoir **one equation for each gas species*.


%% INPUT PARAMETERS

% Assigned Inputs

% constant logicals
is_css = simulation_state.constant_logicals.is_css;
is_gss = simulation_state.constant_logicals.is_gss;
is_sacm = simulation_state.constant_logicals.is_sacm;
is_traceSpecies = simulation_state.constant_logicals.is_traceSpecies;

% dynamic logicals
is_stable = simulation_state.dynamic_logicals.is_stable;

% constant parameters
formula_matrix = simulation_state.constant_params.formulaMatrix;
R = simulation_state.constant_params.R;
stickingCoefficients = simulation_state.constant_params.stickingCoeffs;
molecularMasses_kg_mol = simulation_state.constant_params.molMasses_kg_mol;
substanceIdentifier = simulation_state.constant_params.substanceIdentifier;

% dynamic parameters
n_species = simulation_state.dynamic_params.n_species;
mu_0_species = simulation_state.dynamic_params.mu_0_species;
mu_species = simulation_state.dynamic_params.mu_species;

T = simulation_state.dynamic_params.temperature_K;
volume_gss = simulation_state.dynamic_params.volumeGas_cL;
particleArea = simulation_state.dynamic_params.particleArea_m2;

% numberical tolerances
tolerance_energyRxn = ...
    simulation_state.numerical_tolerances.tolerance_energyRxn;
    % tolerance_relativeFlux = ...
    %     simulation_state.numerical_tolerances.tolerance_relativeFlux;

% Hardcoded Inputs
step_size_reducer = 0.5;


%% SET UP CALCULATION

s_gss = stickingCoefficients(is_gss);   % MOVE THIS TO BETTER LOCATION
m_gss = molecularMasses_kg_mol(is_gss); % MOVE THIS TO BETTER LOCATION

%%% Create Indices

traceSubstanceIdentifier = substanceIdentifier(is_traceSpecies);
uniqueTraceSubstanceIdentifier = unique(traceSubstanceIdentifier);
numUniqueTraceSubstances = length(uniqueTraceSubstanceIdentifier);
is_majorSpecies = ~ is_traceSpecies;

% Designate for which species chemical potentials will be held constant:

% Create logical for fixed (held constant) major species
is_fixedMajorSpecies = is_majorSpecies & is_stable & is_css ;
if any( is_fixedMajorSpecies & is_sacm )
    is_fixedMajorSpecies = is_fixedMajorSpecies & is_sacm ;
end

% Create logical for fixed (held constant) trace species
is_fixedTraceSpecies = false(size(n_species)); 
for i = 1:numUniqueTraceSubstances 
    
    idx_thisSubstance = uniqueTraceSubstanceIdentifier(i);
    is_thisSubstance = idx_thisSubstance == substanceIdentifier;

    % if trace substance is stable in more than 1 condensed species, then
    % select and fix only one of the condensed species

    % select one (first) species consisting of this trace substance
    idx_chosenSpeciesOfSubstance = ...    
        find( (is_thisSubstance & is_stable & is_css), 1 );
    
    is_fixedTraceSpecies(idx_chosenSpeciesOfSubstance) = true;
end

% Create logical for fixed (held constant) species - major and trace
is_fixed = is_fixedMajorSpecies | is_fixedTraceSpecies;


% Create list indices: idx1 indices correspond to input species vector
idx1_gss = find(is_gss);  
idx1_fixed = find(is_fixed);
idx1 = [idx1_gss; idx1_fixed];


%%% Create formula matrix and stoichiometric matrix for select system

% Create formula matrix for selected system 'A'
A = formula_matrix(:,idx1);

% Create stoichiometric matrix
V = null(A,'rational');
% Commented out NOV 13, 2023 ] loop to enforce sign convention that 
% evaporating condensed species (as well as H2) is '-'
    % for j = 1:size(V,2)       % [PROBABLY NOT NEEDED!, 
    %     if V(1,j) > 0         % [PROBABLY NOT NEEDED! ]
    %        V(:,j) = - V(:,j); % [PROBABLY NOT NEEDED! ]
    %     end                   % [PROBABLY NOT NEEDED! ]
    % end                       % [PROBABLY NOT NEEDED! ]

% Number of reactions
num_rxns = size(V,2);


%%% Split stoichiometric matrix and chemical potential vector into parts...
% for free (gas) and fixed (condensed) species

% Create list indices: idx2 indices correspond to species of select system

num_gss = sum(is_gss);
num_fixed = sum(is_fixed);

idx2_gss = 1:num_gss;
idx2_fixed = (num_gss+1):(num_gss+num_fixed);

% Split stoichiometric matrix into free (gas) and fixed (condensed) parts:
V_gss = V(idx2_gss,:);
V_fixed = V(idx2_fixed,:);

% Split standard chemical potential vector into part for free(gas) species:
mu_0_species_gss = mu_0_species(idx1_gss);

%% PERFORM CALCULATION

%%% Initialization (i.e., iteration #0)

iteration = 0;

n_species_gss = n_species(idx1_gss);
n_gss_ambient = n_species_gss;

n_total_gss = sum(n_species_gss);

% chemical potential vectors
mu_species_gss = mu_species(is_gss);
mu_species_fixed = mu_species(is_fixed);

% Free energy of rxns: Below, mu and J_mu expressions are for constant V (not P)
energyRxns_gssPart = transpose(V_gss) * mu_species_gss;
energyRxns_fixedPart = transpose(V_fixed) * mu_species_fixed;
energyRxns = energyRxns_gssPart + energyRxns_fixedPart;

% Jacobian of mu w.r.t. n species (for constant gas volume)
J_mu_wrt_n_gss = R*T*diag(1./n_species_gss); 

% Jacobian of 'energy of rxns' w.r.t. n species (for constant gas volume)
J_energyRxns_wrt_n = transpose(V_gss)*J_mu_wrt_n_gss;

% Jacobian of 'energy of rxns' w.r.t. psi (where psi is extent of rxn ...
% 'xi' per unit time 't'; i.e., psi = xi / t; unit of time is here 1 s)
J_energyRxns_wrt_psi = zeros(num_rxns);

% Net flux into / out from (infinitesimal) atmosphere around particle

z1 = particleArea * (R*T/(2*pi))^(1/2) / (volume_gss*1e-5);    % a scalar         %%%%%% volume used here should be m3 !!!!!!!!!!!!!!!!!!!
z2 = s_gss./(m_gss.^(1/2));                                 % a vector
z3 = z1 * z2;                                               % a vector

% psi_vector = zeros(num_rxns,1); % use 0 as initial guess?

evaporationFluxes = z3 .* (n_species_gss - n_gss_ambient);

psi_vector = V_gss \ evaporationFluxes;

productionFluxes = V_gss * psi_vector;

netFluxesAtmosphere = productionFluxes - evaporationFluxes;

% Jacobian of 'net fluxes' w.r.t. n species (for constant gas volume)
J_netFluxesAtmosphere_wrt_n = -diag(z3); % will be constant

% Jacobian of 'net fluxes' w.r.t. psi 
J_netFluxesAtmosphere_wrt_psi = V_gss; % will be constant

% Residuals Function
residualsFun = [ energyRxns; netFluxesAtmosphere];

% Jacobian of Residuals Function
J_residualsFun = [ J_energyRxns_wrt_n, J_energyRxns_wrt_psi; ...
    J_netFluxesAtmosphere_wrt_n, J_netFluxesAtmosphere_wrt_psi];

% Calculate step direction and magnitude
delta_n_and_psi = J_residualsFun \ ( - residualsFun);


%%% Iterate until convergence:

% convergence criteria
shouldContinue1 = max( abs( energyRxns/(R*T) ) ) > tolerance_energyRxn;
% shouldContinue2 =  max( abs( netFluxesAtmosphere./n_species_gss ) ) > ...
%     tolerance_relativeFlux;
% shouldContinue2 = max(abs(productionFluxes./evaporationFluxes-1)) > ...
%         tolerance_relativeFlux;
%     shouldExclude = isnan(evaporationFluxes) | (evaporationFluxes == 0);
%     shouldContinue2 = max(abs(productionFluxes(~shouldExclude)./evaporationFluxes(~shouldExclude)-1)) > ...
%         tolerance_relativeFlux;
    shouldContinue2 = false;

while shouldContinue1 || shouldContinue2 

    iteration = iteration + 1;

    w = 1; % reset step-size (w) to 1

    % Update species abundance and psi vectors (provisionally)
    n_species_gss_prov = n_species_gss + w * delta_n_and_psi(1:num_gss);
    psi_vector_prov = psi_vector + w * delta_n_and_psi(num_gss+1:end);

    % Enforce Condition #1: Nonnegativity constraint for each time step:
    while any(n_species_gss_prov < 0)
        w = step_size_reducer * w;
    n_species_gss_prov = n_species_gss + w * delta_n_and_psi(1:num_gss);
    psi_vector_prov = psi_vector_prov + w * delta_n_and_psi(num_gss+1:end);    
    end

    % Update all values - no longer provisional:

    n_species_gss = n_species_gss_prov;
    psi_vector = psi_vector_prov;

    n_total_gss = sum(n_species_gss);

    rho = n_total_gss / volume_gss;                                             %%%%%% volume used here should be cL !!!!!!!!!!!!!!!!!!!
    P = R*T*rho;
    mu_species_gss = ...
        mu_0_species_gss + R*T*log(n_species_gss./n_total_gss*P);

    energyRxns_gssPart = transpose(V_gss)*mu_species_gss; 
    energyRxns = energyRxns_gssPart + energyRxns_fixedPart;

    % Update - commented out lines where values are constant

    J_mu_wrt_n_gss = R*T*diag(1./n_species_gss); 
    J_energyRxns_wrt_n = transpose(V_gss)*J_mu_wrt_n_gss;
    % J_energyRxns_wrt_psi = zeros(num_rxns);
    
    productionFluxes = V_gss * psi_vector;
    % z1 = particleArea * (R*T/(2*pi))^(1/2) / particleVolume;    % a scalar
    % z2 = s_gss./(m_gss.^(1/2));                                 % a vector
    % z3 = z1 * z2;                                               % a vector
    evaporationFluxes = z3 .* (n_species_gss - n_gss_ambient);
    netFluxesAtmosphere = productionFluxes - evaporationFluxes;
    
    % J_netFluxesAtmosphere_wrt_n = -diag(z3); % will be constant
    % J_netFluxesAtmosphere_wrt_psi = V_gss; % will be constant
    
    residualsFun = [ energyRxns; netFluxesAtmosphere];
    
    J_residualsFun(1:num_rxns,1:num_gss) = J_energyRxns_wrt_n;
    
    % Update

    delta_n_and_psi = J_residualsFun \ ( - residualsFun);

    % convergence criteria
    shouldContinue1 = max( abs( energyRxns/(R*T) ) ) > tolerance_energyRxn;
%     shouldContinue2 =  max( abs( netFluxesAtmosphere./n_species_gss ) ) ...
%         > tolerance_relativeFlux;
%     shouldContinue2 = max( abs( netFluxesAtmosphere./z3) )> ...
%         tolerance_relativeFlux;

% % %      productionFluxes./evaporationFluxes

%     max(abs(productionFluxes./evaporationFluxes-1))
%     shouldExclude = isnan(evaporationFluxes) | (evaporationFluxes == 0);
%     shouldContinue2 = max(abs(productionFluxes(~shouldExclude)./evaporationFluxes(~shouldExclude)-1)) > ...
%         tolerance_relativeFlux
        shouldContinue2 = false;


% abs( productionFluxes(~shouldExclude)./evaporationFluxes(~shouldExclude)-1 )
% evaporationFluxes(~shouldExclude)

end


%% Organize Output

sz = size(n_species);
num_species = length(n_species);
num_rxns = size(V,2);

n_sat_vector = nan(sz);
n_sat_vector(is_gss) = n_species_gss;

mu_sat_vector = nan(sz);
mu_sat_vector(is_gss) = mu_species_gss;

V_output = zeros(num_species,num_rxns);
V_output(is_gss,:) = V_gss;
V_output(is_fixed,:) = V_fixed;


evaporationFluxes_all = nan(sz);
evaporationFluxes_all(is_gss) = evaporationFluxes(is_gss);

rxnRates = psi_vector;

productionFluxes_all = V_output * psi_vector;
productionFluxes = productionFluxes_all;


end %% End of function

