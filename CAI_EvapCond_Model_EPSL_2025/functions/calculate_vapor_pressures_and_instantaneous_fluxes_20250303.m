function [ n_saturation, mu_saturation, V_output, rxnRates, ... 
    evaporationFluxes, productionFluxes, is_newlyStable, P_envelope ] = ...
    calculate_vapor_pressures_and_instantaneous_fluxes_20250303( ...
    simulation_state )
% CALCULATE_VAPOR_PRESSURES_AND_INSTANTANEOUS_FLUXES
%
% This function is adapted from calcPressuresFluxesInstant5
% 
% This function includes the function "calculate_vapor_pressure_fluxes...",
% which is the main part of the vapor pressure / flux calculation.
%
% This function runs "calculate_vapor_pressure_fluxes...", then checks to
% see if any additional condensed (solid) phases, besides those already in 
% the condensed system, are "newly stable".  That is, after calculating the
% vapor pressures of all gas species ("calculate_vapor_pressure_fluxes.."),
% this function checks if any additional phases will be thermodynamically
% inclined to condense.  If so, it adds that phase to the condensed system,
% and repeats the main calculation (again using 
% "calculate_vapor_pressure_fluxes..")
%
%
% Description of outputs:
% (1) n_saturation —> species abundances in gas that is in equilibrium with
% condsensed system
% (2) mu_saturation -> chemical potential of species in gas that is in
% equilibrium with condensed system 
% (3) V_output —> stoichiometric vectors for relevant
% evaporation-condensation reactions
% (4) evaporationFluxes
% (5) productionFluxes (should be equal to evaporation fluxes)
% (6) is_newlyStable —> logical indicating which condensed species, if any,
% are newly stabilized


%% INPUT PARAMETERS

% Unpack required inputs:

% constant logicals
is_cond = simulation_state.constant_logicals.is_css;
is_gas = simulation_state.constant_logicals.is_gss;
is_sol = simulation_state.constant_logicals.is_sol;
is_traceSpecies = simulation_state.constant_logicals.is_traceSpecies;

% dynamic logicals
is_stable = simulation_state.dynamic_logicals.is_stable;

% constant parameters
formulaMatrix = simulation_state.constant_params.formulaMatrix;

% dynamic parameters
mu_species = simulation_state.dynamic_params.mu_species;

% (Note also: Entire input structure will be passed once to ... 
% function "calculate_vapor_pressure_fluxes...")

%% CALCULATE VAPOR PRESSURES AND FLUXES FOR PRESCRIBED CONDENSED SYSTEM

[ n_saturation, mu_saturation, V_output, rxnRates, evaporationFluxes, ...
  productionFluxes, P_envelope ] = ...
  calculate_vapor_pressures_fluxes_prescribed_condensate_20250303( ...
  simulation_state);


%% CHECK FOR NEWLY STABLE CONDENSED SPECIES

% Create stoichiometric matrix, expressed as all possible
% evaporation-condensation reactions (Slightly inefficient to run this
% repeatedly rather than just once in main script, outside this function,
% but oh well)
[EvaporationStoichiometricMatrix] = ... 
    calcEvaporationCondensationMatrix_v2( formulaMatrix,is_gas, is_cond);


% Check for newly stable condensed species (only solids for now)

% Organize vector of chemical potentials for all species in system,
% neglecting trace species (as they do not determine if phases are stable
% or not)
mu_particleAndAtmosphere = nan(size(mu_species));
mu_particleAndAtmosphere(is_cond & ~is_traceSpecies ) = ...
    mu_species(is_cond & ~is_traceSpecies );
mu_particleAndAtmosphere(is_gas & ~is_traceSpecies ) = ...
    mu_saturation(is_gas & ~is_traceSpecies );

% Calculate energies of reactions:

% reducedEnergyRxns = transpose(EvaporationStoichiometricMatrix) * ...
%     mu_particleAndAtmosphere ./ (R*T);
energyRxns = nan(sum(is_cond),1); 
energyRxns(~is_traceSpecies(is_cond)) = ... 
   transpose(EvaporationStoichiometricMatrix( ...
   ~is_traceSpecies,~is_traceSpecies(is_cond))) * ...
   mu_particleAndAtmosphere(~is_traceSpecies); 

% identify where energies of reactions are such (>0) that condensation
% will occur:

is_newlyStableSolid = false(size(is_stable));
is_newlyStableSolid(is_cond) = (energyRxns > 0) & ...
    (~is_stable(is_cond)) & is_sol(is_cond);

% create "is_newlyStable" logical - main output of this function
is_newlyStable = is_newlyStableSolid;

% Make sure only one newly stable species is added
if sum(is_newlyStable) > 1
    % [~,idx] = max( is_newlyStable(is_cond) .* reducedEnergyRxns ) 
    idx = find( is_newlyStable(is_cond) ,1 ) ;
    is_newlyStable(is_cond) = ...
        transpose( 1:length(is_newlyStable(is_cond)) == idx ); 
end


%% IF THERE IS A 'NEWLY STABLE SPECIES', RECALCULATE PRESSURES/FLUXES

if any(is_newlyStable)

    % Create copy of simulation state
    simulation_state_new = simulation_state;

    % Add newly stable species to is_stable logical in new simulation state
    simulation_state_new.dynamic_logicals.is_stable = ...
        is_stable | is_newlyStable;

    % Calculate vapor pressures and fluxes for (updated) prescribed 
    % condensed system
    [ n_saturation, mu_saturation, V_output, rxnRates, ...
      evaporationFluxes, productionFluxes, P_envelope ] = ...
      calculate_vapor_pressures_fluxes_prescribed_condensate_20250303( ...
      simulation_state_new);

end


end %% End of function