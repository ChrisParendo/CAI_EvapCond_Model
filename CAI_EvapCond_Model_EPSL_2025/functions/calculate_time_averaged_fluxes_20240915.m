function [ productionFluxes_timeAveraged] = ...
    calculate_time_averaged_fluxes_20240915( simulation_state ) 

% CALCULATE_TIME_AVERAGED_FLUXES_20240915
%
% This function is adapted from calcPressuresFluxesTimeAveraged5
%                                                                           % Version 4 - limits loss of trace element species
% This function
%
% [ REMINDER: stoichiometricMatrix * rxnsRates = productionFluxes ]


%% INPUT PARAMETERS

% Unpack required inputs:

% constant logicals
is_traceSpecies = simulation_state.constant_logicals.is_traceSpecies;
is_cond = simulation_state.constant_logicals.is_css;

% dynamic logicals
is_stable = simulation_state.dynamic_logicals.is_stable;
is_newlyStable = simulation_state.dynamic_logicals.is_newlyStable;

% constant parameters
dt = simulation_state.constant_params.dt;

substanceIdentifier = simulation_state.constant_params.substanceIdentifier;

% dynamic parameters
n_species = simulation_state.dynamic_params.n_species;

productionRates = simulation_state.dynamic_params.productionFluxesInstant; % dynamic

% copy of simulation state structure
simulation_state_new = simulation_state;


%% Preparation

nan_vector = nan(size(n_species));

is_stable = is_stable  | is_newlyStable;


%% CALCULATE TIME UNTIL CONDENSED SPECIES/PHASES DISAPPEARS
% I.e., calculate time till condensed species/phases fully evaporate, and
% note that this calculation applies to major, not trace species

% Calculate consumption fluxes for condensed species:
consumptionRates = - productionRates;

% Create logical indicating which species are being consumed
isDecreasing = consumptionRates > 0;
is_decreasingSpecies = is_cond & is_stable & isDecreasing;

% Create vector giving time until condensed species will disappear,
% excluding trace species

timeToDisappear_vector = nan_vector;    % Preallocate with NaN

timeToDisappear_vector(is_decreasingSpecies & ~is_traceSpecies) = ...
    n_species(is_decreasingSpecies & ~is_traceSpecies) ./ ...
    consumptionRates(is_decreasingSpecies & ~is_traceSpecies);

% Calculate time until the first species/phase will disappear
[ timeToDisappearance, index_speciesDisappearing ] = ... 
    min( timeToDisappear_vector );


%% CALCULATE TIME-AVERAGED FLUXES FOR SPECIES — EXCLUDING TRACE SPECIES

%%% IF NO SPECIES/PHASE WILL DISAPPEAR OVER NEXT TIME STEP, THEN
%%% TIME-AVERAGED FLUXES SIMPLY EQUAL INSTANTANEOUS FLUXES:
if timeToDisappearance >= dt || isnan(timeToDisappearance)

    productionFluxes_timeAveraged = productionRates;

%%% OTHERWISE (I.E., ELSE), CARRY OUT CALCULATION THAT BREAKS UP TIME STEP
%%% INTO DISCRETE INTERVALS IN ORDER TO CALCULATE TIME-AVERAGED FLUXES
else

    i = 1;  % initialize time-interval counter

    % create vectors that save number of intervals, duration of
    % intervals, and fluxes within intervals:

    saved_intervalCount(i,1) = i;
    saved_time(i,1) = timeToDisappearance;
    saved_productionRates(:,i) = productionRates;

    % while a species is forecast to disappear prior to end of time step:
    while timeToDisappearance < dt  

        i = i + 1; % count time intervals

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Run function to to obtain fluxes for interval *after* a species 
        % disappears

        % Update/calculate species abundance vector for moment after a 
        % condensed species has disappeared:
        n_species(is_cond & is_stable) = n_species(is_cond & is_stable) ...
            + productionRates(is_cond & is_stable) * timeToDisappearance;

        % Update 'is_stable' logical to account for disappeared species
        is_stable(index_speciesDisappearing) = false;

        % Create updated simulation state that will describe system after a
        % condensed species has disappeared:
        simulation_state_new.dynamic_params.n_species = n_species;
        simulation_state_new.dynamic_logicals.is_stable = is_stable;

        % WARNING: Potential for bugs in above/below few lines, if values
        % are not properly overwritten.  Also, note, not updating system
        % completely — e.g., surface area is not updated.

        % Calculate instantaneous fluxes for system after a condensed
        % species has disappeared, by calling function:
        [ ~, ~, ~, ~, ~, productionRates ] = ...
            calculate_vapor_pressures_and_instantaneous_fluxes_20240915(...
            simulation_state_new);

        % ASIDE: Dimensions of flux-related vectors
        % rxnRates(:,i) = % (# of rxns) x (1)
        % evaporationRates(:,i) = % (# of species) x (1)
        % consumptionRates(:,i) = % (# of species) x (1)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Calculate time for condensed species to disappear:

        consumptionRates = - productionRates;

        isDecreasing = consumptionRates > 0;
        is_decreasingSpecies = is_cond & is_stable & isDecreasing;
        
        timeToDisappear_vector = nan_vector;
        timeToDisappear_vector(is_decreasingSpecies & ~is_traceSpecies) = ... 
            n_species(is_decreasingSpecies & ~is_traceSpecies) ./ ...
            consumptionRates(is_decreasingSpecies & ~is_traceSpecies);

        % Calculate time until next species/phase disappearance:

        [ timeToDisappearance, index_speciesDisappearing ] = ... 
            min( timeToDisappear_vector );


        % Update vectors that save number of intervals, duration of
        % intervals, and fluxes within intervals:

        saved_intervalCount(i,1) = i;
        saved_time(i,1) = dt - sum(saved_time);
        saved_productionRates(:,i) = productionRates;

    end

    % Multiple productionFlux Matrix by time vector to obtain time-weighted
    % average of productionFlux for all species:

    productionFluxes_timeAveraged = ...
        ( saved_productionRates * saved_time ) ./ dt ;


    % ASIDE: stoichiometricMatrix * rxnsRates = productionFluxes ...
    % but different stoichiometric matrix for each time interval
    % ASIDE: timeWeighted_rxnRates = ( saved_rxnRates * saved_time ) ./ dt ;

end


%% CALCULATE TIME-AVERAGED FLUXES FOR TRACE SPECIES

% Prepare
traceSubstanceIdentifier = substanceIdentifier(is_traceSpecies & is_cond); 
uniqueTraceSubstanceIdentifier = unique(traceSubstanceIdentifier);
numUniqueTraceSubstances = length(uniqueTraceSubstanceIdentifier);

% For each trace substance: e.g., Nd2O3 in multiple condensed phases
for i = 1:numUniqueTraceSubstances

    idx_thisSubstance = uniqueTraceSubstanceIdentifier(i);
    is_thisSubstance = substanceIdentifier == idx_thisSubstance; 

    % calculate total amount of trace species in condensed system,
    % accounting for its presence within one or more phases:
    n_substanceCss = sum(n_species(is_thisSubstance)); 

    % calculate time-averaged flux of this substance 
    flux_substance = sum(productionFluxes_timeAveraged(is_thisSubstance));

    % calculate if this substance is exhausted - rather than before,
    % substance is considered exhausted if greater than or equal to 50% 
    % (rather than 100%) is removed in time step
    is_thisSubstanceExhausted = ...
        (-flux_substance)*dt > (0.5 * n_substanceCss);

    % HARDCODED 0.5 — Can make this line much more 
    % sophisticated/precise by calculating how much material would remain 
    % if residue was at equilibrium with ambient gas at end of time step  

    % If species in exhausted, calculate flux from 0.5 amount of species in
    % condensed system divided by time step — I.e., flux is defined such
    % that half of species will evaporate in single time step
    if is_thisSubstanceExhausted  
        flux_substanceNew = 0.5 * (-n_substanceCss/dt);
        scalingMultiplier = flux_substanceNew/flux_substance;
        productionFluxes_timeAveraged(is_thisSubstance) = ...
            scalingMultiplier * ...
            productionFluxes_timeAveraged(is_thisSubstance);
    end

end

 
end %% End of function






