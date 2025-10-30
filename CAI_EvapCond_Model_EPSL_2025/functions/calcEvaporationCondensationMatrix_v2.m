function [EvaporationStoichiometricMatrix] = ...
    calcEvaporationCondensationMatrix_v2( formula_matrix,is_gas, is_cond)

% calcEvaporationCondensationMatrix
% creates a stoichiometric matrix for all evaporation(-condensation) rxns,
% one for each condensed species


num_cond = sum(is_cond);
num_species = size(formula_matrix,2);

numberedCondensedSpecies = (cumsum(is_cond) .* is_cond);

%%% Loop over condensed species to generate evaporation rxns for each:

% preallocate
EvaporationStoichiometricMatrix = nan(num_species,num_cond);

for i = 1:num_cond

    is_particularCondensedSpecies = i == numberedCondensedSpecies;

    is_subsystem = is_gas | is_particularCondensedSpecies;

    A_subsystem = formula_matrix(:,is_subsystem);

    rxns = null(A_subsystem,'rational');            %%% Changed 'rxn' to 'rxns' - Nov 13, 2023

    %%% Addition from Nov 13, 2023 - Function renamed 'v2' %%%
    % Required when modifying the code to accomodate multiple gas species,
    % e.g., Nd and NdO
    % In short, the null space of the formulat matrix now will yield more
    % than 1 rxn ... this modification serves to select 1 rxn when more
    % than 1 occur: The chosen rxn must include the condensed species, and
    % only one rxn per condensed species is needed/wanted.  

    % Logical to identify the condensed species in the subsystem
    isCondensedSpeciesAmongSubsystemSpecies = ...
        is_particularCondensedSpecies(is_subsystem);

    % Identify which reaction(s) include the condensed species
    isRxnIncludingCondensedSpecies = ...
        any( rxns & isCondensedSpeciesAmongSubsystemSpecies);

    % In case more than one rxn can be written for condensed species,
    % choose just one (the first)

    % index
    indexRxn = find( isRxnIncludingCondensedSpecies, 1 );

    % reaction
    chosenRxn = rxns(:,indexRxn);

    %%% End of addition from Nov 13, 2023 %%%

    rxn_expanded = zeros(num_species,1);
    rxn_expanded(is_subsystem) = chosenRxn;           %%% Edited - Nov 13, 2023

    % Use sign convention where coefficients of condensed species are '-':
    if rxn_expanded(is_particularCondensedSpecies) > 0 
        rxn_expanded = -(rxn_expanded);
    end

    EvaporationStoichiometricMatrix(:,i) = rxn_expanded;

end

end
