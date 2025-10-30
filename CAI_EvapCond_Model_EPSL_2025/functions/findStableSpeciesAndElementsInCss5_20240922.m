function [is_speciesStableInCSS,is_elementStableInCSS, ...
    is_speciesInCSSComprisedOfUnstableElementInCSS, ...
    is_speciesPermissibleInCSS] = ...
    findStableSpeciesAndElementsInCss5_20240922( ...
    formulaMatrix, ...          % (1)
    n_species, ...              % (2)
    b_elements, ...             % (3)
    significance_threshold, ... % (4)
    is_css, ...                 % (5)
    phaseIdentifier, ...
    is_speciesInMultiSpeciesPhase, ...
    is_traceSpecies)

%findStableSpeciesAndElementsInCSS 
% (1) Calculates an element distribution matrix - a matrix where entries
% provide the fraction of an element in each species (all entries for a
% single element will sum to one)
% (2) Determines which species host a significant amount of an element —
% e.g., species that host greater than 1e-9 (or other small value) of an
% element as a fraction of the total amount of that element in systems.
% Specifically, identifies these in CSS.
% (3) Determines which elements are hosted in significant amount in CSS
% (specifically which elements are hosted in a significant amount by any
% CSS species)

            % DOES THIS WORK FOR MAJOR (SPECIES-FORMING) ELEMENTS AND TRACE ELEMENTS?
            % Think so — but need to consider if any problem for trace
            % elements. - PROBABLY ONLY WANT THIS FUNCTION TO UTILIZE
            % SPECIES-FORMING ELEMENTS, NOT REEs, ALSO NEED TO CHANGE LINES
            % 50 & 63 - NEED TO ID SUBSTANCES (NOT MERELY SPECIES) IN CODE
            % THAT HANDLES TRACE ELEMENTS

% Outputs are logicals of length # of species: entries for gas species will
% always be 0 (i.e., false)




% Calculate element distribution matrix — i.e., the fraction of an element
% in each of the species
eleDistMatrix = formulaMatrix .* transpose(n_species) ./ b_elements;

% Identify species that host any element in significant amounts (above 
% threshold)
is_speciesHostingAnyElementInSigAmt = ...
    transpose(any(eleDistMatrix > significance_threshold, 1));

% Identify species that are both hosting significant amounts of an element 
% and are present in the condensed system
is_speciesSignificantInCSS = is_speciesHostingAnyElementInSigAmt & is_css;

% Identify elements that are hosted in significant amounts in any condensed
% species
is_elementStableInCSS = ...
    any(eleDistMatrix(:, is_css) > significance_threshold, 2);

% Identify elements that are not hosted in significant amounts in any 
% condensed species
is_elementUnstableInCSS = ~ is_elementStableInCSS;

% Identify species that are comprised of elements that are not present in
% significant amounts in the condensed system
is_speciesInCSSComprisedOfUnstableElementInCSS = ...
    any(is_elementUnstableInCSS & formulaMatrix & is_css', 1)';

% Identify condensed species that are permissible — i.e., are comprised
% only of elements that are present in significant amounts in the condensed
% system (Note: species that are "permissible" may or may not be "stable"
is_speciesInCssPermissible = ...
    ~ is_speciesInCSSComprisedOfUnstableElementInCSS;

% Modify "is species permissible" vectory by considering trace species that
% are in solution in host species: the trace species are not permissible if
% the host phase/species is not permissible:

multiSpeciesPhaseIdentifier = phaseIdentifier(is_speciesInMultiSpeciesPhase);   

uniqueMultiSpeciesPhases = unique(multiSpeciesPhaseIdentifier);

numMultiSpeciesPhases = length(uniqueMultiSpeciesPhases);

% for each phase that consists of multiple species - i.e., solutions
for i = 1:numMultiSpeciesPhases

    idx = uniqueMultiSpeciesPhases(i);

    isPhase = idx == phaseIdentifier;

    isHostSpecies = isPhase & ~ is_traceSpecies;
    isHostedSpecies = isPhase & is_traceSpecies;

    % Determine whether the host species is permissible
    isHostSpeciesPermissibleInCSS = ...
        is_speciesInCssPermissible(isHostSpecies);

    % If the host species is not permissible, then mark the "hosted"
    % species — e.g., trace species in solution in a mineral — as also not
    % permissible
    if all(~isHostSpeciesPermissibleInCSS) %changed any() to all() Sept., 2024
        is_speciesInCssPermissible(isHostedSpecies) = false;
    end

end


%
is_speciesPermissibleInCSS = is_css & is_speciesInCssPermissible;




% Identify stable condensed species
is_speciesStableInCSS = is_speciesSignificantInCSS & is_speciesPermissibleInCSS;

% Identify stable elements in condensed system
is_elementStableInCSS = ...
    any( (logical(formulaMatrix) & is_speciesStableInCSS' ), 2 ) ;


end