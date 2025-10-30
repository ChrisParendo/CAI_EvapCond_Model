function [ db_dt_isotopes, x_isotopesCss, x_isotopesGss ] = ...
    calculateEvaporativeIsotopeFluxes20240922( ...
    isElementIsoTracked, formulaMatrix, n_species, n_saturation, ...
    isGas, isotopologueMasses, molarMassSpecies_g_per_mol,     ...
    dn_dt_gss, b_isotopes_gss, b_isotopes_css, b_elements_gss, b_elements_css)

% This function is adapted from calcECFluxIso_v3
%   
% This function calculates evaporative fluxes of isotopologues.  
%
% Inputs:
% (1) isElementIsoTracked
% (2) formulaMatrix
% (3) n_species
% (4) n_saturation
% (5) isGas
% (6) isotopologueMasses
% (7) molarMassSpecies_g_per_mol
% (8) dn_dt_gss
% (9) b_isotopes_gss
% (10) b_isotopes_css
% (11) b_elements_gss
% (12) b_elements_css
%
% Outputs:
% (1) db_dt_isotopes
% (2) x_isotopesCss
% (3) x_isotopesGss


%%% Preparation

numElementsIsoTracked = sum(isElementIsoTracked);
idx_elementsIsoTracked = find(isElementIsoTracked);

numElements = size(formulaMatrix,1);

db_dt_isotopes = nan(numElements,2); % One way to organize

x_isotopesGss = b_isotopes_gss ./ b_elements_gss; % mole fraction of isotopes of each element
x_isotopesCss = b_isotopes_css ./ b_elements_css; % mole fraction of isotopes of each element

%%% Added Nov. 29
isNotInCss = b_elements_css ./ (b_elements_gss + b_elements_css) < 1e-3;                     %%% Temporarily hardcoded
x_isotopesCss(isNotInCss,:) = x_isotopesGss(isNotInCss,:);
%%% Added Nov. 29


%%% Calculation of evaporative fluxes

% For each element that is isotopically tracked, calculate evaporative 
% fluxes of gas isotopologues (e.g.: 40Ca, 44Ca, 40Ca16O, & 44Ca16O), 
% using previous determination of evaporative fluxes of 'generic' gas
% species (e.g., Ca & CaO)

% For each isotopically tracked element
for i = 1:numElementsIsoTracked 

    % Identify gas species (e.g.: Ca & CaO) that contain this element
    % (e.g., Ca)

    % index identifying this element:
    idx_thisElement = idx_elementsIsoTracked(i); 

    % logical identifying this element
    isThisElement = false(numElements,1);
    isThisElement(idx_thisElement) = true;    

    % logical identifying gas species that contain this element
    isGasContainingThisElement = isGas & ...
        transpose( any( formulaMatrix & isThisElement) );    

    % Calculate abundance of each isotopologue (e.g.: 40Ca, 44Ca, ...
    % 40Ca16O, & 44Ca16O) in ambient gas - Output is in form of 2D array, 
    % with species as rows (e.g.: Ca & CaO) and isotopes as columns 
    % (e.g.: 40 & 44)

    n_isotopologuesGss = ...
        n_species(isGasContainingThisElement) * x_isotopesGss(isThisElement,:) ;

    % Calculate abundance of each isotopologue (e.g.: 40Ca, 44Ca, ...
    % 40Ca16O, & 44Ca16O) in gas in equilibrium with particle - Output is,
    % as above, in form of 2D array

    n_isotopologuesSat = ...
        n_saturation(isGasContainingThisElement) * x_isotopesCss(isThisElement,:) ;

    % Calculate evaporative flux of each isotopologue of this element - 
    % Output is, as above, in form of 2D array

    mass_theseIsotopologues = isotopologueMasses(isGasContainingThisElement,:);
    mass_referenceSpecies = molarMassSpecies_g_per_mol(isGasContainingThisElement);

    dn_dt_species = dn_dt_gss(isGasContainingThisElement);
    n_speciesSat = n_saturation(isGasContainingThisElement);
    n_speciesGss = n_species(isGasContainingThisElement);

    dn_dt_isotopologues =  dn_dt_species .* ...
      (n_isotopologuesSat - n_isotopologuesGss) ./ sqrt(mass_theseIsotopologues) ...
      ./ ( (n_speciesSat - n_speciesGss) ./ sqrt(mass_referenceSpecies) );
    is_NotDefined = isnan(dn_dt_isotopologues) | isinf(dn_dt_isotopologues) ;  % to account for zeros in denominator in above
    dn_dt_isotopologues(is_NotDefined) = 0;

    % Evaporative isotopic flux for particular element:
    db_dt_isotopesThisElement = formulaMatrix(isThisElement,isGasContainingThisElement) * dn_dt_isotopologues;

    % Save the isotopic flux in 2D array with values for all isotopes
    db_dt_isotopes( idx_thisElement,1:2 ) = db_dt_isotopesThisElement;   % One way to organize


end


end