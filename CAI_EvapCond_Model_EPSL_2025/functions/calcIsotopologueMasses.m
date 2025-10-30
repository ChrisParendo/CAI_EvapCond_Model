function [ massIsotopologues ] = calcIsotopologueMasses( ...
    formulaMatrix, isElementIsoTracked, T_isotopes, massElements, isGas, elements)

%CalcMassesIsotopologues calculates the molar masses of isotopologues.
%
% Inputs required:
% (1) formula matrix
% (2) logical indicating which elements are isotopically tracked
% (3) table/array with isotopic masses
% (4) molar masses of all elements (i.e., masses from periodic table)
%
% Output:
% (1) an array of dimensions (# of species) x (# of isotopes = 2) that
%    contains the masses of isotopologues (e.g.: the mass of 40Ca and 40CaO
%    in column 1, and the masses of 44Ca and 44CaO in column 2) 


idxElementsIsoTracked = find(isElementIsoTracked);
numElementsIsoTracked = sum(isElementIsoTracked);

% Preallocate:
numSpecies = size(formulaMatrix,2);
numElements = size(formulaMatrix,1);

massIsotopologues = nan(numSpecies,2);

% for each element (e.g.: Ca) that is isotopically tracked
for i = 1:numElementsIsoTracked

    % Identify (gas) species that contain this element

    % index identifying this element
    idxThisElement = idxElementsIsoTracked(i);

    % logical identifying this element
    isThisElement = false(numElements,1);
    isThisElement(idxThisElement) = true;    

    % logical identifying gas species that contain this element
    isGasContainingThisElement = isGas & ...
        transpose( any( formulaMatrix & isThisElement) );   

    % Extract 2-element vector with masses of isotopes of element

    selectRows = ismember( T_isotopes.element, elements(idxThisElement) );

    isoMassesThisElement = T_isotopes.mass(selectRows);

    % for each of two isotopes (e.g.: 40Ca & 44Ca) of this element
    for j = 1:2

        % take mass of one isotope
        massThisIsotope = isoMassesThisElement(j);

        % mass of elements that make up species, using specific isotopes
        % where needed
        massConstituentElements = massElements;
        massConstituentElements(idxThisElement) = massThisIsotope;

        % formula matrix for (gas) species that contain this element
        formulaMatrixTheseSpecies = ...
            formulaMatrix(:,isGasContainingThisElement);

        % masses of (gas) species of this element that contain this isotope
        massTheseIsotopologues = transpose( ...
            sum(formulaMatrixTheseSpecies .* massConstituentElements, 1) );


        % Save for output (2D array): save masses for species of this 
        % element that contain this isotope

        massIsotopologues(isGasContainingThisElement,j) = ...
            massTheseIsotopologues;


    end

end


end





    % thisElementName = element(idxThisElement);
