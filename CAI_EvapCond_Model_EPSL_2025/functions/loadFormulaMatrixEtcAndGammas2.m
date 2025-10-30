function [ formulaMatrix, elements, species, phases, substances, ...
    gammas, isTraceSpecies, stickingCoeffs ] = loadFormulaMatrixEtcAndGammas2( ...
    fileName, sheet, range)


% loadFormulaMatrixEtcAndGammas
%
% loads info from Excel Spreadsheet:

% Load cell array that contains formula matrix and vectors w/ species,
% elements, phases, substances, gammas (activity coefficients, and the
% logical 'isTraceSpecies' indicating which species are trace constitutents
% in solid solution within minerals (e.g., Lu2O3 in hibonite is a trace 
% species in the phase hibonite/CaAl12O19).  

T = readcell(fileName,'Sheet',sheet,'Range',range,'TextType','string');

T_column1 = reshape( [T{:,1}], size(T(:,1)) );

% all elements from periodic table
allElements = { 'H'; 'He'; 'Li'; 'Be'; 'B'; 'C'; 'N'; 'O'; 'F'; 'Ne'; ...
    'Na'; 'Mg'; 'Al'; 'Si'; 'P'; 'S'; 'Cl'; 'Ar'; 'K'; 'Ca'; 'Sc'; ...
    'Ti'; 'V'; 'Cr'; 'Mn'; 'Fe'; 'Co'; 'Ni'; 'Cu'; 'Zn'; 'Ga'; 'Ge'; ...
    'As'; 'Se'; 'Br'; 'Kr'; 'Rb'; 'Sr'; 'Y'; 'Zr'; 'Nb'; 'Mo'; 'Tc'; ...
    'Ru'; 'Rh'; 'Pd'; 'Ag'; 'Cd'; 'In'; 'Sn'; 'Sb'; 'Te'; 'I'; 'Xe'; ...
    'Cs'; 'Ba'; 'La'; 'Ce'; 'Pr'; 'Nd'; 'Pm'; 'Sm'; 'Eu'; 'Gd'; 'Tb'; ...
    'Dy'; 'Ho'; 'Er'; 'Tm'; 'Yb'; 'Lu'; 'Hf'; 'Ta'; 'W'; 'Re'; 'Os'; ...
    'Ir'; 'Pt'; 'Au'; 'Hg'; 'Tl'; 'Pb'; 'Bi'; 'Po'; 'At'; 'Rn'; 'Fr'; ...
    'Ra'; 'Ac'; 'Th'; 'Pa'; 'U'; 'Np'; 'Pu'; 'Am'; 'Cm'; 'Bk'; 'Cf'; ...
    'Es'; 'Fm'; 'Md'; 'No'; 'Lr'; 'Rf'; 'Db'; 'Sg'; 'Bh'; 'Hs'; 'Mt'; ...
    'Ds'; 'Rg'; 'Cn'; 'Nh'; 'Fl'; 'Mc'; 'Lv'; 'Ts'; 'Og' };

% formula matrix
isElementsRow = ismember(T_column1,allElements);

formulaMatrixCell = T(isElementsRow,2:end);
formulaMatrix = ...
    reshape( [formulaMatrixCell{:}], size(formulaMatrixCell) );

% elements vector
elementsCell = T(isElementsRow,1);
elements = reshape( [elementsCell{:}], size(elementsCell(:)) );

% logicals for later use
isSpeciesRow = ismember(T_column1,'species');
isPhasesRow = ismember(T_column1,'phases');
isSubstancesRow = ismember(T_column1,'substances');
isGammasRow = ismember(T_column1,'gammas');
isIsTraceSpeciesRow = ismember(T_column1,'isTraceSpecies');

%
isStickingCoeffsRow = ismember(T_column1,'stickingCoeffs'); %%%%%% Modified August 12, 2024 (version 2)
%
stickingCoeffsCell = T(isStickingCoeffsRow,2:end)   %%%%%% Modified August 12, 2024 (version 2)
stickingCoeffs = reshape( [stickingCoeffsCell{:}], size(stickingCoeffsCell(:)) );   %%%%%% Modified August 12, 2024 (version 2)

% species vector
speciesCell = T(isSpeciesRow,2:end);
species = reshape( [speciesCell{:}], size(speciesCell(:)) );

% phases vector
phasesCell = T(isPhasesRow,2:end);
phases = reshape( [phasesCell{:}], size(phasesCell(:)) );

% substances vector
substancesCell = T(isSubstancesRow,2:end);
substances = reshape( [substancesCell{:}], size(substancesCell(:)) );

% gammas vector
gammasCell = T(isGammasRow,2:end);
gammas = reshape( [gammasCell{:}], size(gammasCell(:)) );

% 'is trace species' logical/vector
isTraceSpeciesCell = T(isIsTraceSpeciesRow,2:end);
isTraceSpecies = logical( ...
    reshape( [isTraceSpeciesCell{:}], size(isTraceSpeciesCell(:)) ) );


end