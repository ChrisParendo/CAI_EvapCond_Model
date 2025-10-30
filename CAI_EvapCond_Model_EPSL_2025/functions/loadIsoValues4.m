function [ T_isotope,  b_isotopes, b_isotopes_css, b_isotopes_gss ] ...
    = loadIsoValues4( ...                                                     % version 2 adds more isotope systems than Ca (44,40) - start with Nd
    elements, b_elements, b_elements_css, b_elements_gss)                   % version 4 adds more isotopes sytems, changing range and sheet for spreadsheet

% loadIsoValues3 Summary of this function goes here
%   Detailed explanation goes here

%%% Load/Create isotope table

% isotope = { '44Ca'; '40Ca'; '146Nd'; '144Nd' };
% element = { 'Ca'; 'Ca' ; 'Nd'; 'Nd' };
% fraction = [ 0.02; 0.97; 0.172; 0.238 ];
% mass = [ 43.9554; 39.9625; 145.9131; 143.9101 ];
% 
% T_isotope = table( isotope, element, fraction, mass );

file = fullfile('data', 'IsotopeInputEvaporationCode.xlsx');
sheet = 'Nov16';
range = 'B2:E18';

T_isotope = readtable(file,'sheet',sheet,'range',range);



isElementIsoTracked = ismember( elements, T_isotope.element );

idxElementsIsoTracked = find(isElementIsoTracked);
numElementsIsoTracked = sum(isElementIsoTracked);

% Preallocate:
b_isotopes = nan( length(elements), 2 );
b_isotopes_css = nan( length(elements), 2);
b_isotopes_gss = nan( length(elements), 2);

% for each element (e.g.: Ca) that is isotopically tracked
for i = 1:numElementsIsoTracked

    % index identifying this element
    idxThisElement = idxElementsIsoTracked(i);  

    % Extract 2-element vector 

    selectRows = ismember( T_isotope.element, elements(idxThisElement) );

    b_isotopes(idxThisElement,:) = b_elements(idxThisElement) * (T_isotope.fraction(selectRows))';

    b_isotopes_css(idxThisElement,:) = b_elements_css(idxThisElement) * (T_isotope.fraction(selectRows))';

    b_isotopes_gss(idxThisElement,:) = b_elements_gss(idxThisElement) * (T_isotope.fraction(selectRows))';

end

end
