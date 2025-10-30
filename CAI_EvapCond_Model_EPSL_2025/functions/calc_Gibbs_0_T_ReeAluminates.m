function [ species_REEaluminates, Gibbs_0_T_species_REEaluminates ] = ...
    calc_Gibbs_0_T_ReeAluminates( T, ...
    species_Berman, Gibbs_0_T_species_Berman, ...
    species_REEsolids, Gibbs_0_T_species_REEsolids )

%   calc_Gibbs_0_T_ReeAluminates calculates Gibbs free energy of REE
%   aluminates
%
%   This function loads data from Table 2 of Wu & Pelton (1992)
%
%   Unlike other functions that I am using that load thermodynamic data and
%   calculate Gibbs energies for a range of temperatures, this one requires
%   more than simply one thermodynamic source.  That is, Wu & Pelton
%   provide enthalpy and entropy of formation from oxides — and hence, for
%   me to calculate useful values requires that I also load my preferred
%   thermodynamic data for these constituent oxides.  For example, I will
%   use data for La2O3 and Al2O3 from other sources, in addition to data
%   from Wu and Pelton, in order to obtain Gibbs energies of LaAlO3.

% Ensure that input vector is treated as row vector:
if iscolumn(T)
    T = transpose(T);
end


%%% Load data and create table

fileName = [ '/Users/chrisparendo/Documents/Hub_Literature_Data/' ...
    'WuPelton1992/WuPelton1992.xlsx'];
sheet = 'Table2';
range = 'B4:J55';

T_RAl = readtable(fileName,'sheet',sheet,'range',range);

numCompounds = size(T_RAl,1);

% Extract data for Al2O3(s) based on Berman dataset
string_Al2O3_s = 'Al2O3_s';
is_Al2O3_s = ismember(species_Berman,string_Al2O3_s);
Gibbs_0_T_Al2O3_s = Gibbs_0_T_species_Berman(is_Al2O3_s,:);

% Preallocate
species_REEaluminates = cell(numCompounds,1);
Gibbs_0_T_species_REEaluminates = nan(numCompounds,length(T));

for i = 1:numCompounds

    % Part 1: Extract Gibbs Energies of Appropriate Oxides
    
    % Identify appropriate R2O3 and extract Gibbs energy for specified T's

    string_thisR2O3_s = [ T_RAl.R2O3{1} '_s' ];
    is_thisR2O3_s = ismember(species_REEsolids,string_thisR2O3_s);
    Gibbs_0_T_R2O3_s = Gibbs_0_T_species_REEsolids(is_thisR2O3_s,:);

    % Stoichiometric coefficients for reactants to make compound

    v1 = 1/2 * T_RAl.x(i); % stoichiometric coefficient for R2O3 reactant
    v2 = 1/2* T_RAl.y(i);  % stoichiometric coefficient for Al2O3 reactant

    % Sum of gibbs energies of reactants

    Gibbs_0_T_Reactants = v1 .* Gibbs_0_T_R2O3_s + v2 .* Gibbs_0_T_Al2O3_s;

    % Part 2: Calc. Gibbs Energy of Formation from Oxides

    H_0_formationfromOxides_calories = T_RAl.H_0_formationOxides(i);
    S_0_formationfromOxides_calories = T_RAl.S_0_formationOxides(i);

    Gibbs_0_formationfromOxides_calories = ...
        H_0_formationfromOxides_calories ...
        - T * S_0_formationfromOxides_calories;

    Gibbs_0_formationfromOxides_joules = ...
        4.184 * Gibbs_0_formationfromOxides_calories;

    % Part 3: Sum

    Gibbs_0_T_species_thisREEaluminate = ...
        Gibbs_0_T_Reactants + Gibbs_0_formationfromOxides_joules;

    % Save:
    
    species_REEaluminates(i) = T_RAl.Compound(i);

    Gibbs_0_T_species_REEaluminates(i,:) = ...
        Gibbs_0_T_species_thisREEaluminate;

end




end