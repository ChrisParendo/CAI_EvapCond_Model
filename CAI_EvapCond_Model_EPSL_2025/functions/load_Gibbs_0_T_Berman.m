function [species_Berman,Gibbs_0_T_species_Berman] = load_Gibbs_0_T_Berman(T)
% load_Gibbs_0_T_Berman — Loads data and computes Gibbs_0_T for species at
% temperature(s) prescribed by the input temperature vector, T
%
% outputs:
% (1) a column vector that lists species
% (2) an array of values of standard Gibbs energy - rows are species,
% columns are temperatures


%%% Ensure that input vector is treated as row vector:
if iscolumn(T)
    T = transpose(T);
end


%%% Load data and create table

file_name = fullfile('data', 'My_Berman_Thermo_SpreadSheet.xlsx');
sheet = 'Berman_All';
range = 'B3:K31';

T_ThermoChemData_Berman = readtable(file_name,'sheet',sheet,'range',range);
T_ThermoChemData_Berman.Properties.RowNames = T_ThermoChemData_Berman.Species;


%%% Create vector of species

species_Berman = T_ThermoChemData_Berman.Species;


%%% Prepare Heat Capacity Eqn. terms used in calculating Enthalpy + Entropy

Heat_Capacity_Eqn_Coefficients = ...
    T_ThermoChemData_Berman{:,{'A','C','D','F'}};

T298 = 298.15;

T_terms_Enthalpy = ...
    [ T-T298; -(T.^-1-T298.^-1); 2*(T.^(1/2)-T298.^(1/2)); log(T/T298) ];

T_terms_Entropy = [ log(T/T298); -1/2*(T.^-2-T298.^-2); ...
    -2*(T.^(-1/2)-T298.^(-1/2)); -(T.^-1-T298.^-1)];


%%% Compute Enthalpy of species as function of temperature

% Enthalpy of formation at 298 K (in J/mol)
Enthalpy_f_0_298 = 1e3 * T_ThermoChemData_Berman.Enthalpy_f_0_298;

% Enthalpy of species — standard elemental reference (SER) model — (in J/mol):
Enthalpy_0_T = ...
    Enthalpy_f_0_298 + Heat_Capacity_Eqn_Coefficients * T_terms_Enthalpy;


%%% Compute Entropy of species as function of temperature

% Entropy at 298 K (in J*K/mol)
Entropy_0_298 = T_ThermoChemData_Berman.Entropy_0_298;

% Entropy of species (J*K/mol)
Entropy_0_T = ...
    Entropy_0_298 + Heat_Capacity_Eqn_Coefficients * T_terms_Entropy;


%%% Calculate standard Gibbs energy of species (H-SER model) as function of T

Gibbs_0_T_species_Berman = Enthalpy_0_T - T .* Entropy_0_T;

end

