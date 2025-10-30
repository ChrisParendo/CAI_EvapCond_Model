function [ species_TiMinerals, Gibbs_0_T_species_TiMinerals ] = ...
    load_Gibbs_0_T_TiMinerals(T)
%load_Gibbs_0_T_TiPhases Loads data and computes Gibbs_0_T for Ti-bearing
%species (Perovskite, possibly others) at temperatures prescribed by the
%input vector, T


% Ensure that input vector is treated as row vector:
if iscolumn(T)
    T = transpose(T);
end


% Load data and create table

fileName = fullfile('data', 'My_Perovskite_Thermo_SpreadSheet.xlsx');
sheet = 'Perovskite';
range = 'B3:K5';

T_TiMinerals = readtable(fileName, 'sheet', sheet, 'range', range);


% Create vector of species

species_TiMinerals = T_TiMinerals.Species;


% Prepare empirical fit terms to be used in calculating H and S

heatCapCoefficients = T_TiMinerals{:,{'k0', 'k1_x_1eMinus2', ...
    'k2_x_1eMinus5', 'k3_x_1eMinus7'}} .* [1, 1e2, 1e5, 1e7];


% Prepare terms to be multiplied by fit parameters

T298 = 298.15;

T_terms_Enthalpy = [ T-T298; 2*(T.^(1/2)-T298.^(1/2)); ...
    -1*(T.^(-1)-T298.^(-1)); -1/2*(T.^(-2)-T298.^(-2))  ]; 

T_terms_Entropy = [ log(T./T298); -2*(T.^(-1/2)-T298.^(-1/2)); ...
    -1/2*(T.^(-2)-T298.^(-2)); -1/3*(T.^(-3)-T298.^(-3)) ];


% Compute Enthalpy (in J/mol) of species as a function of temperature

%%% Compute Enthalpy (in J/mol) of species as function of temperature

Enthalpy_f_0_298_kJ_per_mol = T_TiMinerals.Enthalpy_f_0_298;
Enthalpy_f_0_298_J_per_mol = 1e3 * Enthalpy_f_0_298_kJ_per_mol;

Enthalpy_0_T = ...
    Enthalpy_f_0_298_J_per_mol + heatCapCoefficients * T_terms_Enthalpy;


%%% Compute Entropy (in J/(K*mol) of species as function of temperature

Entropy_0_298_J_per_K_mol = T_TiMinerals.Entropy_0_298;

Entropy_0_T = Entropy_0_298_J_per_K_mol + heatCapCoefficients * T_terms_Entropy;


%%% Calc. standard Gibbs energy of species (H-SER model) as function of T

Gibbs_0_T_species_TiMinerals = Enthalpy_0_T - T .* Entropy_0_T;


end



