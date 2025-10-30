function [species_NIST, Gibbs_0_T_species_NIST] = load_Gibbs_0_T_NIST_v2(T) % version 2 - Nov19, 2023, different 'sheet' called
% load_Gibbs_0_T_NIST - Loads data and computes Gibbs_0_T for species
%   at temperature(s) prescribed by the input temperature vector, T
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

file_name = fullfile('data', 'My_NIST_Thermo_SpreadSheet.xlsx');
sheet = 'NIST_Compiled';
range = 'B3:AL28';

T_ThermoChemData_NIST = readtable(file_name,'sheet',sheet,'range',range);
T_ThermoChemData_NIST.Properties.RowNames = T_ThermoChemData_NIST.Species;


%%% Create vector of species

species_NIST = T_ThermoChemData_NIST.Species;


%%% Create vectors indicating temperature boundaries

T_lim1 = T_ThermoChemData_NIST.T_Int2_LowLim;
T_lim2 = T_ThermoChemData_NIST.T_Int3_LowLim;


%%% Create logical indices - for species by temperature array, indicating
%   which heat capacity equation/coefficients to use for calculating each
%   value

% entries/values to be computed by equation 1:
i1 = T <= T_lim1 | isnan(T_lim1);
% entries/values to be computed by equation 2:
i2 =  T > T_lim1 & ( T <= T_lim2 | isnan(T_lim2) );
% entries/values to be computed by equation 3:
i3 = T > T_lim2 ;


%%% Prepare Shomate terms to be used in calculating Ethalpy and Entropy

Shomate_Parameters_1 = ...
    T_ThermoChemData_NIST{:,{'A1' 'B1' 'C1' 'D1' 'E1' 'F1' 'G1' 'H1'}};

Shomate_Parameters_2 = ...
    T_ThermoChemData_NIST{:,{'A2' 'B2' 'C2' 'D2' 'E2' 'F2' 'G2' 'H2'}};

Shomate_Parameters_3 = ...
    T_ThermoChemData_NIST{:,{'A3' 'B3' 'C3' 'D3' 'E3' 'F3' 'G3' 'H3'}};

% define 't' in accord with Shomate eqn. of NIST webbook
t = T/1e3;

t_terms_Enthalpy = [t; t.^2/2; t.^3/3; t.^4/4; -1./t; ...
    ones(size(t)); zeros(size(t)); -ones(size(t))];

t_terms_Entropy = [log(t); t; t.^2/2; t.^3/3; -1./(2*t.^2); ...
    zeros(size(t)); ones(size(t)); zeros(size(t))];


%%% Compute Enthalpy of species as function of temperature

% Enthalpy of formation at 298 K (in kJ/mol)
Enthalpy_f_0_298_kJ = T_ThermoChemData_NIST.Enthalpy_f_0_298; 

% Enthalpy of species — standard elemental reference (SER) model — (in J/mol):
Enthalpy_0_T_Eqn1 = 1e3 * ...
    ( Enthalpy_f_0_298_kJ + (Shomate_Parameters_1 * t_terms_Enthalpy) );
Enthalpy_0_T_Eqn2 = 1e3 * ...
    ( Enthalpy_f_0_298_kJ + (Shomate_Parameters_2 * t_terms_Enthalpy) );
Enthalpy_0_T_Eqn3 = 1e3 * ...
    ( Enthalpy_f_0_298_kJ + (Shomate_Parameters_3 * t_terms_Enthalpy) );
% 
Enthalpy_0_T = nan(length(species_NIST),length(T));
Enthalpy_0_T(i1) = Enthalpy_0_T_Eqn1(i1);
Enthalpy_0_T(i2) = Enthalpy_0_T_Eqn2(i2);
Enthalpy_0_T(i3) = Enthalpy_0_T_Eqn3(i3);


%%% Compute Entropy of species as function of temperature

% Entropy of species - calculated using equation for each temp interval
Entropy_0_T_Eqn1 =  Shomate_Parameters_1 * t_terms_Entropy ;
Entropy_0_T_Eqn2 =  Shomate_Parameters_2 * t_terms_Entropy ;
Entropy_0_T_Eqn3 =  Shomate_Parameters_3 * t_terms_Entropy ;
% 
Entropy_0_T = nan(length(species_NIST),length(T));
Entropy_0_T(i1) = Entropy_0_T_Eqn1(i1);
Entropy_0_T(i2) = Entropy_0_T_Eqn2(i2);
Entropy_0_T(i3) = Entropy_0_T_Eqn3(i3);


%%% Calculate standard Gibbs energy of species (H-SER model) as function of T

Gibbs_0_T_species_NIST = Enthalpy_0_T - T .* Entropy_0_T;

end

