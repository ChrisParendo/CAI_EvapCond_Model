 %%% CAI EVAPORATION-CONDENSATION (DIS)EQUILIBRIUM MODEL 

clear; close all; clc 

addpath('functions');

% Recommend using warning("off") to suppress MATLAB warning about matrix
% having high condition number; message is expected, and repeated warnings 
% will slightly slow down calculation
warning("off") 


%% INPUTS: INITIAL CONDITIONS AND PARAMETERS

oxygen_multiplier = 0.5; % Factor by which oxygen abundance is mutiplied
% Set to 0.5 if carbon (C) is not explicitly included in the system,
% as approximately 0.5 of oxygen will bind with C to form CO in gas of
% solar composition [Generally, this should be = 0.5]

dust_to_gas_ratio = 1; % Factor by which abundances of elements other 
% than hydrogen (H) will be multiplied; set to 1 for solar composition
% [Generally, this should be = 1]

% Prescribe initial temperature and pressure:

P = 1e-4;   % Pressure (bar)
T_initial = 1445; 

% Other Possible Temperatures:
% T_0 = 1642.34; % (corresponds to 0.01 Ca condensed for P = 1e-4 bar)
% T_0 = 1445;    % (corresponds to 0.01 Ca evaporated for P = 1e-4 bar)

R = 8.314;  % Gas Constant (J/mol*K or m^3*Pa/mol*K)

% % Linear heating Rate
heating_rate_k_per_hour = 100; % Heating rate (K/hour)
heating_rate_k_per_second = heating_rate_k_per_hour / 3600; % (K/s)

% Duration of simulation
t_max = 2e4 ;               % Total simulation time (seconds)
dt = 1e1 ;                  % Time step (seconds)

% % Linear temperature change
t = 0:dt:t_max;                                 % time vector (s)
T = T_initial + heating_rate_k_per_second * t;  % Temperature vector (K)

% % Stepwise temperature change
% t = 0:dt:t_max;
% T = nan(size(t));
% T(1) = T_initial;
% T(2:end) = 1773 + 20;

% Prescribe particle size and density
diameter_particle_mm = 0.1;          % Particle diameter (mm)

% Density assumed to be the same for all condensed phases
density_particle_g_cm3 = 3;          % Particle density (g/cm^3)

num_time_steps = length(t); % number of time steps


%% NUMERICAL PARAMETERS

% Prescribe sequence of values for tau, used in interior point algorithm
tau_0_exp = 5; % exponent for defining inital tau
tau_f_exp = -35; % exopnenent for defining final tau
tau_sequence = ...
    logspace(tau_0_exp, tau_f_exp, 1 * (tau_0_exp - tau_f_exp + 1) );

% Prescribe threshold for significance
significance_threshold = 1e-12;

% Tolerances
tolerance_energyRxn = 1e-12;
tolerance_relativeFlux = 1e-2; % not necessarily used (only in select versions of code)


%% LOAD FORMULA MATRIX, SPECIES INFORMATION

% Load/Prescribe from Spreadsheet:
% (1) formula matrix
% (2) elements
% (3) species
% (4) phases
% (5) substances
% (6) gammas
% (7) logical identifying trace species
% (8) sticking coefficients

% Specify file/data location
file_name = fullfile('data', 'FormulaMatricesTrim.xlsx');   %%% path updated
% sheet = 'Active2024x08x12_stickingCoeffs';
sheet = 'Active2025';
range = 'B4:AV22';

% Load formula matrix and species data
[ formula_matrix, elements, species, phases, substances, ...
    gamma_species, is_traceSpecies, sticking_coeffs ] = ...
    loadFormulaMatrixEtcAndGammas2( file_name, sheet, range );


%% CREATE IDENTIFIER VECTORS FOR PHASES AND SUBSTANCES

% Identify unique phases and assign identifiers
unique_phases = unique(phases,'stable');
num_phases = length(unique_phases);
phase_identifier = nan(length(phases),1); % Initialize with NaN
for i = 1:num_phases
    is_phase = ismember( phases, unique_phases(i) );
    phase_identifier(is_phase) = i;
end

% Identify unique phases and assign identifiers
unique_substances = unique(substances,'stable');
num_substances = length(unique_substances);
substance_identifier = nan(length(substances),1); % Initialize with NaN
for i = 1:num_substances
    is_substance = ismember( substances, unique_substances(i) );
    substance_identifier(is_substance) = i;
end


%% LOAD PERIODIC TABLE, CALCULATE MOLAR MASSES OF SPECIES

% Create table with periodic table (tbl = table)
pt_file = fullfile('data', 'PubChemElements_all.csv');  %%% updated path
tbl_periodic_table = readtable(pt_file);                %%% updated path
tbl_periodic_table.Properties.RowNames = tbl_periodic_table.Symbol;

% Load molar masses of elements 
molar_mass_elements_g_mol = tbl_periodic_table.AtomicMass(elements);

% Calculate molar masses of species
molar_mass_species_g_per_mol = ... 
    transpose( sum(formula_matrix .* molar_mass_elements_g_mol, 1) );
molar_mass_species_kg_per_mol = 1e-3 * molar_mass_species_g_per_mol;


%% LOAD SOLAR SYSTEM COMPOSITION

% Load elemental abundances of solar system (Lodders, 2003)
tbl_solar_system_abundances = load_SolarSys_Abundances_Lodders2003();

% Create vector with mole fractions of select elements in solar system
element_mole_fractions_solar_system = nan(size(elements));
for i = 1:length(elements)
    index = ismember(tbl_solar_system_abundances.Element,elements(i));
    element_mole_fractions_solar_system(i) = ...
        tbl_solar_system_abundances.MoleFraction(index);
end


%% LOAD THERMODYNAMIC DATA AND ACTIVITY MODEL PARAMETERS

% Load & calculate Gibbs Energies of species at all temperatures
mu_0_species = loadGibbs0Tspecies_v4_Dec03(species,T);  

% Load activity model parameters for SACM melt (Berman, 1983) 
[W_sacm, Q_m_sacm] = calc_Berman1983_W_parameters_and_Qm(T);


%% CREATE LOGICAL VECTORS AND INDEXING VECTORS

%%% CREATE INDEXING VECTORS

is_gas = endsWith(species,'_g');
i_gas = find(is_gas);

is_sacm = ismember(species, {'SiO2_l','Al2O3_l','CaO_l','MgO_l'});
i_sacm = [ find(ismember(species,'SiO2_l')); ...
    find(ismember(species,'Al2O3_l')); ...
    find(ismember(species,'CaO_l')); ...
    find(ismember(species,'MgO_l'))]; % i_sacm is ordered Si,Al,Ca,Mg

is_sol = endsWith(species,'_s') | contains(species,'_s');
i_sol = find(is_sol);

is_cond = ~ endsWith(species,'_g') ;
i_cond = find(is_cond);

[~, idx] = unique(phase_identifier) ;
is_speciesInUniSpeciesPhase = ( ismember(1:length(phase_identifier),idx) )';
is_speciesInMultiSpeciesPhase = ~is_speciesInUniSpeciesPhase;


isAl = find(elements == 'Al');
isCa = find(elements == 'Ca');
isNd = find(elements == 'Nd');
isLu = find(elements == 'Lu');

is_K2O_melt = ismember(species, 'K2O_l');
is_H2_g = ismember(species, 'H2_g');
is_H2O_g = ismember(species, 'H2O_g');


%% EQUILIBRATE INITIAL FULL SYSTEM: PARTICLE AND GAS

% Specify relative elemental abundances in vector ("b" is used as symbol
% for element abundance vectors)
b_element_proportions = element_mole_fractions_solar_system;
% Adjust oxygen abundance if necessary
b_element_proportions(ismember(elements,'O')) = ...
    oxygen_multiplier * element_mole_fractions_solar_system(ismember(elements,'O'));

% identify elements in dust: defined here as all elements but hydrogen (H)
is_dust_element = ~ismember(elements,'H');
% Adjust element abudnances if necessary 
b_element_proportions(is_dust_element) = ...
    b_element_proportions(is_dust_element) * dust_to_gas_ratio;
                                

% Find a feasible point to use as initial guess for gibbs minimization
[ n_species_0, rel_violation ] = find_feasible_n_species_v5( ...
    formula_matrix, b_element_proportions);

% Run gibbs minimization algorithm:

% This is all true if all elements Si, Al, Ca, Mg are in system:
is_sacmSpeciesIncluded = true(4,1); 

%%% Prepare input parameters to equilibrium calculation (Added Sept. 28, 2024)
equilibrium_params = struct();
equilibrium_params.n_species = n_species_0;  % n_species is a dynamic variable
equilibrium_params.gamma_species = gamma_species;  % gammas are a dynamic variable
equilibrium_params.formula_matrix = formula_matrix;  
equilibrium_params.b_elements = b_element_proportions;
equilibrium_params.tau_sequence = tau_sequence;
equilibrium_params.thermodynamics = struct(...
    'R', R, ...
    'T', T(1), ...
    'P', P, ...
    'mu_0_species', mu_0_species(:,1));  
equilibrium_params.cmas_melt = struct(...
    'W_sacm', W_sacm(:,1), ...
    'Q_m_sacm', Q_m_sacm);
equilibrium_params.indices = struct(...
    'i_gas', i_gas, ...
    'i_sacm', i_sacm, ...
    'i_sol', i_sol);
equilibrium_params.logicals = struct(...
    'is_sacmSpeciesIncluded', is_sacmSpeciesIncluded, ...
    'is_traceSpecies', is_traceSpecies, ...
    'is_K2O_melt', is_K2O_melt);
equilibrium_params.identifiers = struct(...
    'phase_identifier', phase_identifier);  


[n_species, gammas_sacm, mu_species, log_act_species, zeta_vector] = ...
        minimizeGibbsByInteriorPoint20240928(equilibrium_params);


%% SCALE SYSTEM

% calcualte b_elements for condensed and gas subsystems
b_css_proportions = formula_matrix(:,is_cond) * n_species(is_cond);
b_gss_proportions = formula_matrix(:,is_gas) * n_species(is_gas);

[b_elements, scaling_factor, mass_particle_g] = scale_system_v3( ...
    diameter_particle_mm, density_particle_g_cm3, b_element_proportions, ...
    b_css_proportions, molar_mass_elements_g_mol);

b_css_elements = b_css_proportions * scaling_factor;
b_gss_elements = b_gss_proportions * scaling_factor;

% Change
n_species_scaled = n_species * scaling_factor;


%% DISPLAY SELECT DETAILS OF INITAL SYSTEM

% inital species abunduances
[ species n_species_scaled ]

% initial fraction elements condensed
[ elements b_css_elements./b_elements ] 


%% LOAD ISOTOPE COMPOSITIONS AND ASSIGN TO SUBSYSTEMS

[ T_isotopeNEW, b_isotopesNEW, b_isotopes_css_initialNEW, b_isotopes_gss_initialNEW ] ...
    = loadIsoValues4( elements, b_elements, b_css_elements, b_gss_elements);

isElementIsoTrackedNEW = ismember( elements, T_isotopeNEW.element );

isotopologueMassesNEW = calcIsotopologueMasses( formula_matrix, ...
    isElementIsoTrackedNEW, T_isotopeNEW, molar_mass_elements_g_mol, ...
    is_gas, elements);


%% PREALLOCATE + FILL IN INITIAL VALUE

% Prepare

numSpecies = length(species);
numElements = length(elements);
numSteps = length(T);

numSpecies_css = sum(is_cond);
numSpecies_gss = sum(is_gas);

% Preallocate:

b_css_array = nan(numElements,numSteps);
b_gss_array = nan(numElements,numSteps);

b_isotopes_css = nan(numSteps,numElements,2);
b_isotopes_gss = nan(numSteps,numElements,2);


n_species_css_array = nan(numSpecies_css,numSteps);
n_species_gss_array = nan(numSpecies_gss,numSteps);
n_species_array = nan(numSpecies,numSteps);

massParticle_g_array = nan(numSteps,1);
radiusParticle_m_array = nan(numSteps,1);
areaParticle_m2_array = nan(numSteps,1);

mu_species_css_array = nan(numSpecies_css,numSteps);
mu_species_gss_array = nan(numSpecies_gss,numSteps);
mu_species_array = nan(numSpecies,numSteps);

n_saturation_array = nan(numSpecies,numSteps);

dn_dt_gss_array = nan(numSpecies,numSteps);

gammas_sacm_array = nan(4,numSteps);

mu_sat_vector_array = nan(numSpecies,numSteps);

savedProductionFluxes = nan(numSpecies,numSteps);
savedEvaporationFluxes = nan(numSpecies,numSteps);

volume_gas_cL_array = nan(numSteps,1);
P_envelope_array = nan(numSteps,1);

%% ASSIGN INITAL VALUES TO VECTORS

n_species = n_species_scaled;
dn = zeros(size(n_species)); 

b_isotopes_cssNEW = b_isotopes_css_initialNEW;
b_isotopes_gssNEW = b_isotopes_gss_initialNEW;
db_isotopes_cssNEW = zeros(length(elements),2);
db_isotopes_gssNEW = zeros(length(elements),2);


%% %%%%%%%%%  START FORWARD FINITE DIFFERENCE (ITERATIVE) LOOP  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic % start timer

for i = 1:num_time_steps


%% DISPLAY

% Progress of simulation
fprintf('Progress: %.2f%%\n', 100 * i / num_time_steps);
% { i/num_time_steps }


%% UPDATE VALUES: SPECIES, ELEMENT, ISOTOPE ABUNDANCE VECTORS

n_species_0 = n_species + dn;

n_species_css_0 = n_species_0(is_cond);
n_species_gss_0 = n_species_0(is_gas);

b_elements_css = formula_matrix(:,is_cond) * n_species_0(is_cond);
b_elements_gss = b_elements - b_elements_css;

x_isotopes_css = NaN; % delete later, not needed
x_isotopes_gss = NaN; % delete later, not needed

b_isotopes_cssNEW = b_isotopes_cssNEW + db_isotopes_cssNEW;
b_isotopes_gssNEW = b_isotopes_gssNEW + db_isotopes_gssNEW;


%% BREAK CONDITION

    % isLu = find(elements == 'Lu');    isAl = find(elements == 'Al');
    % if ( b_elements_css(isLu)/b_elements(isLu) < 1e-3 ) && ...
    %         (b_elements_css(isAl)/b_elements(isAl) < 1e-3 )
    %     break;
    % end


%% Re(calculate) Equilibrated Condensed SubSystems (CSS)

% Determine which species + elements to include in free energy minimization

[ ~, is_elementStableInCSS, ~, is_speciesPermissibleInCSS ] = ...
findStableSpeciesAndElementsInCss5_20240922( ...
    formula_matrix, ...                       % (1) formula matrix
    n_species_0, ...                          % (2) n_species
    b_elements, ...                           % (3) b_elements
    significance_threshold, ...               % (4) significance threshold
    is_cond, ...                              % (5) is_css
    phase_identifier, ...
    is_speciesInMultiSpeciesPhase, ...
    is_traceSpecies);

is_elements_css_excerpt = ~ ismember(elements,{'H','O'}) & is_elementStableInCSS;

% Inputs to free energy minimization function for condensed subsystem (CSS)


mu_0_species_css_excerpt = mu_0_species(is_speciesPermissibleInCSS,i);
A_css_excerpt = formula_matrix(is_elements_css_excerpt,is_speciesPermissibleInCSS);

b_css_excerpt = formula_matrix( ...
    is_elements_css_excerpt,...
    (is_speciesPermissibleInCSS | (is_traceSpecies & is_cond) ) ) ...
    * n_species_0(is_speciesPermissibleInCSS | (is_traceSpecies & is_cond) );      
% MODIFED - to allow major species in css to go to low numbers but not be
% set to zero

i_gas_css = [];
i_sacm_css = find(is_sacm(is_speciesPermissibleInCSS));  % 
i_sol_css = find(is_sol(is_speciesPermissibleInCSS));

% logical indicated which of 4 species in melt (SiO,Al2O3,CaO,MgO) to incl.
is_sacmSpeciesIncluded = is_speciesPermissibleInCSS(is_sacm);


% Find a feasible point to use as initial guess for gibbs minimization
[ n_species_css_0_excerpt, rel_violation ] = ...
    find_feasible_n_species_v5( ...
        A_css_excerpt, b_css_excerpt);



%%% Prepare input parameters to equilibrium calculation (MODIFIED)
equilibrium_params = struct();
equilibrium_params.n_species = n_species_css_0_excerpt;  % n_species is a dynamic variable
equilibrium_params.gamma_species = gamma_species(is_speciesPermissibleInCSS);  % gammas are a dynamic variable
equilibrium_params.formula_matrix = A_css_excerpt;  
equilibrium_params.b_elements = b_css_excerpt;
equilibrium_params.tau_sequence = tau_sequence;
equilibrium_params.thermodynamics = struct(...
    'R', R, ...
    'T', T(i), ...
    'P', P, ...
    'mu_0_species', mu_0_species_css_excerpt);  
equilibrium_params.cmas_melt = struct(...
    'W_sacm', W_sacm(:,i), ...
    'Q_m_sacm', Q_m_sacm);
equilibrium_params.indices = struct(...
    'i_gas', i_gas_css, ...
    'i_sacm', i_sacm_css, ...
    'i_sol', i_sol_css);
equilibrium_params.logicals = struct(...
    'is_sacmSpeciesIncluded', is_sacmSpeciesIncluded, ...
    'is_traceSpecies', is_traceSpecies(is_speciesPermissibleInCSS), ...
    'is_K2O_melt', is_K2O_melt(is_speciesPermissibleInCSS));
equilibrium_params.identifiers = struct(...
    'phase_identifier', phase_identifier(is_speciesPermissibleInCSS)); 



[n_species_css_excerpt, gammas_sacm_array(:,i), mu_species_cond_excerpt, ...
    log_act_species_cond_cut, zeta_vector_cond] = ...
        minimizeGibbsByInteriorPoint20240928( equilibrium_params );


%%% RECONCONSTRUCT SPECIES ABUNDANCE VECTOR "N_SPECIES_CSS" AND CHEMICAL
%%% POTENTIAL VECTOR "MU_SPECIES_CSS":

% Create species abundance vector of size num. of possible condensed
% species, and fill in values for "permissible" species — i.e., any
% condensed species (stable or not) that is comprised of elements that are
% present in the condensed system (i.e., not fully evaporated)
n_species_css = nan(sum(is_cond),1);
n_species_css(is_speciesPermissibleInCSS(is_cond)) = n_species_css_excerpt;

% For condensed species that are not "permissible" —
% i.e., that consist of elements that are fully evaporated — Set abundance
% equal to same value so that it won't chance within time step.  In other
% words, allow very small — practically zero — numbers to persist as such.
n_species_css(~is_speciesPermissibleInCSS(is_cond) & ~is_traceSpecies(is_cond)) = ...
    n_species_css_0(~is_speciesPermissibleInCSS(is_cond) & ~is_traceSpecies(is_cond));

% For condensed ***trace*** species, if they're not
% "permissible" — i.e., they occur in a phase that is not "permissible —
% then simply set to zero
n_species_css(~is_speciesPermissibleInCSS(is_cond) & is_traceSpecies(is_cond)) = 0;

% MODIFED - changed along above lines - to allow major 
% species in css to go to low numbers but not be set to zero
b_css_excerpt = ...
    formula_matrix( is_elements_css_excerpt,(is_speciesPermissibleInCSS ...
    | (is_traceSpecies & is_cond) ) ) ...
    * n_species_0(is_speciesPermissibleInCSS | (is_traceSpecies & is_cond) );       

% MODIFIED - changed NAN to Zeros - should be cleaned up
mu_species_css = zeros(sum(is_cond),1); % move into preallocate area            
mu_species_css(is_speciesPermissibleInCSS(is_cond)) = mu_species_cond_excerpt;

mu_species_css(~is_speciesPermissibleInCSS(is_cond) & is_sol(is_cond) ...
    & ~is_traceSpecies(is_cond)) = ...
    mu_0_species(is_cond & ~is_speciesPermissibleInCSS & is_sol & ...
    ~is_traceSpecies,i); % mu_0 used for mu for not-permissible species in CSS — true for pure phases, not true for melt; may want to revisit for melt


%% Equilibrate Gas SubSystem (GSS)

% Express formula matrix for GSS

A_gss = formula_matrix(:,is_gas);

% Find a feasible initial point for GSS


n_species_gss_0;  
% Prepare inputs for Gibbs minimization algorithm for GSS

mu_0_species_gss = mu_0_species(is_gas,i);
i_gas_gss = transpose(1:sum(is_gas));

% Run Gibbs minimization algorithm for GSS
[ n_species_gss, ~, mu_species_gss, ~, zeta_vector_gss ] = ...
    minimizeGibbsByInteriorPoint20240924( ...
        n_species_gss_0, ...
        A_gss, ...
        b_elements_gss, ...
        mu_0_species_gss, ...
        R, T(i), P, ...
        i_gas_gss, [], [], ...
        W_sacm, Q_m_sacm, ...
        tau_sequence(end), ...
        [], phase_identifier, [], is_traceSpecies(is_gas));

%%% RECONCONSTRUCT SPECIES ABUNDANCE VECTOR "N_SPECIES" AND CHEMICAL
%%% POTENTIAL VECTOR "MU_SPECIES"

n_species(is_cond) = n_species_css;
n_species(is_gas) = n_species_gss;

mu_species(is_cond) = mu_species_css;
mu_species(is_gas) = mu_species_gss;


%% CALCULATE SIZE OF PARTICLE AND GAS

%%% CONDENSED SYSTEM

% Mass (in grams) of elements in condensed system
mass_grams_css_elements = b_elements_css .* molar_mass_elements_g_mol;

% Total mass (in grams) of condensed system
mass_grams_css = sum(mass_grams_css_elements);

% Volume (in cm^3) of spherical particle
volume_cm3_css = mass_grams_css / density_particle_g_cm3;

% Radius (in cm) of spherical particle
radius_cm_css = (3/4 * volume_cm3_css / pi)^(1/3);

% Radius (in m) of spherical particle
radius_m_css = 1e-2 * radius_cm_css;

% Surface area (in m^2) of particle
surface_area_m2_css = 4 * pi * radius_m_css ^ 2;

%%% GAS SYSTEM

% Volume (in cL) of gas; NOTE: 1 cL = 1/100 L
volume_gas_cL = sum(n_species(is_gas)) * R * T(i) / P;         

% Volume (in m^3) of gas
volume_gas_m3 = volume_gas_cL * 1e-5; 

radius_gas_m = (3/4 * volume_gas_m3 / pi)^(1/3);     


%% CALCULATE SATURATION LEVELS AND FLUXES

[ is_speciesStableInCSS,~, ...
    ~, ...
    ~ ] = ...
    findStableSpeciesAndElementsInCss5_20240922( ...
    formula_matrix, ...                 % (1) formula matrix
    n_species, ...                      % (2) n_species
    b_elements, ...                     % (3) b_elements
    significance_threshold, ...         % (4) significance threshold
    is_cond, ...                        % (5) is_css
    phase_identifier, ...
    is_speciesInMultiSpeciesPhase, ...
    is_traceSpecies);    

%%% Calculate equilibrium (vapor) pressures and evaporative fluxes

% Inputs to equilibrium pressure and evaporative fluxes function

%%% Create Main Structure 'S', containing critical info for each time step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These logicals do not update with each iteration
constant_logicals.is_css = is_cond;
constant_logicals.is_gss = is_gas;
constant_logicals.is_sacm = is_sacm;
% constant_logicals.is_cationInSpecies = is_cationInSpecies;
constant_logicals.is_sol = is_sol;
constant_logicals.is_traceSpecies = is_traceSpecies;                           

% These logicals do (or may) update with each iteration
dynamic_logicals.is_stable = is_speciesStableInCSS;

% These are numerical tolerances
numerical_tolerances.tolerance_energyRxn = tolerance_energyRxn;
numerical_tolerances.tolerance_relativeFlux = tolerance_relativeFlux;

% These parameters do not update with each iteration
constant_params.formulaMatrix = formula_matrix;
constant_params.species = species;
constant_params.elements = elements;
constant_params.stickingCoeffs = sticking_coeffs;
constant_params.molMasses_kg_mol = molar_mass_species_kg_per_mol;
constant_params.R = R;
constant_params.dt = dt;
constant_params.substanceIdentifier = substance_identifier;  
constant_params.b_elements = b_elements;

% These parameters do (or may) update with each iteration
dynamic_params.pressure_bar = P;
dynamic_params.temperature_K = T(i);
dynamic_params.volumeGas_cL = volume_gas_cL;
dynamic_params.particleArea_m2 = surface_area_m2_css;
dynamic_params.n_species = n_species;
dynamic_params.mu_0_species = mu_0_species(:,i);
dynamic_params.mu_species = mu_species;
dynamic_params.n_saturation = n_saturation_array(:,i); % NaN
dynamic_params.mu_saturation = mu_sat_vector_array(:,i); % NaN

% The "simulation state" structure comprises the preceding structures
simulation_state.constant_logicals = constant_logicals;
simulation_state.dynamic_logicals = dynamic_logicals;
simulation_state.numerical_tolerances = numerical_tolerances;
simulation_state.constant_params = constant_params;
simulation_state.dynamic_params = dynamic_params;


%% CALCULATE VAPOR PRESSURES AND INSTANTANEOUS AND TIME-INTEGRATED FLUXES

% Calculate Instantaneous Fluxes
[ simulation_state.dynamic_params.n_saturation, ...
    simulation_state.dynamic_params.mu_saturation, ~, ...
    simulation_state.dynamic_params.rxnRatesInstant, ...
    simulation_state.dynamic_params.evaporationFluxesInstant, ...
    simulation_state.dynamic_params.productionFluxesInstant, ...
    simulation_state.dynamic_logicals.is_newlyStable, ...
    simulation_state.dynamic_logicals.P_envelope] ...
    = calculate_vapor_pressures_and_instantaneous_fluxes_20250303( ...
    simulation_state );  
% Calculate Instantaneous Fluxes
[ simulation_state.dynamic_params.productionFluxesTimeAveraged] = ...
    calculate_time_averaged_fluxes_20240915( simulation_state );


%% DISPLAY

showThese = is_gas | simulation_state.dynamic_logicals.is_stable;
[ species(showThese) n_species(showThese) ...
    simulation_state.dynamic_params.productionFluxesInstant(showThese) ]


%% CALCULATE ISOTOPIC EVAPORATIVE FLUXES

productionFluxesTimeAveraged = ...
    simulation_state.dynamic_params.productionFluxesTimeAveraged;
n_saturation_vector = simulation_state.dynamic_params.n_saturation;

[ db_dt_isotopes, ~, ~ ] = calculateEvaporativeIsotopeFluxes20240922( ...
    isElementIsoTrackedNEW, formula_matrix, n_species, n_saturation_vector, ...
    is_gas, isotopologueMassesNEW, molar_mass_species_g_per_mol,     ...
    productionFluxesTimeAveraged, b_isotopes_gssNEW, b_isotopes_cssNEW, ...
    b_elements_gss, b_elements_css);


%% Calculate Change in Species and Element Abundances For Time Step

dn_dt = simulation_state.dynamic_params.productionFluxesTimeAveraged;
dn = dn_dt * dt;
dn_css = dn(is_cond);
dn_gss = dn(is_gas);
db_css = formula_matrix(:,is_cond) * dn_css;
db_gss = formula_matrix(:,is_gas) * dn_gss;

db_isotopes_cssNEW = -(db_dt_isotopes) * dt; % NEW - NOVEMBER 10,2023
db_isotopes_gssNEW = db_dt_isotopes * dt; % NEW - NOVEMBER 10,2023


%% Save Values to Arrays

b_css_array(:,i) = b_elements_css;
b_gss_array(:,i) = b_elements_gss;

b_isotopes_css(i,:,:) =  b_isotopes_cssNEW;
b_isotopes_gss(i,:,:) =  b_isotopes_gssNEW;

% 
massParticle_g_array(i) = mass_grams_css;
radiusParticle_m_array(i) = radius_m_css;
areaParticle_m2_array(i) = surface_area_m2_css;

%
n_species_css_array(:,i) = n_species_css;
mu_species_css_array(:,i) = mu_species_css;

n_species_gss_array(:,i) = n_species_gss;
mu_species_gss_array(:,i) = mu_species_gss;

volume_gas_cL_array(i) = volume_gas_cL;

%
n_species_array(is_cond,i) = n_species_css_array(:,i);
n_species_array(is_gas,i) = n_species_gss_array(:,i);

%
mu_species_array(:,i) = mu_species;

%
savedProductionFluxes(:,i) = ...
    simulation_state.dynamic_params.productionFluxesTimeAveraged;
savedEvaporationFluxes(:,i) = ...
    simulation_state.dynamic_params.evaporationFluxesInstant;

dn_dt_gss_array(:,i) = ...
    simulation_state.dynamic_params.productionFluxesTimeAveraged;

n_saturation_array(:,i) = simulation_state.dynamic_params.n_saturation;
mu_sat_vector_array(:,i) = simulation_state.dynamic_params.mu_saturation;

P_envelope_array(i) = simulation_state.dynamic_logicals.P_envelope;


end 

toc

%% %%%%%%%%  END OF FORWARD FINITE DIFFERENCE (ITERATIVE) LOOP  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Plot model curves

% Simple figure 1
figure

% is_Ca_g = ismember(species, 'Ca_g');
isCa = find(elements == 'Ca');
isAl = find(elements == 'Al');

x = t;
y = ( b_gss_array(isAl,:) ./ b_elements(isAl) ) ./ ...
        ( b_gss_array(isCa,:) ./ b_elements(isCa) );

plot(x,y)
xlabel('time (s)','FontSize',12) 
ylabel('(Al/Ca)_{CI} of Gas','FontSize',12) 

% Simple figure 2
figure

x = t;
y = 1e3 * (( b_isotopes_gss(:,isCa,2) ./ b_isotopes_gss(:,isCa,1) ) ./ ...
    ( b_isotopesNEW(isCa,2) ./ b_isotopesNEW(isCa,1) ) - 1);

plot(x,y)
xlabel('time (s)','FontSize',12) 
ylabel('{\delta}^{44/40}Ca of Gas','FontSize',12) 
