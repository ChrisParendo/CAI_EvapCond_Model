function [Gibbs_0_T_species] = loadGibbs0Tspecies_v4_Dec03(species,T)   % Nov. 19 version uses load_Gibbs_0_T_NIST_v2
%loadGibbs0Tspecies Loads data (by calling additional functions) and organizes ...
%   Function calls additional functions (e.g., loadGibbs0Tnist) in order to
%   load data from varous sources (e.g., NIST), and subsequently organizes
%   data into single array with standard Gibbs energies of species (rows)
%   at various temperature(s) (columns)


%%% Load data

% Load/calculate Gibbs energies for species from NIST webbook
[species_NIST, Gibbs_0_T_species_NIST]  = load_Gibbs_0_T_NIST_v2(T);
            
% Load/calculate Gibbs energies for species from Berman (1983) thesis
[species_Berman, Gibbs_0_T_species_Berman] = load_Gibbs_0_T_Berman(T);

% % % Placeholder for manuscript submission

species_REEsolids = [];                 % Placeholder for manuscript submission
Gibbs_0_T_species_REEsolids = [];       % Placeholder for manuscript submission
species_REE_liquids = [];               % Placeholder for manuscript submission
Gibbs_0_T_species_REE_liquids = [];     % Placeholder for manuscript submission


% % % Placeholder for manuscript submission

species_REEAtomsGas = [];               % Placeholder for manuscript submission
Gibbs_0_T_REEAtomsGas = [];             % Placeholder for manuscript submission
species_REEOxides_g = [];               % Placeholder for manuscript submission
Gibbs_0_T_REEOxides_g = [];             % Placeholder for manuscript submission



% Load/calculate Gibbs energies for Ti bearing species

[ species_TiMinerals, Gibbs_0_T_species_TiMinerals ] = ...
    load_Gibbs_0_T_TiMinerals(T);


% % % Placeholder for manuscript submission

species_REEaluminates = [];                  % Placeholder for manuscript submission
Gibbs_0_T_species_REEaluminates = [];        % Placeholder for manuscript submission


%%% Organize Gibbs energies of species data
%   to create Gibbs_0_T_species array

% Concatenate data from all sources

species_AllSources = [species_NIST; species_Berman; species_REEsolids; species_REE_liquids; ...
    species_REEAtomsGas; species_REEOxides_g; species_TiMinerals; species_REEaluminates];

Gibbs_0_T_species_AllSources = ...
    [Gibbs_0_T_species_NIST; Gibbs_0_T_species_Berman; ...
    Gibbs_0_T_species_REEsolids; Gibbs_0_T_species_REE_liquids; Gibbs_0_T_REEAtomsGas; Gibbs_0_T_REEOxides_g; Gibbs_0_T_species_TiMinerals; Gibbs_0_T_species_REEaluminates ];

% Create array of Gibbs_0_T_species for system of interest —
% with the species as rows and the temperatures as columns
    
% Pre-allocate space in array prior to 'for loop'
Gibbs_0_T_species = nan(length(species),length(T));

                                % Fill array with values
                                for i = 1:length(species)
                                    % index = ismember(species_AllSources,species(i));      
                                    % may want to change the following in the future
                                    index = false(size(species_AllSources));
                                    for j = 1:length(species_AllSources)
                                        index(j) = ismember(species(i),species_AllSources(j));      % changed "contains to ismember Sept. 29,2024 /// confusing! may cause bugs in old versions of code ... set to "ismember" from "contains" on Jan 29, 2025
                                    end
                                    Gibbs_0_T_species(i,:) = Gibbs_0_T_species_AllSources(index,:);
                                end    


end


