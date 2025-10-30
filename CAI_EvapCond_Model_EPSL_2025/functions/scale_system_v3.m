function [b_elements_scaled, scaling_factor, mass_particle_g] = ...
    scale_system_v3( ...
    diameter_particle_mm, density_particle_g_cm3, ...
    b_elements, b_elements_css, ...
    molar_mass_elements_g_mol)
%scale_system



%%% Calculate moles of each element per gram of condensed material

% element abundances (moles) in condensed material (prior to scaling): 
moles_elements = b_elements_css;

% element mass abundances (grams) in condensed material (prior to scaling): 
grams_elements = molar_mass_elements_g_mol .* moles_elements;

% total mass of condensed material (prior to scaling):
grams_condensate = sum(grams_elements);

% % (moles of each element, vector) / (grams of condensed material, scalar):
% molesElements_gramCondensate = moles_elements / grams_condensate;



%%% Calculate mass of particle in grams (at prescribed scale)

% calculate particle volume (cm^3) - assuming it's spherical
diameter_particle_cm = 1e-1 * diameter_particle_mm;
radius_particle_cm = 1/2 * diameter_particle_cm;
volume_particle_cm3 = 4/3 * pi * (radius_particle_cm)^3 ;

% calculate particle mass (g)
mass_particle_g = density_particle_g_cm3 * volume_particle_cm3;



%%% Calculate scaling factor

% scaling factor:
scaling_factor = mass_particle_g / grams_condensate;



%%% Scale element and species abundance vectors (applies to both gas and
%%% condensed materials)
% b_elements_scaled = b_elements * scaling_factor;
b_elements_scaled = b_elements * scaling_factor;

end