function [T_SolarSys_Abundances] = ...
            load_SolarSys_Abundances_Lodders2003()

% Load solar system elemental abundances from Lodders (2003)
file = fullfile('data', 'Lodders_2003_Table6.xlsx');
sheet = 'Lodders_2003_Table6_2';
range = 'A1:H84';

T_SolarSys_Abundances = readtable(file,'Sheet',sheet,'Range',range);
T_SolarSys_Abundances.Properties.RowNames = T_SolarSys_Abundances.Element;

end