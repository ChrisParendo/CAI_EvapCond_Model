function [W_array,Q_m_array] = ...
    calc_Berman1983_W_parameters_and_Qm(T)
%calc_Berman1983_W_parameters_and_Qm

%%% Load W_H, W_S, and (less importantly) table with indices for polynomial
%   terms

% Load Margules parameters: W_H and W_S
file_name = ['/Users/chrisparendo/Documents/Hub_Literature_Data/' ...
    'Berman1983/Berman_Margules_Params.xlsx'];
sheet_name = 'Wparams_TableXI';
range_name = 'A1:F32';
T_W_params = readtable(file_name,'sheet',sheet_name,'range',range_name);

% Assign values for Margules Parameters
W_H = T_W_params.W_H;
W_S = T_W_params.W_S;


%%% Calculate W_{Gibbs} vector, simply represented as W, as function of T
%   That is, calculate margules parameters (W's) as function of T -
%    an array comprising W column vectors for each input temperautre

W_array = W_H - T .* W_S;


%%% Create Q_m vectors for each of species in SACM melt: Place vectors into
%   single array: yields a 4x31 array, where each column vector is Q_m as
%   defined in Berman1983 for one of the 4 components (S A C M)

indices = T_W_params{:,["S", "A", "C", "M"]};

Q_m_array = nan(31,4);
for m = 1:4
    Q_m_array(:,m) = sum( (indices == m) , 2);
end


end
