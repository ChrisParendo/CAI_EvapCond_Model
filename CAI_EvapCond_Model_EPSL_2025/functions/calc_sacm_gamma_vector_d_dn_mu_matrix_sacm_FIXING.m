function [gamma_sacm, R_T_log_gamma_sacm,d_dn_mu_sacm] = ...
    calc_sacm_gamma_vector_d_dn_mu_matrix_sacm_FIXING( ...
        W_sacm, Q_m_sacm, X_sacm_input, d_dn_invX, xxxx_sacm, d_dn_xxxx_sacm, R,T, nt_sacm)
%calc_sacm_gamma_vector_d_dn_mu_matrix_sacm
%   
num_sacm = 4;
% rename:
X = X_sacm_input;
% d_dn_X = d_dn_x_sacm;
XXXX = xxxx_sacm;
d_dn_XXXX = d_dn_xxxx_sacm;

% preallocate:
R_T_log_gamma_sacm = nan(4,1);
d_dn_R_T_log_gamma_sacm = nan(4);

for m = 1:4

    % calculate R*T*log(gamma) vector for SACM species:
    
    temp_array0 = ( Q_m_sacm(:,m) * 1/X(m) - 3 ) .* XXXX;

    R_T_log_gamma_m = transpose(W_sacm) * temp_array0; % true matrix mult.

    R_T_log_gamma_sacm(m) = R_T_log_gamma_m;

    % calculate Jacobian for chemical potential vector for SACM species —
    % i.e., obtain 4x4 derivative matrix for mu vector with respect to n_i
    % species:

    temp_array1 = Q_m_sacm(:,m) .* d_dn_invX(m,:) .* XXXX;

    temp_array2 = ( Q_m_sacm(:,m) * 1/X(m) - 3 ) .* d_dn_XXXX;

    temp_array3 = temp_array1 + temp_array2;

    d_dn_R_T_log_gamma_m = transpose(W_sacm) * temp_array3;

    d_dn_R_T_log_gamma_sacm(m,:) = d_dn_R_T_log_gamma_m;



% %     temp_array1 = ( Q_m_sacm(:,m) + X(m) ) .* d_dn_XXXX;
% % 
% %     temp_array2 = (-3) * d_dn_X(m,:) .* XXXX; % -3 = 1-p
% % 
% %     temp_array3 = temp_array1 + temp_array2;
% % 
% %     d_dn_mu_m = R*T * ( transpose(W_sacm) * temp_array3 ); % true matrix mult.
% % 
% %     d_dn_mu_sacm(m,:) = d_dn_mu_m;

    
    

end

% calculate gamma vector for SACM species:
gamma_sacm = exp( (R*T)^-1 * R_T_log_gamma_sacm);


% d_dn_R_T_log_x_sacm = R*T*(repmat(-1/sum(n_species(i_sacm)), num_sacm, num_sacm) + diag(1./n_species(i_sacm)));
d_dn_R_T_log_x_sacm = R*T*( repmat(-1/nt_sacm, num_sacm, num_sacm) + diag(1./(X_sacm_input*nt_sacm)) );


d_dn_mu_sacm = d_dn_R_T_log_x_sacm + d_dn_R_T_log_gamma_sacm;

end
