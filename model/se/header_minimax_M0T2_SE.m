% Header File -- Min-max CRLB analysis for T2 estimation from SE sequences
% Worst-case M0/T2 CRLB over T1/T2 ROI, varying echo times TE
% 
% Assumptions:
% 1) TR values >> T1 (complete recovery)
% 2) TR values fixed and equal to each other
% 3) Inversion pulse perfectly 180 degrees and not spatially varying
% 
% Written by: Gopal Nataraj
% Originally created: 2015-05-08

% Set imaging constraints
TR_tot = 6000;                          % ms 
TE_min = 10;                            % ms
TE_res = 5;                             % ms

% Implicit parameters, including a T2/kappa range of interest
T2 = linspace(50, 120, 5);              % WM/GM voxels of interest
kap = linspace(0.9, 1.1, 5);            % a.u. (10% variation)
wf = 0;                                 % rad/ms

% Outer robustness criterion: check minima over loose bounds
T2_rob  = linspace(40,  200,  9);       % ms
kap_rob = 2 .^ linspace(-1, 1, 5);      % a.u. (factor-of-two variation)

% Other constant declarations
a_ex = pi/2;                            % excitation pulse flip (rad)
noise_var_1coil = 2.62e-07;             % a.u.; computed at 1x1x5 resolution
noise_var_ssos = noise_var_1coil/2;     % High-SNR regime approximation

w = [0 1]';                             % [M0 T2]' relative weighting
tol = 0.01;                             % Global minima tolerance
time_comp = 0;                          % Toggle time compensation on/off
pr = 0;                                 % Toggle print on/off
savedat = 0;                            % Toggle saving .mat files on/off
T2range = [0 1];                        % sigmaT2 plot range

% Choose number of SE scans
n_se = 2;                               % Number of SE scans
profile = sprintf('%u SE', n_se);

% Introduce the experiment and start timer
fprintf('\n\nScan Profile: (%u SE)\n', n_se);
fprintf('Total TR constraint: %1.2f ms\n', TR_tot);
fprintf('Weightings: (%1.1f,%1.1f)\n', w(1), w(2));
tic;

switch profile
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scan Sequence One: 2 SE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '2 SE'
        % Controllable parameter grid: varying TE
        % Allow TE to go up to no more than 2*T2max of interest
        TR = (TR_tot/2) * ones(2, 1);
        TE_max = 2*max(T2);
        TE1_se = [TE_min : TE_res : TE_max];
        TE2_se = [TE_min : TE_res : TE_max];
        flip_ex = a_ex * ones(2, 1);
        
        % Grid definitions
        sig_T2_2    = inf(length(TE1_se), length(TE2_se));
        sigw_T2_2   = inf(length(TE1_se), length(TE2_se));
        Psi_2       = inf(length(TE1_se), length(TE2_se));
        
        % Covariance matrix: scaled identity with D = 2
        Sigma_inv = (1 ./ noise_var_ssos) * speye(2);
        
        %% Step One: Min-max CRLB analysis over tight parameter grid
        TE_se_2 = [];
        for i_TE1 = 1:length(TE1_se)
        for i_TE2 = 1:length(TE2_se)
            TE_se_2 = [TE1_se(i_TE1) TE2_se(i_TE2)]';
            
            % Inner maximization: store max SDs over tight (T2,kap) range
            worst_sigM0 = NaN(length(T2), length(kap));
            worst_sigT2 = NaN(length(T2), length(kap));
    
            for i_t2 = 1:length(T2)
                for i_kap = 1:length(kap)
                    [~, worst_sigM0(i_t2, i_kap), worst_sigT2(i_t2, i_kap)] = ...
                        norm_crlb_SE_2parm(T2(i_t2), kap(i_kap)*flip_ex,...
                        TR, TE_se_2, Sigma_inv, time_comp);
                end
            end
            
            % Store the worst-case sig_T2 value
            sigw_T2_2(i_TE1, i_TE2) = max(worst_sigT2(:));
            
            % Store the worst-case Psi and its corresponding T2 index
            worst_Psi = w(1)*worst_sigM0 + w(2)*worst_sigT2;
            [Psi_2(i_TE1, i_TE2), idx_Psiw] = max(worst_Psi(:));
            
            % Extract this index and use it to save corresponding sigmas
            [T2_idx, kap_idx] = ind2sub(size(worst_Psi), idx_Psiw);
            sig_T2_2(i_TE1, i_TE2) = worst_sigT2(T2_idx, kap_idx);
        end
        end
        
        % Find the indices of (multiple) Psi_2 minima, to within a tolerance
        Psimin_idx_2 = find( (Psi_2 - min(Psi_2(:))) ./ min(Psi_2(:)) <= tol);
        num_min_2 = length(Psimin_idx_2);
    
        %% Step Two: Select one minimum based on robustness over broad (T2, kap) range
        % Indices of each of the "num_min_2" total minima
        TE1_idx_2 = NaN(num_min_2, 1);
        TE2_idx_2 = NaN(num_min_2, 1);
        
        sigM0_worst_2 = NaN(num_min_2, length(T2_rob), length(kap_rob));
        sigT2_worst_2 = NaN(num_min_2, length(T2_rob), length(kap_rob));
        
        for i_min = 1:num_min_2
            % Convert the 1D index to ND-grid indices
            [TE1_idx_2(i_min), TE2_idx_2(i_min)] = ...
                ind2sub(size(Psi_2), Psimin_idx_2(i_min));
            
            % Store the minimizing parameters
            TE1_min_2 = TE1_se(TE1_idx_2(i_min));
            TE2_min_2 = TE2_se(TE2_idx_2(i_min));
            
            % Evaluate CRLB at minimizers, over (T2_rob, kap_rob) ranges
            for i_t2 = 1:length(T2_rob)
                for i_kap = 1:length(kap_rob)
                    [~, sigM0_worst_2(i_min, i_t2, i_kap),...
                        sigT2_worst_2(i_min, i_t2, i_kap)] = ...
                        norm_crlb_SE_2parm(T2_rob(i_t2), kap_rob(i_kap)*flip_ex,...
                        TR, [TE1_min_2 TE2_min_2]', Sigma_inv, time_comp);
                end
            end
        end
        
        % Store the worst-case sigT2
        sigM0_min_2 = max(reshape(sigM0_worst_2,...
            [num_min_2 length(T2_rob)*length(kap_rob)]), [], 2);
        sigT2_min_2 = max(reshape(sigT2_worst_2,...
            [num_min_2 length(T2_rob)*length(kap_rob)]), [], 2);
        
        % Compute the worst-case Psi_2 value over the wider range
        Psi_min_2 = w(1)*sigM0_worst_2 + w(2)*sigT2_worst_2;
        [Psi_maxmin_2, wide_idx_2] = max(reshape(Psi_min_2,...
            [num_min_2 length(T2_rob)*length(kap_rob)]), [], 2);
        
        % Extract the parameter indices that minimize Psi_maxmin_2 over minima
        % Psi_star - min(Psi(:)) is the degradation from narrow->wide range
        [Psi_star_2, i_star_2] = min(Psi_maxmin_2);
        sigT2_star_2 = sigT2_min_2(i_star_2);
        
        TE1star_2 = TE1_idx_2(i_star_2);
        TE2star_2 = TE2_idx_2(i_star_2);
        
        % Display results
        fprintf('Selected scan design: (TE1, TE2) = (%1.1f ms, %1.1f ms)\n',...
            TE1_se(TE1star_2), TE2_se(TE2star_2));
        fprintf('\t Worst-case sigT2: %6.2f\n', sigw_T2_2(TE1star_2, TE2star_2));
        fprintf('\t Robust-range sigT2: %6.2f\n', sigT2_star_2);
        fprintf('\t Worst-case Psi: %6.2f\n', Psi_2(TE1star_2, TE2star_2));
        fprintf('\t Robust-range Psi: %6.2f\n', Psi_star_2);
        fprintf('\t Selected from %d candidates, at a %d%% tolerance, in %5.0f seconds.\n\n',...
            num_min_2, 100*tol, toc); 
    
        %% Graphical Output: 2D Heatmaps for 2 SE Scans
        % Minimizer location(s)
        xmin = [TE1_se(TE1star_2) TE2_se(TE2star_2)];
        ymin = [TE2_se(TE2star_2) TE1_se(TE1star_2)];
        
        % T2 sigma vs. (TE1/TE2)
        figure; hold on; imagesc(TE1_se, TE2_se, sigw_T2_2', T2range);...
            scatter(xmin, ymin, 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');...
            xlabel('SE TE_1 (ms)', 'FontSize', 14);...
            ylabel('SE TE_2 (ms)', 'FontSize', 14);...
            if(~pr), title('Worst-case T2 Standard Deviation (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(TE1_se)' minmax(TE2_se)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('sigw_t2_2', suffix', '.eps')), end;
        
        % Worst-case Psi vs. (TE1/TE2)
        figure; hold on; imagesc(TE1_se, TE2_se, Psi_2', T2range);...
            scatter(xmin, ymin, 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');...
            xlabel('SE TE_1 (ms)', 'FontSize', 14);...
            ylabel('SE TE_2 (ms)', 'FontSize', 14);...
            if(~pr), title('Worst-case Psi (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(TE1_se)' minmax(TE2_se)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('psiw_2', suffix', '.eps')), end;
        
    otherwise
        warning('Unexpected scan profile!');
end     