% Header File -- Min-max CRLB analysis for T1 estimation from RF SE IR sequences
% Worst-case M0s/T1 CRLB over T1/T2 ROI, varying inversion times TI
% In this version, hold TR values fixed and equal
%
% Written by: Gopal Nataraj
% Originally created:   2015-05-07
% Last modified:        2015-06-25

% Set imaging constraints 
TR_tot = 6000;                          % ms (TR_ir_min = 3000ms)
TI_min = 50;                            % ms
TI_res = 50;                            % ms

% Implicit parameters, including a T1/kappa range of interest
T1 = linspace(800, 1400, 5);            % WM/GM voxels of interest
kap = linspace(0.9, 1.1, 5);            % a.u. (10% variation)
wf = 0;                                 % rad/ms

% Outer robustness criterion: check minima over loose bounds
T1_rob  = linspace(400, 2000, 9);       % ms
kap_rob = 2 .^ linspace(-1, 1, 5);      % a.u. (factor-of-two variation)

% Other constant declarations
TE = 14;                                % ms (minimum, fixed)
a_inv = pi;                             % inversion pulse flip (rad)
a_ex = pi/2;                            % excitation pulse flip (rad)
noise_var_1coil = 2.62e-07;             % a.u.; computed at 1x1x5 resolution
noise_var_ssos = noise_var_1coil/2;     % High-SNR regime approximation

w = [0 1]';                             % [M02 T1]' relative weighting
tol = 0.01;                             % Global minima tolerance
time_comp = 0;                          % Toggle time compensation on/off
pr = 0;                                 % Toggle print on/off
savedat = 0;                            % Toggle saving .mat files on/off
T1range = [0 10];                       % sigmaT1 plot range

% Choose number of IR scans
n_ir = 2;                               % Number of IR scans
profile = sprintf('%u IR', n_ir);

% Introduce the experiment and start timer
fprintf('\n\nScan Profile: (%u IR)\n', n_ir);
fprintf('Total TR constraint: %1.2f ms\n', TR_tot);
fprintf('Weightings: (%1.1f,%1.1f)\n', w(1), w(2));
tic;

switch profile
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scan Sequence One: 2 IR %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '2 IR'
        % Controllable parameter grid: varying TI
        % Allow TI to go up to no more than TR-2*TE
        TR = (TR_tot/2)  * ones(2, 1);
        TI_max = (TR_tot/2) - ceil(2*TE/10)*10; 
        TI1_ir = [TI_min : TI_res : TI_max];
        TI2_ir = [TI_min : TI_res : TI_max];
        flip_inv = a_inv * ones(2, 1);
        flip_ex  = a_ex  * ones(2, 1);
        
        % Grid definitions
        sig_T1_2    = inf(length(TI1_ir), length(TI2_ir));
        sigw_T1_2   = inf(length(TI1_ir), length(TI2_ir));
        Psi_2       = inf(length(TI1_ir), length(TI2_ir));
        
        % Covariance matrix: scaled identity with D = 2
        Sigma_inv = (1 ./ noise_var_ssos) * speye(2);
        
        %% Step One: Min-max CRLB analysis over tight parameter grid
        TI_ir_2 = [];
        for i_TI1 = 1:length(TI1_ir)
        for i_TI2 = 1:length(TI2_ir)
            TI_ir_2 = [TI1_ir(i_TI1) TI2_ir(i_TI2)]';
            
            % Inner maximization: store max SDs over tight (T1,kap) range
            worst_sigM02 = NaN(length(T1), length(kap));
            worst_sigT1  = NaN(length(T1), length(kap));
            
            for i_t1 = 1:length(T1)
                for i_kap = 1:length(kap)
                    [~, worst_sigM02(i_t1, i_kap), worst_sigT1(i_t1, i_kap)] = ...
                        norm_crlb_IR_2parm(T1(i_t1), wf, kap(i_kap)*flip_inv,...
                        kap(i_kap)*flip_ex, TR, TI_ir_2, TE, Sigma_inv, time_comp);
                end
            end
            
            % Store the worst-case sig_T1 value
            sigw_T1_2(i_TI1, i_TI2) = max(worst_sigT1(:));
            
            % Store the worst-case Psi and its corresponding T1 index
            worst_Psi = w(1)*worst_sigM02 + w(2)*worst_sigT1;
            [Psi_2(i_TI1, i_TI2), idx_Psiw] = max(worst_Psi(:));
            
            % Extract this index and use it to save corresponding sigmas
            [T1_idx, kap_idx] = ind2sub(size(worst_Psi), idx_Psiw);
            sig_T1_2(i_TI1, i_TI2) = worst_sigT1(T1_idx, kap_idx);
        end
        end
        
        % Find the indices of (multiple) Psi_2 minima, to within a tolerance
        Psimin_idx_2 = find( (Psi_2 - min(Psi_2(:))) ./ min(Psi_2(:)) <= tol);
        num_min_2 = length(Psimin_idx_2);
        
        %% Step Two: Select one minimum based on robustness over broad (T1, kap) range
        % Indices of each of the "num_min_2" total minima
        TI1_idx_2 = NaN(num_min_2, 1);
        TI2_idx_2 = NaN(num_min_2, 1);
        
        sigM02_worst_2 = NaN(num_min_2, length(T1_rob), length(kap_rob));
        sigT1_worst_2  = NaN(num_min_2, length(T1_rob), length(kap_rob));
        
        for i_min = 1:num_min_2
            % Convert the 1D index to ND-grid indices
            [TI1_idx_2(i_min), TI2_idx_2(i_min)] = ...
                ind2sub(size(Psi_2), Psimin_idx_2(i_min));
            
            % Store the minimizing parameters
            TI1_min_2 = TI1_ir(TI1_idx_2(i_min));
            TI2_min_2 = TI2_ir(TI2_idx_2(i_min));
            
            % Evaluate CRLB at minimizers, over (T1_rob, kap_rob) ranges
            for i_t1 = 1:length(T1_rob)
                for i_kap = 1:length(kap_rob)
                    [~, sigM02_worst_2(i_min, i_t1, i_kap),...
                        sigT1_worst_2(i_min, i_t1, i_kap)] = ...
                        norm_crlb_IR_2parm(T1_rob(i_t1), wf, kap_rob(i_kap)*flip_inv,...
                        kap_rob(i_kap)*flip_ex, TR, [TI1_min_2 TI2_min_2]',...
                        TE, Sigma_inv, time_comp);
                end
            end
        end
        
        % Store the worst-case sigT1
        sigM02_min_2 = max(reshape(sigM02_worst_2,...
            [num_min_2 length(T1_rob)*length(kap_rob)]), [], 2);
        sigT1_min_2  = max(reshape(sigT1_worst_2,...
            [num_min_2 length(T1_rob)*length(kap_rob)]), [], 2);
        
        % Compute the worst-case Psi_2 value over the wider range
        Psi_min_2 = w(1)*sigM02_worst_2 + w(2)*sigT1_worst_2;
        [Psi_maxmin_2, wide_idx_2] = max(reshape(Psi_min_2,...
            [num_min_2 length(T1_rob)*length(kap_rob)]), [], 2);
        
        % Extract the parameter indices that minimize Psi_maxmin_2 over minima
        % Psi_star - min(Psi(:)) is the degradation from narrow->wide range
        [Psi_star_2, i_star_2] = min(Psi_maxmin_2);
        sigT1_star_2 = sigT1_min_2(i_star_2);
        
        TI1star_2 = TI1_idx_2(i_star_2);
        TI2star_2 = TI2_idx_2(i_star_2);
        
        % Display results
        fprintf('Selected scan design: (TI1, TI2) = (%1.1f ms, %1.1f ms)\n',...
            TI1_ir(TI1star_2), TI2_ir(TI2star_2));
        fprintf('\t Worst-case sigT1: %6.2f\n', sigw_T1_2(TI1star_2, TI2star_2));
        fprintf('\t Robust-range sigT1: %6.2f\n', sigT1_star_2);
        fprintf('\t Worst-case Psi: %6.2f\n', Psi_2(TI1star_2, TI2star_2));
        fprintf('\t Robust-range Psi: %6.2f\n', Psi_star_2);
        fprintf('\t Selected from %d candidates, at a %d%% tolerance, in %5.0f seconds.\n\n',...
            num_min_2, 100*tol, toc); 
        
        %% Graphical Output: 2D Heatmaps for 2 IR Scans
        % Minimizer location(s)
        xmin = [TI1_ir(TI1star_2) TI2_ir(TI2star_2)];
        ymin = [TI2_ir(TI2star_2) TI1_ir(TI1star_2)];
        
        % T1 sigma vs. (TI1/TI2)
        figure; hold on; imagesc(TI1_ir, TI2_ir, sigw_T1_2', T1range);...
            scatter(xmin, ymin, 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');...
            xlabel('IR TI_1 (ms)', 'FontSize', 14);...
            ylabel('IR TI_2 (ms)', 'FontSize', 14);...
            if(~pr), title('Worst-case T1 Standard Deviation (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(TI1_ir)' minmax(TI2_ir)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('sigw_t1_2', suffix', '.eps')), end;
        
        % Worst-case Psi vs. (TI1/TI2)
        figure; hold on; imagesc(TI1_ir, TI2_ir, Psi_2', T1range);...
            scatter(xmin, ymin, 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');...
            xlabel('IR TI_1 (ms)', 'FontSize', 14);...
            ylabel('IR TI_2 (ms)', 'FontSize', 14);...
            if(~pr), title('Worst-case Psi (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(TI1_ir)' minmax(TI2_ir)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('psiw_2', suffix', '.eps')), end;
        
    otherwise
        warning('Unexpected scan profile!');
end     