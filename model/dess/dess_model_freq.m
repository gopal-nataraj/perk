%% Single-Isochromat Frequency Analysis for DESS
% Written by: Gopal Nataraj

% Constant declarations
M0 = 1;             
T1 = 500;               % ms
T2 = 50;                % ms
alpha = 45*(pi/180);    % rad
TR = 20;                % ms
TE = 5;                 % ms

% Intermediate constants
E1R = exp(-TR./T1);
E2R = exp(-TR./T2);
p = (1 - E1R.*cos(alpha)) - (E1R - cos(alpha)).*(E2R.^2);
q = E2R .* (1-E1R) .* (1 + cos(alpha));

% Compute M+ and M- over different values of beta (off-resonance effects)
beta = (-180:180)' * (pi/180);
Qbeta = (M0*(1-E1R)) ./ (p - q*cos(beta));
Mp = (+1i*Qbeta*sin(alpha)) .* (1 - E2R*exp(1i*beta));
Mm = (-1i*Qbeta*sin(alpha)) .* (E2R^2 - E2R*exp(-1i*beta));

% Project onto the subspace [1 cos(beta) sin(beta)] to visually
% see if most of the energy is in the zeroth and first fundamentals.
basis = [ones(length(beta), 1) cos(beta) sin(beta)];
coef_p = (1/pi) * [0.5 1 1] .* trapz(beta, repmat(conj(Mp), [1 3]) .* basis); 
coef_m = (1/pi) * [0.5 1 1] .* trapz(beta, repmat(conj(Mm), [1 3]) .* basis); 
proj_p = basis * coef_p';
proj_m = basis * coef_m'; 

% Plot real and imaginary components 
figure; hold on;
    plot(beta*180/pi, real(Mp), 'b'); 
    plot(beta*180/pi, imag(Mp), 'r'); 
    plot(beta*180/pi, real(proj_p), 'b--');
    plot(beta*180/pi, imag(proj_p), 'r--'); hold off;
    set(gca, 'Xtick', [-180:45:180]); axis([-180 180 -0.1 0.2]);
    xlabel('Cumulative off-resonance, beta (deg)');
    ylabel('First echo isochromat signal, Mp_x_y(beta)');
    title('Single-isochromat First-Echo Signal and Fund. Frequency Projection');
    legend('Real(Mp_x_y(beta))', 'Imag(Mp_x_y(beta))',...
        'Real(FF projection)', 'Imag(FF projection)');
%     print('-depsc', 'Mp_beta.eps');
    
figure; hold on;
    plot(beta*180/pi, real(Mm), 'b'); 
    plot(beta*180/pi, imag(Mm), 'r'); 
    plot(beta*180/pi, real(proj_m), 'b--');
    plot(beta*180/pi, imag(proj_m), 'r--'); hold off;
    set(gca, 'Xtick', [-180:45:180]); axis([-180 180 -0.15 0.1]);
    xlabel('Cumulative off-resonance, beta (deg)');
    ylabel('Second echo isochromat signal, Mm_x_y(beta)');
    title('Single-isochromat Second-Echo Signal and Fund. Frequency Projection');
    legend('Real(Mm_x_y(beta))', 'Imag(Mm_x_y(beta))',...
        'Real(FF projection)', 'Imag(FF projection)');
%     print('-depsc', 'Mm_beta.eps');
    