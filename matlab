clc; clear; close all;

% System Parameters
lambda = 1550e-9;    
k = 2 * pi / lambda; 
Cn2 = 5e-15;          % Moderate turbulence

D_T = 0.1;           
D_R = 0.2;           
theta = 1e-3;        

atten_dB_per_km = 5; 
atten_linear = @(L) 10.^(-atten_dB_per_km*L/1000/10);

B = 1e9;             
N0 = 1e-19;          
M = 16;              
num_samples = 1e4;   

% Distance range
distance = 50:100:1500; 
BER = zeros(length(distance), 4); 

% SNR reference scaling
Target_SNR_dB_at_100m = 13;  
Target_SNR_linear = 10^(Target_SNR_dB_at_100m/10);
Ptx = 1e-3;  

for i = 1:length(distance)
    L = distance(i);
    
    % Calculate turbulence strength (Rytov variance)
    sigmaR2 = 1.23 * Cn2 * k^(7/6) * L^(11/6);
    
    % Gamma-Gamma parameters
    alpha = (exp(0.49*sigmaR2/(1+1.11*sigmaR2^(6/5))) - 1)^(-1);
    beta  = (exp(0.51*sigmaR2/(1+0.69*sigmaR2^(6/5))) - 1)^(-1);
    
    % Total losses
    geom_loss = (D_R / (D_T + theta*L))^2;
    total_loss = atten_linear(L) * geom_loss;
    Pr = Ptx * total_loss;
    
    % Normalize SNR
    Pr_ref = 1e-3 * atten_linear(100) * (D_R / (D_T + theta*100))^2;
    SNR_scale = Target_SNR_linear / (Pr_ref / (N0*B));
    SNR_samples = (Pr .* gamrnd(alpha,1/alpha,[num_samples 1]) .* gamrnd(beta,1/beta,[num_samples 1])) * SNR_scale / (N0*B);
    
    % Reference BER (without coding)
    ber_ref = (4/log2(M)) .* (1-1/sqrt(M)) .* qfunc(sqrt(3*SNR_samples/(M-1)));
    BER(i,1) = mean(ber_ref);
    
    % LDPC model: ~4 dB gain
    SNR_LDPC = SNR_samples * 10^(4/10);
    ber_ldpc = (4/log2(M)) .* (1-1/sqrt(M)) .* qfunc(sqrt(3*SNR_LDPC/(M-1)));
    BER(i,2) = mean(ber_ldpc);
    
    % Gray coding model: ~2 dB gain
    SNR_gray = SNR_samples * 10^(2/10);
    ber_gray = (4/log2(M)) .* (1-1/sqrt(M)) .* qfunc(sqrt(3*SNR_gray/(M-1)));
    BER(i,3) = mean(ber_gray);
    
    % Hybrid Gray+LDPC model: ~6 dB gain
    SNR_hybrid = SNR_samples * 10^(6/10);
    ber_hybrid = (4/log2(M)) .* (1-1/sqrt(M)) .* qfunc(sqrt(3*SNR_hybrid/(M-1)));
    BER(i,4) = mean(ber_hybrid);
end

% Plotting
figure; hold on; grid on; box on;

semilogy(distance, BER(:,1), '-o','LineWidth', 3, 'MarkerSize',7, 'Color','b');
semilogy(distance, BER(:,2), '--s','LineWidth', 3, 'MarkerSize',7, 'Color','r');
semilogy(distance, BER(:,3), '-.d','LineWidth', 3, 'MarkerSize',7, 'Color','g');
semilogy(distance, BER(:,4), ':^','LineWidth', 3, 'MarkerSize',7, 'Color','k');

%yline(1e-2,'r--','LineWidth',2);
%text(100, 2e-4, 'Threshold', 'FontSize',12,'FontWeight','bold')

xlabel('Distance (m)','FontWeight','bold');
ylabel('BER','FontWeight','bold');
%title('BER vs Distance under Moderate Turbulence','FontWeight','bold');
legend('Reference Model [15]','Proposed Model with LDPC','Proposed Model with Gray Coding','Proposed Model with Gray+LDPC','Location','southwest');
ylim([1e-5 1e-1]);
set(gca,'YScale','log','FontWeight','bold','FontSize',12,'LineWidth',1.5);
set(gcf,'Position',[100 100 800 600]);

% High quality export
print(gcf,'BER_vs_Distance_Turbulence_Final','-dpng','-r600');
