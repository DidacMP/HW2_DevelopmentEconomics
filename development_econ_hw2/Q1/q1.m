%%%%%%%%%%%%% DEVELOPMENT ECONOMICS - HOMEWORK 2: QUESTION 1 %%%%%%%%%%%%%%
% Author: Didac Marti Pinto (CEMFI)
% Date: February 8th, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
rng(1234);

%% SET-UP %%

% Set parameters
N = 1000;
T = 12 * 40; % Total periods = 12 months * 40 years
beta = 0.99^(1/12); % Since anual beta (b^12) is 0.99
sigma_sq_e = 0.2;
sigma_e = sqrt(sigma_sq_e);
sigma_sq_u = 0.2; 
sigma_u = sqrt(sigma_sq_u);
eta_1 = 1;
eta_2 = 2;
eta_4 = 4;

% Deterministic seasonal component
gm_low = [-0.073, -0.185, 0.071, 0.066, 0.045, 0.029, 0.018, 0.018, ... 
    0.018, 0.001, -0.017, -0.041];
gm_mid = [-0.147, -0.370, 0.141, 0.131, 0.090, 0.058, 0.036, 0.036, ...
    0.036, 0.002, -0.033, -0.082];
gm_high = [-0.293, -0.739, 0.282, 0.262, 0.180, 0.116, 0.072, 0.072, ...
    0.072, 0.004, -0.066, -0.164];

% Stochastic seasonal component
sigma_sq_m_low = [ 0.043, 0.034, 0.145, 0.142, 0.137, 0.137, 0.119, ...
    0.102, 0.094, 0.094, 0.085, 0.068];
sigma_sq_m_mid = [0.085, 0.068, 0.290, 0.283, 0.273, 0.273, 0.239, ...
    0.205, 0.188, 0.188, 0.171, 0.137];
sigma_sq_m_high = [0.171, 0.137, 0.580, 0.567, 0.546, 0.546, 0.478, ...
    0.410, 0.376, 0.376, 0.341, 0.273];


%% PART 1: SET UP %%

% Create the individual component z_i 
ln_u = mvnrnd(zeros(N,1),eye(N) * sigma_sq_u).'; % Nx1 matrix with the individual ln_u
z = exp(-sigma_sq_u/2) * exp(ln_u); % Nx1 column with individual components Z_i's
Z = z * ones(1,T); % NxT matrix, where each row contains the individual component 

% Create NxT matrices with seasonal components
S_low = exp( kron(ones(N,40),gm_low) );
S_mid = exp( kron(ones(N,40),gm_mid) );
S_high = exp( kron(ones(N,40),gm_high) );

% Create NxT matrix with individual errors for any period
ln_e = zeros(N,T);
for i = 1:N
    for j = 0:39
        ln_e(i,(1+12*j):((j+1)*12)) = normrnd(0,sigma_e);
    end
end
E = exp(-sigma_sq_e/2) * exp(ln_e);

% Calculate matrix of individual consumptions (NxT each matrix)
Cl_s_or = Z .* S_low .* E; % Individual shock X Low seasonal component (common) 
                           % X Other individual and yearly-varying shock
Cm_s_or = Z .* S_mid .* E;
Ch_s_or = Z .* S_high .* E;
Cl_s = Z .* S_low; % Consumption without individual time-varying shock
Cm_s = Z .* S_mid;
Ch_s = Z .* S_high;
C_or = Z .* E; % Consumption without seasonal component
C = Z ; % Consumption without seasonal component and without non-seasonal
        % risk. Just individual shock

% Discounting matrix
beta_m = zeros(1,12);
beta_y = zeros(1,40);
 for i = 1:12
     beta_m(1,i) = beta.^(i-1);
 end
 for i = 1:40 
     beta_y(1,i) = beta.^(12*i);
 end
Betas = ones(N,1) * kron(beta_y,beta_m);


%% PART 1.A. : WELFARE GAINS REMOVING SEASONALITY %%

% Welfare gains removing seasonal component (eta = 1)
gl1_s = zeros(N,1); % Here I store individuals "g"
gm1_s = zeros(N,1);
gh1_s = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* log( Cl_s_or(i,:)*(1+gl))) .' - ...
          (Betas(i,:) .* log(C_or(i,:))) .').' );
gl1_s(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* log( Cm_s_or(i,:)*(1+gm))) .' - ...
          (Betas(i,:) .* log(C_or(i,:))) .').' );
gm1_s(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* log( Ch_s_or(i,:)*(1+gh))) .' - ...
          (Betas(i,:) .* log(C_or(i,:))) .').' );
gh1_s(i,1) = fminbnd(funh,-2,5);
end

    % Results:
Results = [mean(gl1_s), mean(gm1_s), mean(gh1_s)];
disp(' RESULTS PART 1.A. Mean welfare gains of removing seasonality (eta=1)')
disp(Results)
disp('Each column: Mean "g" for low, mid, and high seasonality')
disp(' ')
disp(' ')

%% PART 1.B. : WELFARE GAINS REMOVING NON-SEASONAL RISK (EPSILONS) %%

% Welfare gains removing non-seasonal consumption risk (eta = 1)
gl1_or = zeros(N,1);
gm1_or = zeros(N,1);
gh1_or = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* log( Cl_s_or(i,:)*(1+gl))) .' - ...
          (Betas(i,:) .* log(Cl_s(i,:))) .').' );
gl1_or(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* log( Cm_s_or(i,:)*(1+gm))) .' - ...
          (Betas(i,:) .* log(Cm_s(i,:))) .').' );
gm1_or(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* log( Ch_s_or(i,:)*(1+gh))) .' - ...
          (Betas(i,:) .* log(Ch_s(i,:))) .').' );
gh1_or(i,1) = fminbnd(funh,-2,5);
end

% Alternative: shut down seasonal risk, and measure the gains then (eta = 1)
g1_or_alt = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* log( C_or(i,:)*(1+gl))) .' - ...
          (Betas(i,:) .* log(C(i,:))) .').' );
g1_or_alt(i,1) = fminbnd(funl,-2,5);
end

Results = [prctile(gl1_or,10), prctile(gm1_or,10), prctile(gh1_or,10), prctile(g1_or_alt,10); ...
    prctile(gl1_or,50), prctile(gm1_or,50), prctile(gh1_or,50), prctile(g1_or_alt,50); ...
    prctile(gl1_or,90), prctile(gm1_or,90), prctile(gh1_or,90), prctile(g1_or_alt,90); ...
    mean(gl1_or), mean(gm1_or), mean(gh1_or), mean(g1_or_alt); ...
    std(gl1_or), std(gm1_or), std(gh1_or), std(g1_or_alt)];

% Results
disp(' RESULTS PART 1.B. Welfare gains of removing non-seasonal consumption risk (eta=1)')
disp(Results)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid, high, and no seasonality')
disp(' ')
disp(' ')

% Graph to compare distribution of welfare gains of removing non-seasonal consumption risk
hold on
hist(g1_or_alt);
xlabel('Individual g')
ylabel('Number of households')
legend('eta=1')
title({'Q1 Part1 B', 'Welfare gains of removing non-seasonal consumption risk'})
print('Q1_Part1_B','-dpng')


%% PART 1.D. : RE-DO FOR ETA = 2 & ETA = 4 %%

% Welfare gains removing seasonal component (eta = 2)
gl2_s = zeros(N,1);
gm2_s = zeros(N,1);
gh2_s = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* ((Cl_s_or(i,:)*(1+gl))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((C_or(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gl2_s(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* ((Cm_s_or(i,:)*(1+gm))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((C_or(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gm2_s(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* ((Ch_s_or(i,:)*(1+gh))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((C_or(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gh2_s(i,1) = fminbnd(funh,-2,5);
end

% Welfare gains removing seasonal component (eta = 4)
gl4_s = zeros(N,1);
gm4_s = zeros(N,1);
gh4_s = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* ((Cl_s_or(i,:)*(1+gl))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((C_or(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gl4_s(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* ((Cm_s_or(i,:)*(1+gm))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((C_or(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gm4_s(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* ((Ch_s_or(i,:)*(1+gh))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((C_or(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gh4_s(i,1) = fminbnd(funh,-2,5);
end

% Results:
Results = [mean(gl2_s), mean(gm2_s), mean(gh2_s); ...
    mean(gl4_s), mean(gm4_s), mean(gh4_s)];
disp(' RESULTS PART 1.D. Mean welfare gains of removing seasonality (eta=2 and 4)')
disp(Results)
disp('Row 1: eta=2. Row 2: eta=4')
disp('Each column: Mean "g" for low, mid, and high seasonality')
disp(' ')
disp(' ')

% Welfare gains removing non-seasonal consumption risk (eta = 2)
gl2_or = zeros(N,1);
gm2_or = zeros(N,1);
gh2_or = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* ((Cl_s_or(i,:)*(1+gl))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((Cl_s(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gl2_or(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* ((Cm_s_or(i,:)*(1+gm))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((Cm_s(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gm2_or(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* ((Ch_s_or(i,:)*(1+gh))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((Ch_s(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gh2_or(i,1) = fminbnd(funh,-2,5);
end

% Alternative: shut down seasonal risk, and measure the gains then (eta = 2)
g2_or_alt = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* ((C_or(i,:)*(1+gl))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((C(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
g2_or_alt(i,1) = fminbnd(funl,-2,5);
end

% Welfare gains removing non-seasonal consumption risk (eta = 4)
gl4_or = zeros(N,1);
gm4_or = zeros(N,1);
gh4_or = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* ((Cl_s_or(i,:)*(1+gl))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((Cl_s(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gl4_or(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* ((Cm_s_or(i,:)*(1+gm))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((Cm_s(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gm4_or(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* ((Ch_s_or(i,:)*(1+gh))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((Ch_s(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gh4_or(i,1) = fminbnd(funh,-2,5);
end

% Alternative: shut down seasonal risk, and measure the gains then (eta = 4)
g4_or_alt = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* ((C_or(i,:)*(1+gl))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((C(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
g4_or_alt(i,1) = fminbnd(funl,-2,5);
end

Results2 = [prctile(gl2_or,10), prctile(gm2_or,10), prctile(gh2_or,10); ...
    prctile(gl2_or,50), prctile(gm2_or,50), prctile(gh2_or,50); ...
    prctile(gl2_or,90), prctile(gm2_or,90), prctile(gh2_or,90); ...
    mean(gl2_or), mean(gm2_or), mean(gh2_or); ...
    std(gl2_or), std(gm2_or), std(gh2_or)];

Results4 = [prctile(gl4_or,10), prctile(gm4_or,10), prctile(gh4_or,10); ...
    prctile(gl4_or,50), prctile(gm4_or,50), prctile(gh4_or,50); ...
    prctile(gl4_or,90), prctile(gm4_or,90), prctile(gh4_or,90); ...
    mean(gl4_or), mean(gm4_or), mean(gh4_or); ...
    std(gl4_or), std(gm4_or), std(gh4_or)];

% Results (eta=2)
disp(' RESULTS PART 1.D. Welfare gains of removing non-seasonal consumption risk (eta=2)')
disp(Results2)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid and high seasonality')
disp(' ')
disp(' ')

% Results (eta=4)
disp(' RESULTS PART 1.D. Welfare gains of removing non-seasonal consumption risk (eta=4)')
disp(Results4)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid and high seasonality')
disp(' ')
disp(' ')

% Graph to compare distribution of welfare gains of removing non-seasonal consumption risk
figure 
hold on
histogram(g1_or_alt);
hold on
histogram(g2_or_alt);
hold on
histogram(g4_or_alt);
xlabel('Individual g')
ylabel('Number of households')
legend('eta=1','eta=2','eta=4')
title({'Q1 Part1 D', 'Welfare gains of removing non-seasonal consumption risk'})
print('Q1_Part1_D','-dpng')


%% PART 2: SET-UP %%


% Create a NxT matrix with seasonal components
Sr_low = zeros(N,T);
for i = 1:N
    for j = 0:39
      Sr_low(i,(1+j*12):(j+1)*12) = exp(-sigma_sq_m_low/2) .* ...
      exp( mvnrnd(zeros(12,1), ones(12,1) * sigma_sq_m_low .* eye(12) ) );
      % Each one is a row of 12 seasonal components used for each indiv.
    end
end
Sr_mid = zeros(N,T);
for i = 1:N
    for j = 0:39
      Sr_mid(i,(1+j*12):(j+1)*12) = exp(-sigma_sq_m_mid/2) .* ...
      exp( mvnrnd(zeros(12,1), ones(12,1) * sigma_sq_m_mid .* eye(12) ) );
      % Each one is a row of 12 seasonal components used for each indiv.
    end
end
Sr_high = zeros(N,T);
for i = 1:N
    for j = 0:39
      Sr_high(i,(1+j*12):(j+1)*12) = exp(-sigma_sq_m_high/2) .* ...
      exp( mvnrnd(zeros(12,1), ones(12,1) * sigma_sq_m_high .* eye(12) ) );
      % Each one is a row of 12 seasonal components used for each indiv.
    end
end

% Calculate matrix of individual consumptions (NxT each matrix)
C_sm_srl_or = Z .* S_mid .*  Sr_low .* E; 
C_sm_srm_or = Z .* S_mid .* Sr_mid .* E;
C_sm_srh_or = Z .* S_mid .* Sr_high .* E;

C_srl_or = Z .* Sr_low .* E; 
C_srm_or = Z .* Sr_mid .* E;
C_srh_or = Z .* Sr_high .* E;

C_sm_srl = Z .* S_mid .*  Sr_low; 
C_sm_srm = Z .* S_mid .* Sr_mid;
C_sm_srh = Z .* S_mid .* Sr_high;

C_sm_or = Z .* S_mid .* E; 


%% PART 2.A : WELFARE GAINS REMOVING SEASONALITY %%

% Welfare gains removing deterministic seasonal component (eta = 1)
gl1_s2 = zeros(N,1); 
gm1_s2 = zeros(N,1);
gh1_s2 = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* log(C_sm_srl_or(i,:)*(1+gl))) .' - ...
          (Betas(i,:) .* log(C_srl_or(i,:))) .').' );
gl1_s2(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* log(C_sm_srm_or(i,:)*(1+gm))) .' - ...
          (Betas(i,:) .* log(C_srm_or(i,:))) .').' );
gm1_s2(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* log(C_sm_srh_or(i,:)*(1+gh))) .' - ...
          (Betas(i,:) .* log(C_srh_or(i,:))) .').' );
gh1_s2(i,1) = fminbnd(funh,-2,5);
end

% Welfare gains removing stochastic seasonal component (eta = 1)
gl1_sr = zeros(N,1); 
gm1_sr = zeros(N,1);
gh1_sr = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* log(C_sm_srl_or(i,:)*(1+gl))) .' - ...
          (Betas(i,:) .* log(C_sm_or(i,:))) .').' );
gl1_sr(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* log(C_sm_srm_or(i,:)*(1+gm))) .' - ...
          (Betas(i,:) .* log(C_sm_or(i,:))) .').' );
gm1_sr(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* log(C_sm_srh_or(i,:)*(1+gh))) .' - ...
          (Betas(i,:) .* log(C_sm_or(i,:))) .').' );
gh1_sr(i,1) = fminbnd(funh,-2,5);
end
  
% Welfare gains removing both seasonal components (eta = 1)
gl1_s_sr = zeros(N,1); 
gm1_s_sr = zeros(N,1);
gh1_s_sr = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* log(C_sm_srl_or(i,:)*(1+gl))) .' - ...
          (Betas(i,:) .* log(C_or(i,:))) .').' );
gl1_s_sr(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* log(C_sm_srm_or(i,:)*(1+gm))) .' - ...
          (Betas(i,:) .* log(C_or(i,:))) .').' );
gm1_s_sr(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* log(C_sm_srh_or(i,:)*(1+gh))) .' - ...
          (Betas(i,:) .* log(C_or(i,:))) .').' );
gh1_s_sr(i,1) = fminbnd(funh,-2,5);
end

% Results removing deterministic component:
Results = [mean(gl1_s2), mean(gm1_s2), mean(gh1_s2)];
disp(' RESULTS PART 2.A. Welfare gains of removing deterministic seasonal component (eta=1)')
disp(Results)
disp('Rows: Means')
disp('Each column: "g" for low, mid and high seasonality')
disp(' ')
disp(' ')

% Results removing stochastic component:
Results = [prctile(gl1_sr,10), prctile(gm1_sr,10), prctile(gh1_sr,10); ...
    prctile(gl1_sr,50), prctile(gm1_sr,50), prctile(gh1_sr,50); ...
    prctile(gl1_sr,90), prctile(gm1_sr,90), prctile(gh1_sr,90); ...
    mean(gl1_sr), mean(gm1_sr), mean(gh1_sr); ...
    std(gl1_sr), std(gm1_sr), std(gh1_sr)];
disp(' RESULTS PART 2.A. Welfare gains of removing stochastic seasonal components (eta=1)')
disp(Results)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid and high seasonality')
disp(' ')
disp(' ')

% Results removing both stochastic and deterministic component:
Results = [prctile(gl1_s_sr,10), prctile(gm1_s_sr,10), prctile(gh1_s_sr,10); ...
    prctile(gl1_s_sr,50), prctile(gm1_s_sr,50), prctile(gh1_s_sr,50); ...
    prctile(gl1_s_sr,90), prctile(gm1_s_sr,90), prctile(gh1_s_sr,90); ...
    mean(gl1_s_sr), mean(gm1_s_sr), mean(gh1_s_sr); ...
    std(gl1_s_sr), std(gm1_s_sr), std(gh1_s_sr)];
disp(' RESULTS PART 2.A. Welfare gains of removing both seasonality components (eta=1)')
disp(Results)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid and high seasonality')
disp(' ')
disp(' ')

% Graph of the distribution of welfare gains of removing the seasonal stochastic component
figure 
hold on
histogram(gl1_sr);
hold on
histogram(gm1_sr);
hold on
histogram(gh1_sr);
xlabel('Individual g')
ylabel('Number of households')
legend('Low seasonality dispersion','Medium seasonality dispersion','High seasonality dispersion')
title({'Q1 Part2 A1', 'Welfare gains of removing stochastic seasonal component'})
print('Q1_Part2_A1','-dpng')

% Graph of the distribution of welfare gains of removing both seasonal components
figure 
hold on
histogram(gl1_s_sr);
hold on
histogram(gm1_s_sr);
hold on
histogram(gh1_s_sr);
xlabel('Individual g')
ylabel('Number of households')
legend('Low seasonality dispersion','Medium seasonality dispersion','High seasonality dispersion')
title({'Q1 Part2 A2', 'Welfare gains of removing both seasonal components'})
print('Q1_Part2_A2','-dpng')

  
%% PART 2.B. : WELFARE GAINS REMOVING NON-SEASONAL RISK  %%

% Welfare gains removing non-seasonal consumption risk (eta = 1)
gl1_or_q2 = zeros(N,1);
gm1_or_q2 = zeros(N,1);
gh1_or_q2 = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* log( C_sm_srl_or(i,:)*(1+gl))) .' - ...
          (Betas(i,:) .* log(C_sm_srl(i,:))) .').' );
gl1_or_q2(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* log(C_sm_srm_or(i,:)*(1+gm))) .' - ...
          (Betas(i,:) .* log(C_sm_srm(i,:))) .').' );
gm1_or_q2(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* log(C_sm_srh_or(i,:)*(1+gh))) .' - ...
          (Betas(i,:) .* log(C_sm_srh(i,:))) .').' );
gh1_or_q2(i,1) = fminbnd(funh,-2,5);
end

% Results
Results = [prctile(gl1_or_q2,10), prctile(gm1_or_q2,10), prctile(gh1_or_q2,10); ...
    prctile(gl1_or_q2,50), prctile(gm1_or_q2,50), prctile(gh1_or_q2,50); ...
    prctile(gl1_or_q2,90), prctile(gm1_or_q2,90), prctile(gh1_or_q2,90); ...
    mean(gl1_or_q2), mean(gm1_or_q2), mean(gh1_or_q2); ...
    std(gl1_or_q2), std(gm1_or_q2), std(gh1_or_q2)];

disp(' RESULTS PART 2.B. Welfare gains of removing non-seasonal consumption risk (eta=1)')
disp(Results)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid and high seasonality')
disp(' ')
disp(' ')

% Graph to compare distribution of welfare gains of removing non-seasonal consumption risk
figure
hold on
hist(gl1_or_q2);
xlabel('Individual g')
ylabel('Number of households')
legend('eta=1')
title({'Q1 Part2 B', 'Welfare gains of removing non-seasonal consumption risk'})
print('Q1_Part2_B','-dpng')



%% PART 2.D. : RE-DO FOR ETA = 2 & ETA = 4 %%

% Welfare gains removing deterministic seasonal component (eta = 2)
gl2_s2 = zeros(N,1); 
gm2_s2 = zeros(N,1);
gh2_s2 = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* ((C_sm_srl_or(i,:)*(1+gl))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((C_srl_or(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gl2_s2(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* ((C_sm_srm_or(i,:)*(1+gm))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((C_srm_or(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gm2_s2(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* ((C_sm_srh_or(i,:)*(1+gh))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((C_srh_or(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gh2_s2(i,1) = fminbnd(funh,-2,5);
end


% Welfare gains removing deterministic seasonal component (eta = 4)
gl4_s2 = zeros(N,1); 
gm4_s2 = zeros(N,1);
gh4_s2 = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* ((C_sm_srl_or(i,:)*(1+gl))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((C_srl_or(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gl4_s2(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* ((C_sm_srm_or(i,:)*(1+gm))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((C_srm_or(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gm4_s2(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* ((C_sm_srh_or(i,:)*(1+gh))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((C_srh_or(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gh4_s2(i,1) = fminbnd(funh,-2,5);
end

% Results removing deterministic component (eta = 2)
Results = [mean(gl2_s2), mean(gm2_s2), mean(gh2_s2)];
disp(' RESULTS PART 2.D. Welfare gains of removing deterministic seasonal component (eta=2)')
disp(Results)
disp('Rows: Means')
disp('Each column: "g" for low, mid and high seasonality')
disp(' ')
disp(' ')

% Results removing deterministic component (eta = 4)
Results = [mean(gl4_s2), mean(gm4_s2), mean(gh4_s2)];
disp(' RESULTS PART 2.D. Welfare gains of removing deterministic seasonal component (eta=4)')
disp(Results)
disp('Rows: Means')
disp('Each column: "g" for low, mid and high seasonality')
disp(' ')
disp(' ')



% Welfare gains removing stochastic seasonal component (eta = 2)
gl2_sr = zeros(N,1); 
gm2_sr = zeros(N,1);
gh2_sr = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* ((C_sm_srl_or(i,:)*(1+gl))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((C_sm_or(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gl2_sr(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* ((C_sm_srm_or(i,:)*(1+gm))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((C_sm_or(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gm2_sr(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* ((C_sm_srh_or(i,:)*(1+gh))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((C_sm_or(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gh2_sr(i,1) = fminbnd(funh,-2,5);
end


% Welfare gains removing stochastic seasonal component (eta = 4)
gl4_sr = zeros(N,1); 
gm4_sr = zeros(N,1);
gh4_sr = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* ((C_sm_srl_or(i,:)*(1+gl))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((C_sm_or(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gl4_sr(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* ((C_sm_srm_or(i,:)*(1+gm))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((C_sm_or(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gm4_sr(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* ((C_sm_srh_or(i,:)*(1+gh))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((C_sm_or(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gh4_sr(i,1) = fminbnd(funh,-2,5);
end

% Results removing stochastic component (eta=2)
Results = [prctile(gl2_sr,10), prctile(gm2_sr,10), prctile(gh2_sr,10); ...
    prctile(gl2_sr,50), prctile(gm2_sr,50), prctile(gh2_sr,50); ...
    prctile(gl2_sr,90), prctile(gm2_sr,90), prctile(gh2_sr,90); ...
    mean(gl2_sr), mean(gm2_sr), mean(gh2_sr); ...
    std(gl2_sr), std(gm2_sr), std(gh2_sr)];
disp(' RESULTS PART 2.D. Welfare gains of removing stochastic seasonal components (eta=2)')
disp(Results)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid and high seasonality')
disp(' ')
disp(' ')

% Results removing stochastic component (eta=4)
Results = [prctile(gl4_sr,10), prctile(gm4_sr,10), prctile(gh4_sr,10); ...
    prctile(gl4_sr,50), prctile(gm4_sr,50), prctile(gh4_sr,50); ...
    prctile(gl4_sr,90), prctile(gm4_sr,90), prctile(gh4_sr,90); ...
    mean(gl4_sr), mean(gm4_sr), mean(gh4_sr); ...
    std(gl4_sr), std(gm4_sr), std(gh4_sr)];
disp(' RESULTS PART 2.D. Welfare gains of removing stochastic seasonal components (eta=4)')
disp(Results)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid and high seasonality')
disp(' ')
disp(' ')
  






% Welfare gains removing both seasonal components (eta = 2)
gl2_s_sr = zeros(N,1); 
gm2_s_sr = zeros(N,1);
gh2_s_sr = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* ((C_sm_srl_or(i,:)*(1+gl))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((C_or(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gl2_s_sr(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* ((C_sm_srm_or(i,:)*(1+gm))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((C_or(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gm2_s_sr(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* ((C_sm_srh_or(i,:)*(1+gh))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((C_or(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gh2_s_sr(i,1) = fminbnd(funh,-2,5);
end

% Welfare gains removing both seasonal components (eta = 4)
gl4_s_sr = zeros(N,1); 
gm4_s_sr = zeros(N,1);
gh4_s_sr = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* ((C_sm_srl_or(i,:)*(1+gl))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((C_or(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gl4_s_sr(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* ((C_sm_srm_or(i,:)*(1+gm))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((C_or(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gm4_s_sr(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* ((C_sm_srh_or(i,:)*(1+gh))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((C_or(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gh4_s_sr(i,1) = fminbnd(funh,-2,5);
end

% Results removing both stochastic and deterministic component (eta=2):
Results = [prctile(gl2_s_sr,10), prctile(gm2_s_sr,10), prctile(gh2_s_sr,10); ...
    prctile(gl2_s_sr,50), prctile(gm2_s_sr,50), prctile(gh2_s_sr,50); ...
    prctile(gl2_s_sr,90), prctile(gm2_s_sr,90), prctile(gh2_s_sr,90); ...
    mean(gl2_s_sr), mean(gm2_s_sr), mean(gh2_s_sr); ...
    std(gl2_s_sr), std(gm2_s_sr), std(gh2_s_sr)];
disp(' RESULTS PART 2.D. Welfare gains of removing both seasonality components (eta=2)')
disp(Results)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid and high seasonality')
disp(' ')
disp(' ')

% Results removing both stochastic and deterministic component (eta=4):
Results = [prctile(gl4_s_sr,10), prctile(gm4_s_sr,10), prctile(gh4_s_sr,10); ...
    prctile(gl4_s_sr,50), prctile(gm4_s_sr,50), prctile(gh4_s_sr,50); ...
    prctile(gl4_s_sr,90), prctile(gm4_s_sr,90), prctile(gh4_s_sr,90); ...
    mean(gl4_s_sr), mean(gm4_s_sr), mean(gh4_s_sr); ...
    std(gl4_s_sr), std(gm4_s_sr), std(gh4_s_sr)];
disp(' RESULTS PART 2.D. Welfare gains of removing both seasonality components (eta=4)')
disp(Results)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid and high seasonality')
disp(' ')
disp(' ')





% Welfare gains removing non-seasonal consumption risk (eta = 2)
gl2_or_q2 = zeros(N,1);
gm2_or_q2 = zeros(N,1);
gh2_or_q2 = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* ((C_sm_srl_or(i,:)*(1+gl))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((C_sm_srl(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gl2_or_q2(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* ((C_sm_srm_or(i,:)*(1+gm))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((C_sm_srm(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gm2_or_q2(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* ((C_sm_srh_or(i,:)*(1+gh))).^(1-eta_2) / (1-eta_2)).' - ...
          (Betas(i,:) .* ((C_sm_srh(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
gh2_or_q2(i,1) = fminbnd(funh,-2,5);
end


% Welfare gains removing non-seasonal consumption risk (eta = 4)
gl4_or_q2 = zeros(N,1);
gm4_or_q2 = zeros(N,1);
gh4_or_q2 = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum( (Betas(i,:) .* ((C_sm_srl_or(i,:)*(1+gl))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((C_sm_srl(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gl4_or_q2(i,1) = fminbnd(funl,-2,5);

funm = @(gm) abs(sum( (Betas(i,:) .* ((C_sm_srm_or(i,:)*(1+gm))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((C_sm_srm(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gm4_or_q2(i,1) = fminbnd(funm,-2,5);

funh = @(gh) abs(sum( (Betas(i,:) .* ((C_sm_srh_or(i,:)*(1+gh))).^(1-eta_4) / (1-eta_4)).' - ...
          (Betas(i,:) .* ((C_sm_srh(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
gh4_or_q2(i,1) = fminbnd(funh,-2,5);
end



% Results
Results = [prctile(gl2_or_q2,10), prctile(gm2_or_q2,10), prctile(gh2_or_q2,10); ...
    prctile(gl2_or_q2,50), prctile(gm2_or_q2,50), prctile(gh2_or_q2,50); ...
    prctile(gl2_or_q2,90), prctile(gm2_or_q2,90), prctile(gh2_or_q2,90); ...
    mean(gl2_or_q2), mean(gm2_or_q2), mean(gh2_or_q2); ...
    std(gl2_or_q2), std(gm2_or_q2), std(gh2_or_q2)];

disp(' RESULTS PART 2.D. Welfare gains of removing non-seasonal consumption risk (eta=2)')
disp(Results)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid and high seasonality')
disp(' ')
disp(' ')

% Results
Results = [prctile(gl4_or_q2,10), prctile(gm4_or_q2,10), prctile(gh4_or_q2,10); ...
    prctile(gl4_or_q2,50), prctile(gm4_or_q2,50), prctile(gh4_or_q2,50); ...
    prctile(gl4_or_q2,90), prctile(gm4_or_q2,90), prctile(gh4_or_q2,90); ...
    mean(gl4_or_q2), mean(gm4_or_q2), mean(gh4_or_q2); ...
    std(gl4_or_q2), std(gm4_or_q2), std(gh4_or_q2)];

disp(' RESULTS PART 2.D. Welfare gains of removing non-seasonal consumption risk (eta=4)')
disp(Results)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid and high seasonality')
disp(' ')
disp(' ')


%%% THE END %%%