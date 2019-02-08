%%%%%%%%%%%%% DEVELOPMENT ECONOMICS - HOMEWORK 2: QUESTION 2 %%%%%%%%%%%%%%
% Author: Didac Marti Pinto (CEMFI)
% Date: February 8th, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% kappa wtf is wrong with it?
% Utility function something wrong. it is increasing and convex in hours worked!
%    so far I assumed a negative 

clc; clear; close all;
rng(1234);

%% SET-UP: GENERAL %%

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
nu = 1;

% Calibration for kappa
theta = 0.6; % Labor share of total output
y_over_c = 1/0.25; % Output over consumption
h_month = 28.5 * (30/7); % Hours worked per month (Bick et al 2018)
kappa = theta * y_over_c * (h_month)^(-1-1/nu);


% Deterministic seasonal component
gm_low = [-0.073, -0.185, 0.071, 0.066, 0.045, 0.029, 0.018, 0.018, ... 
    0.018, 0.001, -0.017, -0.041];
gm_mid = [-0.147, -0.370, 0.141, 0.131, 0.090, 0.058, 0.036, 0.036, ...
    0.036, 0.002, -0.033, -0.082];
gm_high = [-0.293, -0.739, 0.282, 0.262, 0.180, 0.116, 0.072, 0.072, ...
    0.072, 0.004, -0.066, -0.164];

% Deterministic seasonal component in negative
gm_low_negative = [+0.073, +0.185, -0.071, -0.066, -0.045, -0.029, -0.018, -0.018, ... 
    -0.018, -0.001, +0.017, +0.041];
gm_mid_negative = [+0.147, +0.370, -0.141, -0.131, -0.090, -0.058, -0.036, -0.036, ...
    -0.036, -0.002, +0.033, +0.082];
gm_high_negative = [+0.293, +0.739, -0.282, -0.262, -0.180, -0.116, -0.072, -0.072, ...
    -0.072, -0.004, +0.066, +0.164];

% Stochastic seasonal component
sigma_sq_m_low = [ 0.043, 0.034, 0.145, 0.142, 0.137, 0.137, 0.119, ...
    0.102, 0.094, 0.094, 0.085, 0.068];
sigma_sq_m_mid = [0.085, 0.068, 0.290, 0.283, 0.273, 0.273, 0.239, ...
    0.205, 0.188, 0.188, 0.171, 0.137];
sigma_sq_m_high = [0.171, 0.137, 0.580, 0.567, 0.546, 0.546, 0.478, ...
    0.410, 0.376, 0.376, 0.341, 0.273];

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

%% SET UP: CONSUMPTION AND LABOR COMMON %%

% Create NxT matrices with deterministic seasonal components
S_low = exp( kron(ones(N,40),gm_low) );
S_mid = exp( kron(ones(N,40),gm_mid) );
S_high = exp( kron(ones(N,40),gm_high) );
S_low_negative = exp( kron(ones(N,40),gm_low_negative) );
S_mid_negative = exp( kron(ones(N,40),gm_mid_negative) );
S_high_negative = exp( kron(ones(N,40),gm_high_negative) );

% Positively correlated: Stochastic seasonal components for both consumption and labor
Sr_low = zeros(N,T);
lab_Sr_low = zeros(N,T);
for k = 1:1000
    for j = 0:39
        for i = 1:12
        varcov = [sigma_sq_m_low(1,i), 0.03; 0.03,sigma_sq_m_low(1,i)];
        ln = mvnrnd(zeros(2,1), varcov);
        Sr_low(k,i+j*12) = exp(-sigma_sq_m_low(1,i)/2) * exp(ln(1,1));
        lab_Sr_low(k,i+j*12) = exp(-sigma_sq_m_low(1,i)/2) * exp(ln(1,2));
        end
    end
end
Sr_mid = zeros(N,T);
lab_Sr_mid = zeros(N,T);
for k = 1:1000
    for j = 0:39
        for i = 1:12
        varcov = [sigma_sq_m_mid(1,i), 0.03; 0.03,sigma_sq_m_mid(1,i)];
        ln = mvnrnd(zeros(2,1), varcov);
        Sr_mid(k,i+j*12) = exp(-sigma_sq_m_mid(1,i)/2) * exp(ln(1,1));
        lab_Sr_mid(k,i+j*12) = exp(-sigma_sq_m_mid(1,i)/2) * exp(ln(1,2));
        end
    end
end
Sr_high = zeros(N,T);
lab_Sr_high = zeros(N,T);
for k = 1:1000
    for j = 0:39
        for i = 1:12
        varcov = [sigma_sq_m_high(1,i), 0.03; 0.03,sigma_sq_m_high(1,i)];
        ln = mvnrnd(zeros(2,1), varcov);
        Sr_high(k,i+j*12) = exp(-sigma_sq_m_high(1,i)/2) * exp(ln(1,1));
        lab_Sr_high(k,i+j*12) = exp(-sigma_sq_m_high(1,i)/2) * exp(ln(1,2));
        end
    end
end

% Negatively correlated: Stochastic seasonal components for both consumption and labor
Sr_low_negative = zeros(N,T);
lab_Sr_low_negative = zeros(N,T);
for k = 1:1000
    for j = 0:39
        for i = 1:12
        varcov = [sigma_sq_m_low(1,i), -0.03; -0.03,sigma_sq_m_low(1,i)];
        ln = mvnrnd(zeros(2,1), varcov);
        Sr_low_negative(k,i+j*12) = exp(-sigma_sq_m_low(1,i)/2) * exp(ln(1,1));
        lab_Sr_low_negative(k,i+j*12) = exp(-sigma_sq_m_low(1,i)/2) * exp(ln(1,2));
        end
    end
end
Sr_mid_negative = zeros(N,T);
lab_Sr_mid_negative = zeros(N,T);
for k = 1:1000
    for j = 0:39
        for i = 1:12
        varcov = [sigma_sq_m_mid(1,i), -0.03; -0.03,sigma_sq_m_mid(1,i)];
        ln = mvnrnd(zeros(2,1), varcov);
        Sr_mid_negative(k,i+j*12) = exp(-sigma_sq_m_mid(1,i)/2) * exp(ln(1,1));
        lab_Sr_mid_negative(k,i+j*12) = exp(-sigma_sq_m_mid(1,i)/2) * exp(ln(1,2));
        end
    end
end
Sr_high_negative = zeros(N,T);
lab_Sr_high_negative = zeros(N,T);
for k = 1:1000
    for j = 0:39
        for i = 1:12
        varcov = [sigma_sq_m_high(1,i), -0.03; -0.03,sigma_sq_m_high(1,i)];
        ln = mvnrnd(zeros(2,1), varcov);
        Sr_high_negative(k,i+j*12) = exp(-sigma_sq_m_high(1,i)/2) * exp(ln(1,1));
        lab_Sr_high_negative(k,i+j*12) = exp(-sigma_sq_m_high(1,i)/2) * exp(ln(1,2));
        end
    end
end

%% SET UP: CONSUMPTION %%

% Create the individual component z_i 
ln_u = mvnrnd(zeros(N,1),eye(N) * sigma_sq_u).'; 
z = exp(-sigma_sq_u/2) * exp(ln_u); 
Z = z * ones(1,T); 

% Create NxT matrix with non-seasonal shocks for any period
ln_e = zeros(N,T);
for i = 1:N
    for j = 0:39
        ln_e(i,(1+12*j):((j+1)*12)) = normrnd(0,sigma_e);
    end
end
E = exp(-sigma_sq_e/2) * exp(ln_e);


%% SET UP: LABOR %%

% Create the individual component lab_z_i 
lab_ln_u = mvnrnd(zeros(N,1),eye(N) * sigma_sq_u).'; 
lab_z = exp(-sigma_sq_u/2) * exp(lab_ln_u); 
lab_Z = lab_z * ones(1,T); 

% Create NxT matrix with non-seasonal shocks for any period
lab_ln_e = zeros(N,T);
for i = 1:N
    for j = 0:39
        lab_ln_e(i,(1+12*j):((j+1)*12)) = normrnd(0,sigma_e);
    end
end
lab_E = exp(-sigma_sq_e/2) * exp(lab_ln_e);


%% SET UP: VALUES OF CONSUMPTION AND LABOR POSITIVELY CORRELATED

% Calculate matrix of individual consumptions (NxT each matrix)
Cl_s_sr_or = Z .* S_low .*  Sr_low .* E; 
Cm_s_sr_or = Z .* S_mid .* Sr_mid .* E;
Ch_s_sr_or = Z .* S_high .* Sr_high .* E;
Cl_s_sr = Z .* S_low .*  Sr_low; 
Cm_s_sr = Z .* S_mid .* Sr_mid;
Ch_s_sr = Z .* S_high .* Sr_high;
Cl_s_or = Z .* S_low .* E; 
Cm_s_or = Z .* S_mid .* E;
Ch_s_or = Z .* S_high .* E;
Cl_s = Z .* S_low; 
Cm_s = Z .* S_mid;
Ch_s = Z .* S_high;
C_or = Z .* E; 
C = Z;

% Calculate matrix of individual labors (NxT each matrix)
Ll_s_sr_or = lab_Z .* S_low .*  lab_Sr_low .* lab_E; 
Lm_s_sr_or = lab_Z .* S_mid .* lab_Sr_mid .* lab_E;
Lh_s_sr_or = lab_Z .* S_high .* lab_Sr_high .* lab_E;
Ll_s_sr = lab_Z .* S_low .*  lab_Sr_low; 
Lm_s_sr = lab_Z .* S_mid .* lab_Sr_mid;
Lh_s_sr = lab_Z .* S_high .* lab_Sr_high;
Ll_s_or = lab_Z .* S_low .* lab_E; 
Lm_s_or = lab_Z .* S_mid .* lab_E;
Lh_s_or = lab_Z .* S_high .* lab_E;
Ll_s = lab_Z .* S_low; 
Lm_s = lab_Z .* S_mid;
Lh_s = lab_Z .* S_high;
L_or = lab_Z .* lab_E; 
L = lab_Z;

%% SET UP: VALUES OF CONSUMPTION AND LABOR NEGATIVELY CORRELATED

% Calculate matrix of individual consumptions (NxT each matrix)
Cl_s_sr_or_negative = Z .* S_low .*  Sr_low_negative .* E; 
Cm_s_sr_or_negative = Z .* S_mid .* Sr_mid_negative .* E;
Ch_s_sr_or_negative = Z .* S_high .* Sr_high_negative .* E;
Cl_s_sr_negative = Z .* S_low .*  Sr_low_negative; 
Cm_s_sr_negative = Z .* S_mid .* Sr_mid_negative;
Ch_s_sr_negative = Z .* S_high .* Sr_high_negative;

% Calculate matrix of individual labors (NxT each matrix)
Ll_s_sr_or_negative = lab_Z .* S_low_negative .*  lab_Sr_low_negative .* lab_E; 
Lm_s_sr_or_negative = lab_Z .* S_mid_negative .* lab_Sr_mid_negative .* lab_E;
Lh_s_sr_or_negative = lab_Z .* S_high_negative .* lab_Sr_high_negative .* lab_E;
Ll_s_sr_negative = lab_Z .* S_low_negative .*  lab_Sr_low_negative; 
Lm_s_sr_negative = lab_Z .* S_mid_negative .* lab_Sr_mid_negative;
Lh_s_sr_negative = lab_Z .* S_high_negative .* lab_Sr_high_negative;


%% PART A : HIGHLY POSITIVELY CORRELATED CONSUMPTION AND LEISURE %%

% Total effects
gl_s_sr = zeros(N,1); 
gm_s_sr = zeros(N,1);
gh_s_sr = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum(...
    (Betas(i,:).*(  log(Cl_s_sr_or(i,:).*(1+gl))  -  kappa.*(Ll_s_sr_or(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:))  -  kappa.*(L_or(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gl_s_sr(i,1) = fminbnd(funl,-5,20);

funm = @(gm) abs(sum(...
    (Betas(i,:).*(  log(Cm_s_sr_or(i,:).*(1+gm))  -  kappa.*(Lm_s_sr_or(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:))  -  kappa.*(L_or(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gm_s_sr(i,1) = fminbnd(funm,-5,20);

funh = @(gh) abs(sum(...
    (Betas(i,:).*(  log(Ch_s_sr_or(i,:).*(1+gh))  -  kappa.*(Lh_s_sr_or(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:))  -  kappa.*(L_or(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gh_s_sr(i,1) = fminbnd(funh,-5,20);
end

% Consumption effects
gl_s_sr_conseff = zeros(N,1); 
gm_s_sr_conseff = zeros(N,1);
gh_s_sr_conseff = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum(...
    (Betas(i,:).*(  log(Cl_s_sr_or(i,:).*(1+gl))  -  kappa.*(Ll_s_sr_or(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:))  -  kappa.*(Ll_s_sr_or(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gl_s_sr_conseff(i,1) = fminbnd(funl,-5,20);

funm = @(gm) abs(sum(...
    (Betas(i,:).*(  log(Cm_s_sr_or(i,:).*(1+gm))  -  kappa.*(Lm_s_sr_or(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:))  -  kappa.*(Lm_s_sr_or(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gm_s_sr_conseff(i,1) = fminbnd(funm,-5,20);

funh = @(gh) abs(sum(...
    (Betas(i,:).*(  log(Ch_s_sr_or(i,:).*(1+gh))  -  kappa.*(Lh_s_sr_or(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:))  -  kappa.*(Lh_s_sr_or(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gh_s_sr_conseff(i,1) = fminbnd(funh,-5,20);
end

% Labor effects
gl_s_sr_labeff = zeros(N,1); 
gm_s_sr_labeff = zeros(N,1);
gh_s_sr_labeff = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum(...
    (Betas(i,:).*(  log(C_or(i,:).*(1+gl))  -  kappa.*(Ll_s_sr_or(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:))  -  kappa.*(L_or(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gl_s_sr_labeff(i,1) = fminbnd(funl,-5,20);

funm = @(gm) abs(sum(...
    (Betas(i,:).*(  log(C_or(i,:).*(1+gm))  -  kappa.*(Lm_s_sr_or(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:))  -  kappa.*(L_or(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gm_s_sr_labeff(i,1) = fminbnd(funm,-5,20);

funh = @(gh) abs(sum(...
    (Betas(i,:).*(  log(C_or(i,:).*(1+gh))  -  kappa.*(Lh_s_sr_or(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:)) -  kappa.*(L_or(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gh_s_sr_labeff(i,1) = fminbnd(funh,-5,20);
end


%% DISPLAY RESULTS PART A %%
 
% Summary statistics total effects
Results = [prctile(gl_s_sr,10), prctile(gm_s_sr,10), prctile(gh_s_sr,10); ...
    prctile(gl_s_sr,50), prctile(gm_s_sr,50), prctile(gh_s_sr,50); ...
    prctile(gl_s_sr,90), prctile(gm_s_sr,90), prctile(gh_s_sr,90); ...
    mean(gl_s_sr), mean(gm_s_sr), mean(gh_s_sr); ...
    std(gl_s_sr), std(gm_s_sr), std(gh_s_sr)];
disp(' RESULTS PART A. Total effects')
disp(Results)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid, high, and no seasonality')
disp(' ')
disp(' ')
 
% Summary statistics consumption effects
Results = [prctile(gl_s_sr_conseff,10), prctile(gm_s_sr_conseff,10), prctile(gh_s_sr_conseff,10); ...
    prctile(gl_s_sr_conseff,50), prctile(gm_s_sr_conseff,50), prctile(gh_s_sr_conseff,50); ...
    prctile(gl_s_sr_conseff,90), prctile(gm_s_sr_conseff,90), prctile(gh_s_sr_conseff,90); ...
    mean(gl_s_sr_conseff), mean(gm_s_sr_conseff), mean(gh_s_sr_conseff); ...
    std(gl_s_sr_conseff), std(gm_s_sr_conseff), std(gh_s_sr_conseff)];
disp(' RESULTS PART A. Consumption effects')
disp(Results)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid, high, and no seasonality')
disp(' ')
disp(' ')

% Summary statistics labor effects
Results = [prctile(gl_s_sr_labeff,10), prctile(gm_s_sr_labeff,10), prctile(gh_s_sr_labeff,10); ...
    prctile(gl_s_sr_labeff,50), prctile(gm_s_sr_labeff,50), prctile(gh_s_sr_labeff,50); ...
    prctile(gl_s_sr_labeff,90), prctile(gm_s_sr_labeff,90), prctile(gh_s_sr_labeff,90); ...
    mean(gl_s_sr_labeff), mean(gm_s_sr_labeff), mean(gh_s_sr_labeff); ...
    std(gl_s_sr_labeff), std(gm_s_sr_labeff), std(gh_s_sr_labeff)];
disp(' RESULTS PART A. Labor effects')
disp(Results)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid, high, and no seasonality')
disp(' ')
disp(' ')

% Graphs
figure
subplot(3,1,1);
hold on
histogram(gl_s_sr,16,'BinWidth',0.01);
hold on
histogram(gl_s_sr_conseff,16,'BinWidth',0.01);
hold on
histogram(gl_s_sr_labeff,16,'BinWidth',0.01);
xlim([-0.05 0.4]);
xlabel('Individual g')
ylabel('Num indiv')
legend('g_{total}','g_{consumption}','g_{labor}')
title('Low seasonality. Positive correlation')
 
subplot(3,1,2);
hold on
histogram(gm_s_sr,16,'BinWidth',0.01);
hold on
histogram(gm_s_sr_conseff,16,'BinWidth',0.01);
hold on
histogram(gm_s_sr_labeff,16,'BinWidth',0.01);
xlim([-0.05 0.4]);
xlabel('Individual g')
ylabel('Num indiv')
legend('g_{total}','g_{consumption}','g_{labor}')
title('Medium seasonality. Positive correlation')
 
subplot(3,1,3);
hold on
histogram(gh_s_sr,16,'BinWidth',0.01);
hold on
histogram(gh_s_sr_conseff,16,'BinWidth',0.01);
hold on
histogram(gh_s_sr_labeff,16,'BinWidth',0.01);
xlim([-0.05 0.4]);
xlabel('Individual g')
ylabel('Num indiv')
legend('g_{total}','g_{consumption}','g_{labor}')
title('high seasonality. Positive correlation')
print('Q2_A','-dpng')


%% PART B : HIGHLY NEGATIVELY CORRELATED CONSUMPTION AND LEISURE %%

% Total effects
gl_s_sr_negative = zeros(N,1); 
gm_s_sr_negative = zeros(N,1);
gh_s_sr_negative = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum(...
    (Betas(i,:).*(  log(Cl_s_sr_or_negative(i,:).*(1+gl))  -  kappa.*(Ll_s_sr_or_negative(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:))  -  kappa.*(L_or(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gl_s_sr_negative(i,1) = fminbnd(funl,-5,20);

funm = @(gm) abs(sum(...
    (Betas(i,:).*(  log(Cm_s_sr_or_negative(i,:).*(1+gm))  -  kappa.*(Lm_s_sr_or_negative(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:))  -  kappa.*(L_or(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gm_s_sr_negative(i,1) = fminbnd(funm,-5,20);

funh = @(gh) abs(sum(...
    (Betas(i,:).*(  log(Ch_s_sr_or_negative(i,:).*(1+gh))  -  kappa.*(Lh_s_sr_or_negative(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:))  -  kappa.*(L_or(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gh_s_sr_negative(i,1) = fminbnd(funh,-5,20);
end

% Consumption effects
gl_s_sr_conseff_negative = zeros(N,1); 
gm_s_sr_conseff_negative = zeros(N,1);
gh_s_sr_conseff_negative = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum(...
    (Betas(i,:).*(  log(Cl_s_sr_or_negative(i,:).*(1+gl))  -  kappa.*(Ll_s_sr_or_negative(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:))  -  kappa.*(Ll_s_sr_or_negative(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gl_s_sr_conseff_negative(i,1) = fminbnd(funl,-5,20);

funm = @(gm) abs(sum(...
    (Betas(i,:).*(  log(Cm_s_sr_or_negative(i,:).*(1+gm))  -  kappa.*(Lm_s_sr_or_negative(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:))  -  kappa.*(Lm_s_sr_or_negative(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gm_s_sr_conseff_negative(i,1) = fminbnd(funm,-5,20);

funh = @(gh) abs(sum(...
    (Betas(i,:).*(  log(Ch_s_sr_or_negative(i,:).*(1+gh))  -  kappa.*(Lh_s_sr_or_negative(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:))  -  kappa.*(Lh_s_sr_or_negative(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gh_s_sr_conseff_negative(i,1) = fminbnd(funh,-5,20);
end

% Labor effects
gl_s_sr_labeff_negative = zeros(N,1); 
gm_s_sr_labeff_negative = zeros(N,1);
gh_s_sr_labeff_negative = zeros(N,1);
for i = 1:N
funl = @(gl) abs(sum(...
    (Betas(i,:).*(  log(C_or(i,:).*(1+gl))  -  kappa.*(Ll_s_sr_or_negative(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:))  -  kappa.*(L_or(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gl_s_sr_labeff_negative(i,1) = fminbnd(funl,-5,20);

funm = @(gm) abs(sum(...
    (Betas(i,:).*(  log(C_or(i,:).*(1+gm))  -  kappa.*(Lm_s_sr_or_negative(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:))  -  kappa.*(L_or(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gm_s_sr_labeff_negative(i,1) = fminbnd(funm,-5,20);

funh = @(gh) abs(sum(...
    (Betas(i,:).*(  log(C_or(i,:).*(1+gh))  -  kappa.*(Lh_s_sr_or_negative(i,:).^(1+1/nu)/(1+1/nu))  )).' - ... 
    (Betas(i,:).*(  log(C_or(i,:)) -  kappa.*(L_or(i,:).^(1+1/nu)/(1+1/nu))                      )).'   ));
gh_s_sr_labeff_negative(i,1) = fminbnd(funh,-5,20);
end

%% DISPLAY RESULTS PART B %%
 
% Summary statistics total effects
Results = [prctile(gl_s_sr_negative,10), prctile(gm_s_sr_negative,10), prctile(gh_s_sr_negative,10); ...
    prctile(gl_s_sr_negative,50), prctile(gm_s_sr_negative,50), prctile(gh_s_sr_negative,50); ...
    prctile(gl_s_sr_negative,90), prctile(gm_s_sr_negative,90), prctile(gh_s_sr_negative,90); ...
    mean(gl_s_sr_negative), mean(gm_s_sr_negative), mean(gh_s_sr_negative); ...
    std(gl_s_sr_negative), std(gm_s_sr_negative), std(gh_s_sr_negative)];
disp(' RESULTS PART B. Total effects')
disp(Results)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid, high, and no seasonality')
disp(' ')
disp(' ')
 
% Summary statistics consumption effects
Results = [prctile(gl_s_sr_conseff_negative,10), prctile(gm_s_sr_conseff_negative,10), prctile(gh_s_sr_conseff_negative,10); ...
    prctile(gl_s_sr_conseff_negative,50), prctile(gm_s_sr_conseff_negative,50), prctile(gh_s_sr_conseff_negative,50); ...
    prctile(gl_s_sr_conseff_negative,90), prctile(gm_s_sr_conseff_negative,90), prctile(gh_s_sr_conseff_negative,90); ...
    mean(gl_s_sr_conseff_negative), mean(gm_s_sr_conseff_negative), mean(gh_s_sr_conseff_negative); ...
    std(gl_s_sr_conseff_negative), std(gm_s_sr_conseff_negative), std(gh_s_sr_conseff_negative)];
disp(' RESULTS PART B. Consumption effects')
disp(Results)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid, high, and no seasonality')
disp(' ')
disp(' ')

% Summary statistics labor effects
Results = [prctile(gl_s_sr_labeff_negative,10), prctile(gm_s_sr_labeff_negative,10), prctile(gh_s_sr_labeff_negative,10); ...
    prctile(gl_s_sr_labeff_negative,50), prctile(gm_s_sr_labeff_negative,50), prctile(gh_s_sr_labeff_negative,50); ...
    prctile(gl_s_sr_labeff_negative,90), prctile(gm_s_sr_labeff_negative,90), prctile(gh_s_sr_labeff_negative,90); ...
    mean(gl_s_sr_labeff_negative), mean(gm_s_sr_labeff_negative), mean(gh_s_sr_labeff_negative); ...
    std(gl_s_sr_labeff_negative), std(gm_s_sr_labeff_negative), std(gh_s_sr_labeff_negative)];
disp(' RESULTS PART B. Labor effects')
disp(Results)
disp('First row: 10th percentile. Second row: Median. Third row: 90th percentile. Fourth row: Means. Fifth row: sd')
disp('Each column: "g" for low, mid, high, and no seasonality')
disp(' ')
disp(' ')

% Graphs
figure
subplot(3,1,1);
hold on
histogram(gl_s_sr_negative,16,'BinWidth',0.01);
hold on
histogram(gl_s_sr_conseff_negative,16,'BinWidth',0.01);
hold on
histogram(gl_s_sr_labeff_negative,16,'BinWidth',0.01);
xlim([-0.05 0.4]);
xlabel('Individual g')
ylabel('Num indiv')
legend('g_{total}','g_{consumption}','g_{labor}')
title('Low seasonality. Negative correlation')
 
subplot(3,1,2);
hold on
histogram(gm_s_sr_negative,16,'BinWidth',0.01);
hold on
histogram(gm_s_sr_conseff_negative,16,'BinWidth',0.01);
hold on
histogram(gm_s_sr_labeff_negative,16,'BinWidth',0.01);
xlim([-0.05 0.4]);
xlabel('Individual g')
ylabel('Num indiv')
legend('g_{total}','g_{consumption}','g_{labor}')
title('Medium seasonality. Negative correlation')
 
subplot(3,1,3);
hold on
histogram(gh_s_sr_negative,16,'BinWidth',0.01);
hold on
histogram(gh_s_sr_conseff_negative,16,'BinWidth',0.01);
hold on
histogram(gh_s_sr_labeff_negative,16,'BinWidth',0.01);
xlim([-0.05 0.4]);
xlabel('Individual g')
ylabel('Num indiv')
legend('g_{total}','g_{consumption}','g_{labor}')
title('high seasonality. Negative correlation')
print('Q2_B','-dpng')