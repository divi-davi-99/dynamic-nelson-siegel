%-------------------------------------------------------------------------
%% Data Import and Transformation
%-------------------------------------------------------------------------
DataTable = readtable('ECBData.xlsx');
Yields = table2array(DataTable(:, 2:18));  % Convert only numerical columns
DataTimeTable = table2timetable(DataTable, 'RowTimes', DataTable{:,1});

dates = DataTimeTable.Time; 
maturities = [3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 48, 60, 72, 84, 96, 108, 120]';
series = {'Y3M', 'Y6M', 'Y9M', 'Y12M', 'Y15M', 'Y18M', 'Y21M', 'Y24M', ...
          'Y30M', 'Y36M', 'Y48M', 'Y60M', 'Y72M', 'Y84M', 'Y96M', 'Y108M', 'Y120M'};

Yields_out = Yields(151:end,:);
dates_out = dates(151:end);

Yields = Yields(1:150,:);                   %% ECB estimation sample
dates = dates(1:150,:);                     %% dates - estimation sample
% Yields_val = Yields(125:end, :);          %% ECB validation sample

%-------------------------------------------------------------------------
%%    1st Step: OLS Estimation
%-------------------------------------------------------------------------
lambda0 = 0.0609;
%lambda0 = 0.0300;
X = [ones(size(maturities)) (1-exp(-lambda0*maturities))./(lambda0*maturities) ...
    ((1-exp(-lambda0*maturities))./(lambda0*maturities)-exp(-lambda0*maturities))];

Beta = zeros(size(Yields,1),3);
Residuals = zeros(size(Yields,1),numel(maturities));

for i = 1:size(Yields,1)
    EstMdlOLS = fitlm(X,Yields(i,:)','Intercept',false);
    Beta(i,:) = EstMdlOLS.Coefficients.Estimate';
    Residuals(i,:) = EstMdlOLS.Residuals.Raw';
end

%-------------------------------------------------------------------------
%%    2nd Step: Three AR(1) processes or one VAR(1)
%-------------------------------------------------------------------------

%% THREE UNRESTRICTED AR(1) PROCESSES:
AR_1 = cell(3, 1);
Mdl_AR = arima(1,0,0); % AR(1) for Beta(:,1:3)
AR_1{1} = estimate(Mdl_AR, Beta(:,1));
AR_1{2} = estimate(Mdl_AR, Beta(:,2));
AR_1{3} = estimate(Mdl_AR, Beta(:,3));

%% VAR(1) PROCESS:
MdlVAR = varm(3,1);
EstMdlVAR = estimate(MdlVAR,Beta);
summarize(EstMdlVAR);
Phi = EstMdlVAR.AR{1};

%% EIGENVALUES OF PHI
% Compute eigenvalues of the Phi matrix
eigenvalues = eig(Phi);

% Display eigenvalues and their moduli
disp('Eigenvalues of the Phi matrix:');
disp(eigenvalues);

disp('Moduli of the eigenvalues:');
moduli = abs(eigenvalues); % Compute modulus of each eigenvalue
disp(moduli);

% Stationarity check
if all(moduli < 1)
    disp('The VAR(1) process is stationary: all eigenvalues have modulus < 1.');
else
    disp('The VAR(1) process is NOT stationary: at least one eigenvalue has modulus >= 1.');
end

% Characteristic polynomial determinant
syms lambda; % Define lambda as symbolic variable
I = eye(size(Phi)); % Identity matrix same size as Phi
characteristic_polynomial = det(Phi - lambda * I); % Determinant of characteristic polynomial

% Display the characteristic polynomial
disp('Characteristic polynomial:');
disp(characteristic_polynomial);

%-------------------------------------------------------------------------
%% Estimation - State-space model (Diebold-Rudebusch-Aruoba)
%-------------------------------------------------------------------------

%% VP0: parameter initialization
cm = 3; % number of latent states
vP0 = ones(36,1); % initialize 36 parameters

% Estimates from the 2-step Diebold-Li analysis
A0 = EstMdlVAR.AR{1};     
vP0(1:9) = A0(:);

vP0(10:12) = mean(Beta)'; % mu0

Q0 = EstMdlVAR.Covariance; 
vP0(13:18) = [sqrt(Q0(1,1)); 0; 0; sqrt(Q0(2,2)); 0; sqrt(Q0(3,3))];

H0 = cov(Residuals);       
vP0(19:35) = sqrt(diag(H0));  

vP0(36) = lambda0; % lambda from Diebold-Li

%-------------------------------------------------------------------------
%% Optimization settings
%-------------------------------------------------------------------------
options = optimset('Display','iter','TolX',1e-8,'TolFun',1e-8,...
                'Diagnostics','off', 'MaxIter',10000, 'MaxFunEvals', 10000,...
                'LargeScale', 'off', 'PlotFcns', @optimplotfval);

%-------------------------------------------------------------------------
%% Constraints
%-------------------------------------------------------------------------
% Lower bound
lb = -inf*ones(36,1);
lb(1) = 0.7;
lb(5) = 0.7;
lb(9) = 0.7;
lb(10:12) = -5.00;
lb(13:18) = -inf;
lb(19:35) = 0.001;
lb(36) = 0.001;   

% Upper bound
ub = +inf*ones(36,1); 
ub(1) = 0.99;
ub(5) = 0.99;
ub(9) = 0.99;
ub(10:12) = 10.00; % mu

%-------------------------------------------------------------------------
%% Maximum likelihood estimation function
%-------------------------------------------------------------------------
[vP,fval] = fmincon(@(vP) -ZKalmanFilter2(vP,Yields,maturities), ...
    vP0,[],[],[],[],lb,ub,[],options);

%-------------------------------------------------------------------------
%% Estimates obtained via Kalman Filter
%-------------------------------------------------------------------------
A = reshape(vP(1:9),[3,3]);
mu = vP(10:12);
a = vP(13:18);
q = triu(ones(3),0);
q(q==1) = a;
q=q'; 
Q=q*q'; 
% Q = [vpEst(13:15)'; 0 vpEst(16:17)'; 0 0 vpEst(18)];
R = diag(vP(19:35).^2);
lambda = vP(36);

loading_L = ones(17, 1);
loading_S = (1 - exp(-lambda * maturities)) ./ (lambda * maturities);
loading_C = loading_S - exp(-lambda * maturities);
C = [loading_L, loading_S, loading_C];

[se, t_stats, p_values] = KalmanInference(vP, Yields, maturities);

loading_C_max = loading_S - exp(-lambda * maturities);
loading_C_DL = (1 - exp(-lambda0 * maturities)) ./ (lambda0 * maturities)...
                    - exp(-lambda0 * maturities);

intercept = C * mu;
DeflatedYields = Yields - intercept';

[smoothed_beta, smoothed_Beta] = KFStateSmoother3(vP, DeflatedYields, maturities);
EstimatedStates = (smoothed_beta + mu)';

EstimatedYields = EstimatedStates * C';
Yields_2step = Beta * X';

figure
hold on
plot(dates, EstimatedYields(:,1), '-', 'LineWidth', 2, 'Color', 'r')
plot(dates, Yields(:,1), ':', 'LineWidth', 2, 'Color', 'k')
plot(dates, Yields_2step(:,1), '-', 'LineWidth', 2, 'Color', 'b')
legend('SSM', 'Observed', '2-step')

MSE_SSM = sum((EstimatedYields - Yields).^2);
MSE_2step = sum((Yields_2step - Yields).^2);
RMSE_SSM = sqrt(MSE_SSM);
RMSE_2step = sqrt(MSE_2step);

figure
hold on
plot(maturities, RMSE_SSM, 'Color', '#A2142F', 'LineWidth', 8)
plot(maturities, RMSE_2step, 'Color', 'k', 'LineWidth', 8)
ax = gca;
ax.FontSize = 16;
ylabel('RMSE')
xlabel('Maturities (in months)')
legend('SSM', '2-step')
grid on


% -------------------------------------------------------------------------
%% FUNCTIONS
% -------------------------------------------------------------------------
function L = ZKalmanFilter2(vP, Yields, maturities)

NoM = 17;                       % Number of maturities
A = reshape(vP(1:9), [3, 3]);   % 3x3 transition matrix
mu = vP(10:12);                 % Vector of state means

% Construction of lower triangular Q
mask = tril(true(3));
Q = zeros(3);
Q(mask) = vP(13:18);
Q = Q * Q';                     % Q is semi-definite positive

H = diag(vP(19:35));            % Diagonal covariance matrix of observation noise, 17x17
lambda = vP(36);                % Nelson-Siegel decay parameter

%% Measurement equation (Nelson-Siegel)
% Construct matrix C (17x3) using input maturities
loading_L = ones(length(maturities), 1);  % Level
loading_S = (1 - exp(-lambda * maturities)) ./ (lambda * maturities);  % Slope
loading_C = loading_S - exp(-lambda * maturities);                     % Curvature
C = [loading_L, loading_S, loading_C];  % 17x3 matrix

%% Kalman Filter implementation
ss = 3;  % Number of latent states (L, S, C)
[T, ~] = size(Yields);  % Yields is T×17 (e.g., 245×17)

% Matrix initialization
beta_pred = zeros(ss, T);
Beta_pred = zeros(ss, ss, T);
dv = zeros(length(maturities), T);  % Innovations (17×T)
loglik = zeros(1, T);

% Initial state (stationary distribution)
prev_beta = mu;
prev_Beta = dlyap(A, Q);  % Solves A*prev_Beta*A' + Q = prev_Beta

for t = 1:T
    %% Step 1: Prediction
    beta_pred(:, t) = A * prev_beta + (eye(ss) - A) * mu;
    Beta_pred(:, :, t) = A * prev_Beta * A' + Q;
    
    %% Step 2: Update
    % Extract current observations (row t of Yields, converted to column)
    y_t = Yields(t, :)';  % Convert 1×17 to 17×1
    
    % Compute innovation
    dv(:, t) = y_t - C * beta_pred(:, t);  % 17×1
    
    % Innovation covariance matrix
    dF = C * Beta_pred(:, :, t) * C' + H;  % 17×17
    
    % Inverse computation
    F_inv = inv(dF);
    
    % Kalman gain
    K = Beta_pred(:, :, t) * C' * F_inv;  % 3×17
    
    % State estimate update
    beta_new(:, t) = beta_pred(:, t) + K * dv(:, t);
    Beta_new(:, :, t) = (eye(ss) - K * C) * Beta_pred(:, :, t);
    
    % Prepare for next iteration
    prev_beta = beta_new(:, t);
    prev_Beta = Beta_new(:, :, t);
    
    % Log-likelihood
    loglik(t) = -0.5 * (log(det(dF)) + dv(:, t)' * F_inv * dv(:, t) + 17 * log(2*pi));
end

L = sum(loglik);
end



function [smoothed_beta, smoothed_Beta] = KFStateSmoother3(vP, DeflatedYields, maturities)
    
    % Inizializzazione
    A = reshape(vP(1:9), [3, 3]);
    mu = vP(10:12);
    mask = tril(true(3));
    Q = zeros(3);
    Q(mask) = vP(13:18);
    Q = Q * Q';
    H = diag(vP(19:35));
    lambda = vP(36);
    
    % Matrice di misura
    loading_L = ones(length(maturities), 1);
    loading_S = (1 - exp(-lambda * maturities)) ./ (lambda * maturities);
    loading_C = loading_S - exp(-lambda * maturities);
    C = [loading_L, loading_S, loading_C];
    
    [T, n_yields] = size(DeflatedYields);
    ss = 3;  % numero di stati latenti
    
    % Inizializzazione delle matrici
    aaf = zeros(ss, T);                     % (beta_{t|t-1})
    aP = zeros(ss, ss, T);                  % (Beta_{t|t-1}) – matrice di covarianza
    av = zeros(n_yields, T);                % Innovazione
    aF = zeros(n_yields, n_yields, T);      % Varianza dell'innovazione
    aK = zeros(ss, n_yields, T);            % Guadagno di Kalman
    
    % COndizioni iniziali
    prev_beta = mu;
    prev_Beta = dlyap(A, Q);
    
   
    for t = 1:T
       
        aaf(:,t) = prev_beta;
        aP(:,:,t) = prev_Beta;
        
        % Innovazione
        y_t = DeflatedYields(t,:)';
        dv = y_t - C * prev_beta;
        dF = C * prev_Beta * C' + H;
        
        % Guadagno di Kalman
        K = A * prev_Beta * C' / dF;
        
        av(:,t) = dv;
        aF(:,:,t) = dF;
        aK(:,:,t) = K;
        
        % Aggiornamento
        prev_beta = A * prev_beta + (eye(ss) - A) * mu + K * dv;
        prev_Beta = A * prev_Beta * A' + Q - K * dF * K';
    end
    
    % Algoritmo di smussamento
    vr = zeros(ss,1); mN = zeros(ss,ss);
    smoothed_beta = zeros(ss, T); smoothed_Beta = zeros(ss, ss, T);
    
    for t = T:-1:1
        dFinv = aF(:,:,t) \ eye(n_yields);
        mL = A - aK(:,:,t) * C;
        vr = C' * dFinv * av(:,t) + mL' * vr;
        mN = C' * dFinv * C + mL' * mN * mL;
        
        prev_beta = aaf(:,t);
        prev_Beta = aP(:,:,t);
        
        smoothed_beta(:,t) = prev_beta + prev_Beta * vr;
        smoothed_Beta(:,:,t) = prev_Beta - prev_Beta * mN * prev_Beta;
    end
end
