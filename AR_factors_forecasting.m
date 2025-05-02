%% FORECAST AR(1) ON FACTORS (DNS)

HrznSet = [1 6 12]; 
  HrznChar = char('step1','step6','step12');    %horizon
  DNS_AR1_STORE = struct('step1',{},'step6',{},'step12',{});

% Construct structure for forecast error storage
  for i_h = 1:3   % 3 forecast horion.
  Hrzn = HrznSet(1,i_h);
  Hrzn_name = strtrim(HrznChar(i_h,:));
  DNS_AR1_STORE(1).(Hrzn_name) = struct([]);
  tmp_output = forecast_DNS_AR1(dates, Yields, Hrzn);
  DNS_AR1_STORE(1).(Hrzn_name) = tmp_output;
  end
% save to result folder  
  save('.\result\AR1_STORE','DNS_AR1_STORE') 

  rmse_results = AR1_RMSE(DNS_AR1_STORE);

  %% RMSE
MatChar = char('m3','y1','y3','y5','y10');
stepChar = char('step1','step6','step12');
RMSE_step = zeros(5,1);
for i_step = 1:3
    step_set = strtrim(stepChar(i_step,:));
    for i_mat = 1:5
        set = strtrim(MatChar(i_mat,:));
        RMSE_step(i_mat,i_step) = rmse_results.(step_set).(set).rmse;  
    end
end

% Salva i risultati RMSE
  save('.\result\RMSE_DNSa', 'rmse_results');

%% STATS
% STATS
set_mat = ["m3", "y1", "y3", "y5", "y10"]; 
steps = ["step1", "step6", "step12"];

% MEDIA E STD.DEV.
% Tutti gli step (1, 6, 12)

mean_values = zeros(3, 5);
std_values = zeros(3, 5);

for t = 1:3
    step = steps(t); % Seleziona lo step attuale
    for i = 1:5
        mat_name = set_mat(i); % Seleziona il nome della variabile come stringa
        stat = DNS_AR1_STORE.(step).(mat_name).error(:,1); % Accesso dinamico ai dati
        
        mean_values(t, i) = mean(stat); % Memorizza la media
        std_values(t, i)  = std(stat);  % Memorizza la deviazione standard
    end
end

% disp('Medie:');
% disp(mean_values');
% disp('Deviazioni standard:');
% disp(std_values')


%                       AUTOCORRELAZIONE
% STEP-1
acf_selected1 = zeros(2,5); % Preallocazione della matrice dei risultati

for i = 1:5
    var_name = set_mat(i); % Seleziona il nome della variabile come stringa
    data = DNS_AR1_STORE.step1.(var_name).error(:,1); % Accesso dinamico ai dati
    
    acf_values = autocorr(data, 'NumLags', 12); % Calcolare fino al lag 12
    acf_selected1(:, i) = acf_values([2, 13]); % Selezionare i valori ai lag 1 e 12
end

% disp('Autocorrelazione ai lag 1 e 12 per ogni serie:');
% disp(acf_selected1');

% STEP-6
acf_selected6 = zeros(2,5); % Preallocazione della matrice dei risultati

for i = 1:5
    var_name = set_mat(i); % Seleziona il nome della variabile come stringa
    data = DNS_AR1_STORE.step6.(var_name).error(:,1); % Accesso dinamico ai dati
    
    acf_values = autocorr(data, 'NumLags', 18); % Calcolare fino al lag 12
    acf_selected6(:, i) = acf_values([7, 19]); % Selezionare i valori ai lag 1 e 12
end

% disp('Autocorrelazione ai lag 6 e 18 per ogni serie:');
% disp(acf_selected6');

% STEP-12
acf_selected12 = zeros(2,5); % Preallocazione della matrice dei risultati

for i = 1:5
    var_name = set_mat(i); % Seleziona il nome della variabile come stringa
    data = DNS_AR1_STORE.step12.(var_name).error(:,1); % Accesso dinamico ai dati
    
    acf_values = autocorr(data, 'NumLags', 24); % Calcolare fino al lag 12
    acf_selected12(:, i) = acf_values([13, 25]); % Selezionare i valori ai lag 1 e 12
end

% disp('Autocorrelazione ai lag 12 e 24 per ogni serie:');
% disp(acf_selected12');



% TABELLA LATEX
%% Out-of-sample 1-month-ahead forecasting results
disp('Out-of-sample 1-month-ahead forecasting results');
disp([mean_values(1,:)' std_values(1,:)' RMSE_step(:,1) acf_selected1']);

%% Out-of-sample 6-month-ahead forecasting results
disp('Out-of-sample 6-month-ahead forecasting results');
disp([mean_values(2,:)' std_values(2,:)' RMSE_step(:,2) acf_selected6']);

%% Out-of-sample 12-month-ahead forecasting results
disp('Out-of-sample 12-month-ahead forecasting results');
disp([mean_values(3,:)' std_values(3,:)' RMSE_step(:,3) acf_selected12']);


%% FUNZIONI

function [ output ] = forecast_DNS_AR1(dates, Yields, Hrzn)
 y = Yields;
 horizon = Hrzn;

 MatSet =[1 4 10 12 17];
 MatChar = char('m3','y1','y3','y5','y10') ;
 [yr,mth,day] = datevec(dates);
 time = yr*10000+mth*100+day;

 id_start = find(time==19940131);    % Diebold-Li
 id_end = find(time==20001229);      % Diebold-Li
 % id_start = find(time==20170301);  % ECB
 % id_end = find(time==20250101);    % ECB

 fDates = dates(id_start: id_end,1);
 bffrDates = dates(1:(id_start-1),1);
 fLength=size(fDates,1);
 bffrLength = size(bffrDates,1);
 
% BETA
 beta_17 = DNSbeta(y,17);
% FATTORI
 fac_5 = DNSloading(5);

% Generate point forecast table
 yHat = zeros(fLength,1,5);
    for i_F = 1 : fLength
         forecastPnt = i_F + bffrLength;
        % Get beta forecasts
         fBeta_17 = Model_forecast_DNS_AR1(beta_17,forecastPnt,horizon);
        % Get yield forecasts
         yHat(i_F,1,:) = fBeta_17 * fac_5';
    end
    % Generate forecast error table: fLength * 1 (diffusion indexes) for 5 'Mat'
    for i_m =1:5
         Mat = MatSet(1,i_m);
         Mat_name = strtrim(MatChar(i_m,:));
         single = Yields(bffrLength+1:end,Mat);
         panel = repmat(single,1,1);
         tmpError = panel-yHat(:,:,i_m);
         output.(Mat_name).error = tmpError;
    end
end

%% MODEL

function [forecastBetas] = Model_forecast_DNS_AR1(beta_17, fPoint, Hrzn)
% Get Beta (NS factor) forecast under an AR(1) specification.
% This function uses an expanding window approach for direct forecasting
% starting from the beginning of the sample up to fPoint-1
    Betas = beta_17;
    forecastBetas = zeros(1,3);
    
    for j = 1:3
        % Predittore lagged (y_{t-Hrzn})
        X = Betas(1:(fPoint-Hrzn-1), j); 
        
        % Variabile dipendente (y_t)
        Y = Betas((Hrzn+1):(fPoint-1), j); 
        
        % Matrice disegno con intercetta
        X_design = [ones(length(X), 1), X]; 
        
        % Stima dei coefficienti AR(1)
        % Opzione 1: usando regress()
        tmpGamma = regress(Y, X_design);
        
        % Opzione 2: usando fitlm() - commentata ma corretta
        % model = fitlm(X, Y, 'linear');
        % tmpGamma = model.Coefficients.Estimate;
        
        % Calcolo della previsione
        forecastBetas(1,j) = tmpGamma(1) + tmpGamma(2)* Betas(fPoint-Hrzn, j);
    end
end


%% DNS-BETA
function [ beta ] = DNSbeta(Yields, yDim)
% Estimate latent factor Beta with 10,6,4-dimension bond yield data.
t = size(Yields,1);
if yDim==17, y_17 = Yields; 
factors_17 = DNSloading(17);     
for i = 1:t
beta_17(i,:) = regress(y_17(i,:)', factors_17)';
end
beta = beta_17;

elseif yDim==6, y_6=horzcat(Yields(:,1),Yields(:,2),Yields(:,3),Yields(:,5),Yields(:,7),Yields(:,10));
factors_6 = DNSloading(6);
for i = 1:t
beta_6(i,:) = regress(y_6(i,:)', factors_6)';
end
beta = beta_6;
elseif yDim==4, y_4=horzcat(Yields(:,1),Yields(:,3),Yields(:,5),Yields(:,10)); 
factors_4 = DNSloading(4);
for i = 1:t
beta_4(i,:) = regress(y_4(i,:)', factors_4)';
end
beta = beta_4;
else
    warning('Invalid yield data dimension.')
end
end

%% DNS-LOADING
function [ factors ] = DNSloading(yDim)
% Estimate latent factor Beta.
% Parameter setup.
lambda = 0.0609;   % Diebold-Li

if yDim == 17
TimeToMat_17 = [3; 6; 9; 12; 15; 18; 21; 24; 30; 36; 48; 60; 72; 84; 96; 108; 120];
factors_17 = [ones(size(TimeToMat_17)) (1 - exp(-lambda*TimeToMat_17))./(lambda*TimeToMat_17) ...
                ((1 - exp(-lambda*TimeToMat_17))./(lambda*TimeToMat_17) - exp(-lambda*TimeToMat_17))];
factors = factors_17;

elseif yDim == 6          
TimeToMat_6 = [12; 24; 36; 60; 84; 120];
factors_6 = [ones(size(TimeToMat_6)) (1 - exp(-lambda*TimeToMat_6))./(lambda*TimeToMat_6) ...
                ((1 - exp(-lambda*TimeToMat_6))./(lambda*TimeToMat_6) - exp(-lambda*TimeToMat_6))];       
factors = factors_6;

elseif yDim == 5          
TimeToMat_5 = [3; 12; 36; 60; 120];
factors_5 = [ones(size(TimeToMat_5)) (1 - exp(-lambda*TimeToMat_5))./(lambda*TimeToMat_5) ...
                ((1 - exp(-lambda*TimeToMat_5))./(lambda*TimeToMat_5) - exp(-lambda*TimeToMat_5))];       
factors = factors_5;

elseif yDim == 4
TimeToMat_4 = [12; 36; 60; 120];
factors_4 = [ones(size(TimeToMat_4)) (1 - exp(-lambda*TimeToMat_4))./(lambda*TimeToMat_4) ...
                ((1 - exp(-lambda*TimeToMat_4))./(lambda*TimeToMat_4) - exp(-lambda*TimeToMat_4))];       
factors = factors_4;

else
     warning('Invalid yield data dimension.')

end

end

%% AR1_RMSE
function results = AR1_RMSE(DNS_AR1_STORE)
% Calcola il Root Mean Square Error (RMSE) per tutti gli orizzonti di previsione e maturità
% Input: DNS_AR_STORE struttura contenente gli errori di previsione
% Output: results struttura contenente RMSE per ogni orizzonte e maturità

    % Inizializza struttura dei risultati
    results = struct();
    
    % Ottieni i nomi degli orizzonti (step1, step3, step12)
    horizonNames = fieldnames(DNS_AR1_STORE);
    
    % Per ogni orizzonte di previsione
    for h = 1:length(horizonNames)
        horizonName = horizonNames{h};
        results.(horizonName) = struct();
        
        % Ottieni i nomi delle maturità (m3, y1, y3, y5, y10)
        maturityNames = fieldnames(DNS_AR1_STORE.(horizonName));
        
        % Per ogni maturità
        for m = 1:length(maturityNames)
            maturityName = maturityNames{m};
            
            % Estrai gli errori di previsione
            errors = DNS_AR1_STORE.(horizonName).(maturityName).error;
            
            % Calcola RMSE per ogni colonna (se ci sono più colonne)
            [rows, cols] = size(errors);
            rmse_values = zeros(1, cols);
            
            for c = 1:cols
                % RMSE = sqrt(mean(errors^2))
                rmse_values(c) = sqrt(mean(errors(:,c).^2));
            end
            
            % Memorizza i risultati
            results.(horizonName).(maturityName).rmse = rmse_values;
            
            % Se ci sono più colonne, calcola anche la media
            if cols > 1
                results.(horizonName).(maturityName).rmse_mean = mean(rmse_values);
            end
        end
    end
    
    % Stampa una tabella di riepilogo
    disp('Root Mean Square Error (RMSE) Summary:');
    disp('======================================');
    
    for h = 1:length(horizonNames)
        disp(['Horizon: ', horizonNames{h}]);
        
        % Prepara una tabella per questo orizzonte
        maturityNames = fieldnames(DNS_AR1_STORE.(horizonNames{h}));
        
        % Intestazioni
        header = '| Maturity |';
        for c = 1:size(DNS_AR1_STORE.(horizonNames{h}).(maturityNames{1}).error, 2)
            header = [header, ' Model ', num2str(c), ' |'];
        end
        if size(DNS_AR1_STORE.(horizonNames{h}).(maturityNames{1}).error, 2) > 1
            header = [header, ' Mean |'];
        end
        disp(header);
        
        separator = repmat('-', 1, length(header));
        disp(separator);
        
        % Valori
        for m = 1:length(maturityNames)
            row = ['| ', maturityNames{m}, ' '.^max(1, 9-length(maturityNames{m})), '|'];
            
            % Valori RMSE per modello
            for c = 1:length(results.(horizonNames{h}).(maturityNames{m}).rmse)
                row = [row, ' ', sprintf('%7.4f', results.(horizonNames{h}).(maturityNames{m}).rmse(c)), ' |'];
            end
            
            % Media (se applicabile)
            if size(DNS_AR1_STORE.(horizonNames{h}).(maturityNames{1}).error, 2) > 1
                row = [row, ' ', sprintf('%7.4f', results.(horizonNames{h}).(maturityNames{m}).rmse_mean), ' |'];
            end
            
            disp(row);
        end
        
        disp(separator);
        disp(' ');
    end
end
