% Based on 'On the Exact Success Rate of Side Channel Analysis in the
% Gaussian Model' by M. Rivain

clear all
close all

% INPUT PARAMETERS---------------------------------------------------------

% size of plaintext X, of key K and of intermediate value S = phi(X,K)
bit_size = 4;
% number of traces N used in the attack
no_traces = 50;
% univariate noise
sigma = 14.1;
% true key ( 0 <= k_star <= 2^bit_size - 1 )
k_star = 11;
% number of attack simulations using the Rivain approximation
no_simulations = 1000;
% number of DPA attacks
no_DPA_attacks = 32;
% number of traces used to estimate mean and variance
no_estimation_traces = 3000;

% estimate mean (1) or assume perfect identity estimation (2) or assume
% perfect Hamming weight estimation (3)
global estimation_choice
estimation_choice = 2;

%--------------------------------------------------------------------------


if estimation_choice == 1
    
    % Using N random inputs (x_1, x_2, ..., x_N)
    % simulate N traces of univariate identity leakage (l_1,l_2,...,l_N)
    % Subsequently compute means and covariance matrices
    [plaintext_for_estimation, traces_for_estimation, m, C, avgC] = generate_traces(bit_size, no_estimation_traces, sigma, k_star);
end
if estimation_choice == 2
    
    % assume perfect estimation with ID model and homoscedasticity
    m = 0:(2^bit_size)-1;
    C = sigma^2 * ones(1,2^bit_size);
    
end
    
if estimation_choice == 3
        
    % assume perfect estimation with HW model and homoscedasticity
    m = 0:bit_size;
    C = sigma^2 * ones(1,2^bit_size);
    
end

% RHO-DOT parameter estimation---------------------------------------------

if estimation_choice == 1
    
% compute occurence ratio tau(x)
tau = occurence_ratio(plaintext_for_estimation, bit_size);

% compute the mean of prediction per key
meanM = zeros(2^bit_size, 1);
for k=0:2^bit_size-1
    
    for x=0:2^bit_size-1
        meanM(k+1) = meanM(k+1) + tau(x+1) * prediction_function(x, k);
    end   
    
end

% compute the variance of prediction per key
varM = zeros(2^bit_size, 1);
for k=0:2^bit_size-1
    
    for x=0:2^bit_size-1
        varM(k+1) = varM(k+1) + tau(x+1) * (prediction_function(x, k) - meanM(k+1))^2;
    end
     
end


% compute the expected value of rho-dot for every key
meanRhoDot = zeros(2^bit_size, 1);
for k=0:2^bit_size-1
    
    for x=0:2^bit_size-1
        meanRhoDot(k+1) = meanRhoDot(k+1) + tau(x+1) * (prediction_function(x,k) - meanM(k+1)) * m(prediction_function(x,k_star) + 1);
    end
    
    meanRhoDot(k+1) = (1/sqrt(varM(k+1))) * meanRhoDot(k+1);
    
end

% compute the variance of rho-dot for every key
varRhoDot = zeros(2^bit_size);
for k1=0:2^bit_size-1
    
    for k2=0:2^bit_size-1
        
        for x=0:2^bit_size-1
            varRhoDot(k1+1, k2+1) = varRhoDot(k1+1, k2+1) + tau(x+1) * (prediction_function(x,k1) - meanM(k1+1)) * (prediction_function(x,k2) - meanM(k2+1)) * C(prediction_function(x,k_star) + 1);
        end
        
        varRhoDot(k1+1, k2+1) = 1/(no_traces * sqrt(varM(k1+1)) * sqrt(varM(k2+1))) * varRhoDot(k1+1, k2+1);
        
    end 
    
end

end

% RHO-DOT-DOT parameter estimation-----------------------------------------

% compute the expected value of rho-dotdot for every key
meanRhoDotDot = zeros(2^bit_size, 1);
for k=0:2^bit_size-1
    
    for x=0:2^bit_size-1
        meanRhoDotDot(k+1) = meanRhoDotDot(k+1) + prediction_function(x,k) * m(prediction_function(x,k_star) + 1);
    end
    
    meanRhoDotDot(k+1) = 1/(2^bit_size) * meanRhoDotDot(k+1);
    
end

% compute the variance of rho-dotdot for every key
varRhoDotDot = zeros(2^bit_size);
for k1=0:2^bit_size-1
    
    for k2=0:2^bit_size-1
        
        for x=0:2^bit_size-1
            varRhoDotDot(k1+1, k2+1) = varRhoDotDot(k1+1, k2+1) + prediction_function(x,k1) * prediction_function(x,k2) * C(prediction_function(x,k_star) + 1);
        end
        
        varRhoDotDot(k1+1, k2+1) = 1/(no_traces * 2^bit_size) * varRhoDotDot(k1+1, k2+1);
        
    end 
    
end


% RANK ESTIMATION----------------------------------------------------------


% instead of the analytical formula (eq. 15) we opt for the simulation
% approach of section 6.2

if estimation_choice == 1
    meanRho = meanRhoDot;
    varRho = varRhoDot;
else
    meanRho = meanRhoDotDot;
    varRho = varRhoDotDot;
end

% simulate a distiguisher vector
d_sample = mvnrnd(meanRhoDotDot, varRhoDotDot, no_simulations);

% compute the comparison vector using the difference between d_k and d_k*
c_sample = zeros(no_simulations, 2^bit_size);
for i=1:2^bit_size
    c_sample(:,i) = d_sample(:,k_star+1) - d_sample(:,i);
end

% compute the ranking by counting how many negative comparison scores exist
% in every simulation
Ranking = zeros(no_simulations,1);
for i=1:no_simulations
   Ranking(i) = sum((c_sample(i,:)<0))+1;
end

% average rank accross simulations
avgRanking = mean(Ranking);




% SIDE-CHANNEL ATTACK------------------------------------------------------

if (estimation_choice == 1) || (estimation_choice == 2) || (estimation_choice == 3)
    
  
% perform DPA attack
RankingDPA = zeros(no_DPA_attacks,1);
for j = 1:no_DPA_attacks
    
% generate traces for DPA    
[plaintext_for_DPA, traces_for_DPA] = generate_traces(bit_size, no_traces, sigma, k_star);
    
correlation = zeros(2^bit_size,1);
for k=0:2^bit_size-1
    
    v = zeros(no_traces,1);
    for i=1:no_traces % use the number of attack traces for the DPA (not the number of estimation traces)
       v(i) = prediction_function(plaintext_for_DPA(i), k);
    end
   
    correlation(k+1) = corr(v, traces_for_DPA);

    
end

[sorted_corr, sorted_index] = sort(abs(correlation), 'descend');
RankingDPA(j) = find(sorted_index == k_star+1);

end

avgRankingDPA = mean(RankingDPA);

end




