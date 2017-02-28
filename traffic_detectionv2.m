%% Road congestion algorithm simulation. Simulation on US-1 South. The maximum speed on US1 south is 65mph. We assume three hypothesis, no congestion - H0 (null hypothesis), medium congestion - H1, high congestion - H2

% H0 mean - 65 std dev - 10
% H1 mean - 45 std dev - 7
% H2 mean - 20 std dev - 12

sample_time = 2; % is the time interval between two gps readings
% assumption is the hypothesis is valid for 1 miles
% find the sample numbers corresponding to a 1 mile stretch
% sample_number_h0 = 56/sample_time;
% sample_number_h1 = 80/sample_time;
% sample_number_h2 = 60/sample_time;

mean1 = 45;
sigma1 = 25;
mean2 = 40;
sigma2 = 25;
sample_number = 500;

rep_num = 1000;
%% Calculating thresholds for the likelihood ratio test

% Pfa21 = 1 - Q((thresh21 - mean1)/sigma1)

Pfa21 = 0.5;
thresh21 = (qfuncinv(1 - Pfa21)*sigma1/sqrt(sample_number)) + mean1; % thresh between H2 and H1 
%%
simulated_prob = 0;
for j = 1:rep_num
    speed_h2 = 35 + sigma2.*randn(1,sample_number);
    
%     for i = 1:sample_number
%         if speed_h2(i)<0
%             speed_h2(i)=0;
%         end
%     end

    distance_from_reference = zeros(1,sample_number);

    distance_from_reference(1) = speed_h2(1)*sample_time;
    for i = 2:sample_number
        distance_from_reference(i) = distance_from_reference(i-1) + speed_h2(i)*sample_time;
    end
    %% Estimating the speeds from the simulated position data.Speed =
    %% (Position2 - Position1)/2
    speed = zeros(1,sample_number);

    speed(1) = distance_from_reference(1)/sample_time;
    for i = 2:sample_number
        speed(i) = (distance_from_reference(i) - distance_from_reference(i-1))/sample_time;
    end

    %% LRT. First compare H2 and H1. If we decide H2 then we are good. If H1 is
    %% then compare H1 and H0. This will give the final decision

    % i = 2 and j = 1
    if (sigma2^2*sum((speed - mean1).^2) - sigma1^2*sum((speed - mean2).^2)) > 2*sigma1^2*sigma2^2*log(((sigma2/sigma1)^sample_number)*thresh21)/log(2.71828)
        simulated_prob = simulated_prob + 1;
    end
end

Pdetection = simulated_prob/rep_num;