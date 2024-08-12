clc; clear all;
%%
%System Parameters
P = db2pow(24 - 30); %24 dBm transmit power.
fc = 3.2e9; %3.2 GHz carrier frequency.
W = 200e6; %System bandwidth.
K = (3e8/(4 * pi * fc))^2; %Path loss/gain constant.
alpha = 2; %Path loss exponent.
R = 100; %Range of interest.
r = 10; %Distance of the typical controlled system from the typical controller.
Noise = db2pow(-176 - 30)*W; %Noise power
gamma = db2pow(-10); %SINR threshold.
n = 20; %Total horizon length
v =  4; %Successful transmissions needed
%p0_vec = linspace(0.2,0.99,30); %Reliability needed
q_vec = linspace(0.1,1,20); %Channel access probability
MD = zeros(1,length(q_vec));
lambda = 5e-4; %Intensity of controller nodes.

iterations = 10000;
for param = 1:length(q_vec)
    q = q_vec(param)
    p0 = 0.9;
    for iter = 1:iterations
        Actual = 0;
        N_points = poissrnd(q*pi*lambda*R^2); %Restless     
        dist_vec = max(1, R*sqrt(rand(1, N_points)));
        Term1 = exp(-gamma*Noise/(K*P*r^(-alpha)));
        Den = gamma*dist_vec.^(-alpha)/r^(-alpha);
        Term2 = prod(1./(1 + Den)); %Uncomment for restless
        Cond_Succ(iter) = Term1.*Term2;
        n_slots = n; %Uncomment for restless
 %%
    for j = 1:floor((n+1)/(v+1))
        Summand(j) = (-1)^(j+1) * (Cond_Succ(iter) + (n - j*v + 1)/j * (1 - Cond_Succ(iter)))...
            * nchoosek(n - j*v, j -1)*Cond_Succ(iter)^(j*v)*(1 - Cond_Succ(iter))^(j-1);
    end
    Cond_CCDF_Burst(iter) = sum(Summand);
    end
    Mean_CCDF_Burst(param) = mean(Cond_CCDF_Burst);
end
Final = ((1 -  (1 - q_vec).^n).*Mean_CCDF_Burst);
q_vec(Final == max(Final))

subplot(1,2,1)
plot(q_vec, smooth(Final))
hold on