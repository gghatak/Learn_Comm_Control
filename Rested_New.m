clc; clear all;

%%
%System Parameters
P = db2pow(24 - 30); %24 dBm transmit power.
fc = 3.2e9; %3.2 GHz carrier frequency.
W = 200e6; %System bandwidth.
K = (3e8/(4 * pi * fc))^2; %Path loss/gain constant.
alpha = 2; %Path loss exponent.
R = 100; %Range of interest.
r = 15; %Distance of the typical controlled system from the typical controller.
Noise = db2pow(-176 - 30)*W; %Noise power
gamma = db2pow(-15); %SINR threshold.
n = 20; %Total horizon length
v =  6; %Successful transmissions needed
%p0_vec = linspace(0.2,0.99,30); %Reliability needed
q_vec = linspace(0.1,1,30); %Channel access probability
MD = zeros(1,length(q_vec));

lambda = 5e-3; %Intensity of controller nodes.

iterations = 10000;
for param = 1:length(q_vec)
    q = q_vec(param)
    p0 = 0.9;
    for iter = 1:iterations
        Actual = 0;
        N_points = poissrnd(pi*lambda*R^2); %Rested
        dist_vec = max(1, R*sqrt(rand(1, N_points)));
        Term1 = exp(-gamma*Noise/(K*P*r^(-alpha)));
        Den = gamma*dist_vec.^(-alpha)/r^(-alpha);
        Term2 = prod(q./(1 + Den) + 1 - q);  %Uncomment for rested        
        Cond_Succ(iter) = Term1.*Term2;        
        n_slots = binornd(n,q); %Uncomment for rested
 %%
 %Uncomment for rested
        for j = v:n_slots
            Actual = Actual + nchoosek(n_slots,j) * (Cond_Succ(iter))^j * (1 - Cond_Succ(iter))^(n_slots - j);
        end     
        if Actual > p0
            MD(param) = MD(param) + 1/iterations;
        end       
   beta_rested = ((p0)/(n - v))^(1/n);   
    end
end
subplot(1,2,1)
plot(q_vec, smooth(MD))
hold on