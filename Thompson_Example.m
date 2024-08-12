clc; clear all;

p_vec = linspace(0.1,1,10); %vector of arms

a = ones(1, length(p_vec));
b = ones(1,length(p_vec));

T = 100000;
T_phase = 20;
pick = zeros(1,T);
v = 6;
for i = 1:T
    if mod(i, T_phase) ==0
        for j = 1:length(p_vec)
            theta(j) = betarnd(a(j), b(j));
        end
        [val, selected_arm] = max(theta);
        best_selection = 6;
        for k = 1:T_phase
            %Restless Reward
            %Reward = (1 -  (1 -  p_vec(selected_arm)).^T_phase)*Restless_Reward(p_vec(selected_arm), T_phase, v);
            
            %Rested Reward
            Reward = Rested_Reward(p_vec(selected_arm), T_phase, v);


            a(selected_arm) = a(selected_arm) + Reward;
            b(selected_arm) = b(selected_arm) + 1 - Reward;
            pick(i) = selected_arm;
        end
        if selected_arm ~= best_selection
            Regret(i) = 1;
        else
            Regret(i) = 0;
        end
    end
end
Frequencies = pick(pick~=0);
[F1 F2] = hist(Frequencies, length(p_vec));
p_vec(F1 == max(F1))
subplot(1,2,2)
hist(Frequencies, length(p_vec), 'BarWidth', 0.5)
hold on
%bar(p_vec, F1)