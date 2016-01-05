%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Instantaneous frequency %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inicialization
peak = 0;
inst = 0;
freq_inst = 0;
inst2 = 0;
index = 0;
curr = 0;
curr2 = 0;
I1 = linspace(0,20,500000);
I2 = linspace(20,0,500001);
I = [I1 I2];

%% Peak detection
for i = 2:1:length(V)
    if V(i,2) > 50 && V(i,2) > V(i+1,2) && V(i,2) > V(i-1,2)
        peak = [peak V(i,2)];
        inst = [inst t(i)];
        index = [index i];
        curr = [curr I(i)];
    end
end

peak = peak(2:end);
inst = inst(2:end);
index = index(2:end);
curr = curr(2:end);

%% Firing rate
for j = 1:1:(length(inst)-1)
    freq_inst = [freq_inst 1000/(inst(j+1)-inst(j))];
    inst2 = [inst2 (inst(j+1)+inst(j))/2]; 
    curr2 = [curr2 (curr(j+1)+curr(j))/2];
end

freq_inst = freq_inst(2:end);
inst2 = inst2(2:end);
curr2 = curr2(2:end);
    
%% Figures
figure
plot(inst,peak,'ok')
grid
xlabel('Time (ms)'), ylabel('Membrane Potential (mV)')
title('Spike times')
axis([0 max(t) min(V(:,2)) max(V(:,2))])

figure
plot(inst2,freq_inst,'.-k')
grid
xlabel('Time (ms)'), ylabel('Firing Rate (Hz)')
xlim([0 max(t)])

% fmean = mean(freq_inst(5:end))

figure
plot(curr2,freq_inst,'ok')
grid
xlabel('Injected Current (nA)'), ylabel('Firing rate (Hz)')