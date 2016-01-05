function current = Iinj_d(t)
%IINJ Injected current into motoneuron dendrite.
% CURRENT = IINJ(T) returns a current value following the waveform
% selected. There are CONSTANT, PULSE, RAMP and SINUSOID waveforms.

d = 0.0;                                                                   %Time delay
% d1 = 0.0;                                                                  %Time delay
% d2 = 0.0;                                                                  %Time delay

%% Constant
if t >= d
    current = 0.0;
else
    current = 0.0;
end

%% Pulse
% pd1 = 0;                                                                   %Pulse duration
% pd2 = 0;                                                                   %Pulse duration
% 
% if t >= d1 && t < d1+pd1
%     current = 0.0;
% elseif t >= d2 && t < d2+pd2
%     current = 0.0;
% else
%     current = 0.0;
% end

%% Ramp
% if t >= d && t < 5000
%     rs = 0.004;
%     current = rs.*(t-d);
% else
%     rs = 0.004;
%     current = 40-rs.*(t-d);
% end

%% Sinusoid
% if t > d
%   f = 0.0;                                                                 %Sinusoid frequency
%   sa = 0.0;                                                                %Sinusoid amplitude
%   current = sa*sin(2*pi*(f/1000)*(t-d));
% else
%   current = 0;
% end
