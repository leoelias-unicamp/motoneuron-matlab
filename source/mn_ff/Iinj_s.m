function current = Iinj_s(t)
%IINJ Injected current into motoneuron soma.
% CURRENT = IINJ(T) returns a current value following the waveform
% selected. There are CONSTANT, PULSE, RAMP and SINUSOID waveforms.

% d = 0.0;                                                                   %Time delay
d1 = 0.0;                                                                  %Time delay
d2 = 5000.0;                                                                  %Time delay

%% Constant
% if t >= d
%     current = 0.0;
% else
%     current = 0.0;
% end

%% Pulse
pd1 = 2000.0;                                                                 %Pulse duration
pd2 = 500.0;                                                                 %Pulse duration

if t >= d1 && t < d1+pd1
    current = 20;
elseif t >= d2 && t < d2+pd2
    current = -5.0;
else
    current = 0.0;
end

%% Ramp
% if t >= d && t < 6000
%     rs = 0.005;
%     current = rs.*(t-d);
% else
%     rs = 0.005;
%     current = 60-rs.*(t-d);
% end

%% Sinusoid
% if t > d
%   f = 10;                                                                  %Sinusoid frequency
%   sa = 0.5;                                                                %Sinusoid amplitude
%   current = 2.5+sa*sin(2*pi*(f/1000)*(t-d));
% else
%   current = 0;
% end
