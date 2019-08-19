function [C_pi] = computeExtendFreq(param, vibrato)
% [C_pi] = computeExtendFreq(incremento, vibrato)
% This function compute the binary matrix for spreading the ground-truth
% transcription (real_clase) considering the possible vibrato of some
% instruments.
%
% Inputs:
%       - param: midi parameters
%       - vibrato: overlap

if nargin<2
    vibrato=0;
end

inc = param.midi_inc;
min = param.midi_min;
max = param.midi_max;
p = max-min+1;
p_inc = p*inc;

C_pi = zeros(p,p_inc);

for ii=1:p
    if ii==1
        C_pi( ii, 1:inc/2+vibrato )=1;
    else
        C_pi( ii, ((ii-1)*inc-vibrato-1):(ii*inc+vibrato-2) )=1;
    end 
end
end    