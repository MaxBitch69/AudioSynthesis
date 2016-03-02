function [ phase ] = princarg( phase_in )
%% Principal argument of general phase

phase = mod(phase_in+pi,-2*pi) + pi;

end

