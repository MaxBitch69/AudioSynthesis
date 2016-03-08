function y = pitchshift(x, shift)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

y = timestretch(x, shift);

[p, q] = rat(shift); % Find fraction corresponding to shift
                     % p = numerator
                     % q = denominator

y = resample(y, q, p);
end

