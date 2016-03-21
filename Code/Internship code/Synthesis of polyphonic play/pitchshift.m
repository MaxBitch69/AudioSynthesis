function y = pitchshift(x, shift)
% Pitch shifting with STFT of signal x using percentage shift
% Does a time stretching first of x, then resample it to obtain the 
% right note.

y = timestretch(x, shift);

[p, q] = rat(shift); % Find fraction corresponding to shift
                     % p = numerator
                     % q = denominator


y = resample(y, q, p);

end

