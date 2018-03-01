function quatTarget = resampleQuats(quatSource, fsSource, fsTarget)
% FUNCTION resampleQuats(quatSource, fsSource, fsTarget) resamples
% quaternions according to arbitrary (floating point) sampling rates. It
% does NOT filter the quaternions before resampling. 
% 
% Author: Andreas Baak
% Date: 20.05.09
%
%
% INPUT: quatSource: (4xN double array)
%        fsSource:   sampling rate of the quats
%        fsTarget:   desired target sampling rate

if size(quatSource, 1) ~= 4
    error('input quaternions should be a 4xN matrix.');
end

nframes = size(quatSource, 2);
slope = fsSource/fsTarget;
samples = 1:slope:nframes;
interpCoeff = samples - floor(samples);
samples = floor(samples);

% interpolate between the right samples of the source signal.

quatTarget = slerpNQuats(quatSource(:,samples), quatSource(:,samples+1), interpCoeff, 10^-6);

