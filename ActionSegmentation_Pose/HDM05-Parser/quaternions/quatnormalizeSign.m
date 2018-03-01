function q = quatnormalizeSign(q)
%% function q = quatnormalizeSign(q)
% INPUT: q, 4xN quaternions
% OUTPUT: q, the input quaternions are modified to lie in the same
%         hemisphere.

if size(q,1) ~= 4 
    error('q must be a 4xN array.');
end

numFrames = size(q,2);

for s=2:numFrames
    if dot(q(:,s-1), q(:,s)) < 0
        q(:,s) = -1 * q(:,s);
    end
end
% the following does not work!!!
% 
% % find the location in the hemisphere (upper hemisphere, s > 0, lower
% % hemisphere, s < 0)
% s = [1 dot(q(:, 1:numFrames-1), q(:, 2:numFrames))];
% s(s<0) = -1;
% s(s>0) = 1;
% sDiff = [0 diff(s)];
% 
% 
% % change the sign of the quats that lie in the lower hemisphere.
% q(:, s<0) = -1*q(:, s<0);