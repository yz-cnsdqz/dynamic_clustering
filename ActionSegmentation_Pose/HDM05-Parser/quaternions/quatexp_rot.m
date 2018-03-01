function P = quatexp_rot(Q)
% P = quatexp_rot(Q)
% Implementation according to Grassia (1998)
%
% input:    Q is a 3xn array of 3-vectors
% output:   P(:,i) is the quaternionic exponential of Q(:,i)

if (size(Q,1)~=3)
    error('Input array: number of rows must be 3!');
end
n = size(Q,2);

thresh = eps^(1/4);

theta = sqrt(sum(Q(:,:).^2)); % euclidean length of the vector part
zero = theta<thresh; %num_zero = length(zero); % identify near-zero vector parts
nonzero = ~zero; % identify nonzero vector parts

vector_parts = zeros(3,n);
%if length(zero>0)
    vector_parts(:,zero) = repmat((0.5 + (1/48)*theta(:, zero).^2),3,1) .* Q(:,zero); % use first two terms of sinc taylor expansion
%end
%if length(nonzero>0)
    vector_parts(:,nonzero) = repmat(sin(0.5*theta(:, nonzero)),3,1).*(Q(:,nonzero)./repmat(theta(:, nonzero),3,1));
%end

real_parts = cos(0.5*theta);

P = [real_parts;vector_parts];