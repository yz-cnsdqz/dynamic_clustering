function P = quatconj(Q)
% input:    Q is a 4xn array of quaternions represented as column vectors
% output:   P(:,i) is the conjugate quaternion of Q(:,i)

if (size(Q,1)~=4)
    error('Input array: number of rows must be 4!');
end
n = size(Q,2);

P = Q .* [ones(1, n);-ones(3,n)];