function P = quatnormalize(Q)
% input:    Q is a 4xn array of quaternions represented as column vectors
% output:   P(:,i) is a unit quaternion in the direction of Q(:,i), 

if (size(Q,1)~=4)
    error('Input array: number of rows must be 4!');
end

nsq = sum(Q.^2);
P = Q./repmat(sqrt(nsq),4,1);