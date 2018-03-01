function Q = slerp(varargin)
% Q = slerp(q1,q2,u,tol)
% Computes the spherical linear interpolation between q1 and q2 for the
% parameter u \in [0,1]
%
% Input:    * q1,q2 are quaternions represented as column vectors.
%             Non-unit quaternions will be normalized prior to
%             interpolation!
%           * u is a vector from [0,1]^n specifying interpolation parameter
%             values at n positions
%           * If tol is specified, it denotes the allowable deviation from unit 
%             length for the input quaternions. Deviations less than tol will
%             lead to re-normalization, deviations larger than tol will lead
%             to re-normalization and a warning. If tol is not specified,
%             all deviations larger than the machine precision will lead to
%             an error.
%
% Output:   Q is a 4xn array of spherically interpolated unit quaternions

switch (nargin)
    case 3
        q1 = varargin{1};
        q2 = varargin{2};
        u = varargin{3};
        tol = 10*eps;
        error_on_deviations = true;
    case 4
        q1 = varargin{1};
        q2 = varargin{2};
        u = varargin{3};
        tol = varargin{4};
        error_on_deviations = false;
    otherwise
        error('Wrong number of arguments.');
end

if (size(q1,1) ~= 4 || size(q2,1)~= 4 || quatnormsq(q1)<eps || quatnormsq(q2)<eps)
    error('Input quaternions have to be given as 4x1 nonzero vectors!');
end;

if (abs(quatnormsq(q1)-1)>tol)
    if (error_on_deviations)
        error('Input quaternion q1 is not of unit length!');
    else
        warning('Input quaternion q1 is not of unit length! q1 will be re-normalized.');
    end
    q1 = quatnormalize(q1);
end

if (abs(quatnormsq(q2)-1)>tol)
    if (error_on_deviations)
        error('Input quaternion q2 is not of unit length!');
    else
        warning('Input quaternion q2 is not of unit length! q2 will be re-normalized.');
    end
    q2 = quatnormalize(q2);
end
    
n = length(u);

Q = zeros(4,n);
for i=1:n
   Q(:,i) = quatmult(q1,quatexp(u(i)*quatlog(quatmult(quatinv(q1),q2),tol)));
%     costheta = transpose(q1)*q2;
%     sintheta = sqrt(1-costheta*costheta);
%     if (abs(sintheta)<eps)
%         continue;
%     end;
%     theta = atan2(sintheta,costheta);
%     Q(:,i) = (sin((1-u(i))*theta) / sintheta) * q1 + (sin(u(i)*theta) / sintheta) * q2;
end