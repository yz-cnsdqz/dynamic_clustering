function Q = slerpNQuats(varargin)
% Q = slerp(Q1,Q2,u,tol)
% Computes the spherical linear interpolation between pairs of Q1 and Q2 for the
% parameters u \in [0,1]
%
% Input:    * Q1,Q2 are 4xN quaternion matrices represented as column vectors.
%             Non-unit quaternions will be normalized prior to
%             interpolation!
%           * u is a vector from [0,1]^N specifying interpolation parameter
%             values for each quaternion pair (Q1(:,n), Q2(:,n))
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
        Q1 = varargin{1};
        Q2 = varargin{2};
        u = varargin{3};
        tol = 10*eps;
        error_on_deviations = true;
    case 4
        Q1 = varargin{1};
        Q2 = varargin{2};
        u = varargin{3};
        tol = varargin{4};
        error_on_deviations = false;
    otherwise
        error('Wrong number of arguments.');
end

if (size(Q1,1) ~= 4 || size(Q2,1)~= 4 || any(quatnormsq(Q1)<eps) || any(quatnormsq(Q2)<eps) || size(Q1, 2) ~= size(Q2, 2))
    error('Input quaternions have to be given as 4xN nonzero vectors!');
end

if (any(abs(quatnormsq(Q1)-1)>tol))
    if (error_on_deviations)
        error('Some input quaternions Q1 are not of unit length!');
    else
        warning('Some input quaternions Q1 are not of unit length! q1 will be re-normalized.');
    end
    Q1 = quatnormalize(Q1);
end

if (any(abs(quatnormsq(Q2)-1)>tol))
    if (error_on_deviations)
        error('Some input quaternions Q2 are not of unit length!');
    else
        warning('Some input quaternion Q2 are not of unit length! Q2 will be re-normalized.');
    end
    Q2 = quatnormalize(Q2);
end
    
N = size(Q1, 2);

%Q = zeros(4,N);
%for i=1:N
   Q = quatmult(Q1,quatexp(u([1;1;1],:).*quatlog(quatmult(quatinv(Q1),Q2),tol)));
%     costheta = transpose(q1)*q2;
%     sintheta = sqrt(1-costheta*costheta);
%     if (abs(sintheta)<eps)
%         continue;
%     end;
%     theta = atan2(sintheta,costheta);
%     Q(:,i) = (sin((1-u(i))*theta) / sintheta) * q1 + (sin(u(i)*theta) / sintheta) * q2;
%end