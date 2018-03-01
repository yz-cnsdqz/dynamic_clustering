% This code belongs to the HDM05 mocap database which can be obtained
% from the website http://www.mpi-inf.mpg.de/resources/HDM05 .
%
% If you use and publish results based on this code and data, please
% cite the following technical report:
%
%   @techreport{MuellerRCEKW07_HDM05-Docu,
%     author = {Meinard M{\"u}ller and Tido R{\"o}der and Michael Clausen and Bernd Eberhardt and Bj{\"o}rn Kr{\"u}ger and Andreas Weber},
%     title = {Documentation: Mocap Database {HDM05}},
%     institution = {Universit{\"a}t Bonn},
%     number = {CG-2007-2},
%     year = {2007}
%   }
%
%
% THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
% KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
% PARTICULAR PURPOSE.

function B = gramschmidt(A)
n = size(A,1);
if (size(A,2) ~= n)
    return;
end;
B = zeros(n,n);

B(:,1) = (1/norm(A(:,1)))*A(:,1);

for i=2:n
    v = A(:,i);
    U = B(:,1:i-1); % subspace basis which has already been orthonormalized
    pc = transpose(U)*v; % orthogonal projection coefficients of v onto U
    p = U * pc; % orthogonal projection vector of v onto U
    v = v - p;
    if (norm(v)<eps) % vectors are not linearly independent!
        return;
    end;
    v = v/norm(v);
    B(:,i) = v;
end;