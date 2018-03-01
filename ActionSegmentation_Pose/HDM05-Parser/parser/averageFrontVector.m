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

function v = averageFrontVector(mot,varargin)

p1_name = 'lhip';
p2_name = 'rhip';
if (nargin>2)
    p2_name = varargin{2};
end

if (nargin>1)
    p1_name = varargin{1};
end

p1 = trajectoryID(mot,p1_name);
p2 = trajectoryID(mot,p2_name);

n = mot.jointTrajectories{p1} - mot.jointTrajectories{p2};
n = n([1 3],:);
n = n./repmat(sqrt(sum(n.^2)),2,1);
v = mean(n,2);
v = v/sqrt(sum(v.^2));
v = [-v(2);v(1)]; % compute front vector as normal of average direction from lhip to rhip