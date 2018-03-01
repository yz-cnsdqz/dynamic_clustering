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

function boundingBox = computeBoundingBox(mot)

boundingBox = [inf;-inf;inf;-inf;inf;-inf];

for k = 1:mot.njoints
    trajectory = mot.jointTrajectories{k};
	boundingBox(1) = min(boundingBox(1),min(trajectory(1,:))); % xmin
	boundingBox(2) = max(boundingBox(2),max(trajectory(1,:))); % xmax
	boundingBox(3) = min(boundingBox(3),min(trajectory(2,:))); % ymin
	boundingBox(4) = max(boundingBox(4),max(trajectory(2,:))); % ymax
	boundingBox(5) = min(boundingBox(5),min(trajectory(3,:))); % zmin
	boundingBox(6) = max(boundingBox(6),max(trajectory(3,:))); % zmax
end