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

function mot = reflectMotYZ(mot)
%%%%%%%% reflect elderwalk about yz plane
for k=1:length(mot.jointTrajectories)
    mot.jointTrajectories{k,1}([1],:) = -mot.jointTrajectories{k,1}([1],:);
    mot.jointTrajectories{k,1} = mot.jointTrajectories{k,1} - repmat(mot.rootTranslation(:,1) - mot.rootTranslation(:,1),1,mot.nframes);
end
mot.rootTranslation([1],:) = -mot.rootTranslation([1],:);
mot.rootTranslation = mot.rootTranslation - repmat(mot.rootTranslation(:,1) - mot.rootTranslation(:,1),1,mot.nframes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
