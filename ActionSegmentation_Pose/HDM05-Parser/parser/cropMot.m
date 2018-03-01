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

function mot_out = cropMot(mot_in,range)

if isempty(range)
    mot_out = mot_in;
    return;
end

mot_out = emptyMotion;

mot_out.njoints = mot_in.njoints;
mot_out.nframes = length(range);
mot_out.frameTime = mot_in.frameTime;
mot_out.samplingRate = mot_in.samplingRate;
mot_out.jointNames = mot_in.jointNames;
mot_out.boneNames = mot_in.boneNames;
mot_out.nameMap = mot_in.nameMap;
mot_out.animated = mot_in.animated;
mot_out.unanimated = mot_in.unanimated;

if not(isempty(mot_in.rootTranslation))
    mot_out.rootTranslation = mot_in.rootTranslation(:,range);
end
for k = 1:size(mot_in.jointTrajectories,1)
    if ~isempty(mot_in.jointTrajectories{k,1})
        mot_out.jointTrajectories{k,1} = mot_in.jointTrajectories{k,1}(:,range);
    end
end
for k = 1:size(mot_in.rotationEuler,1)
    if (~isempty(mot_in.rotationEuler{k}))
        mot_out.rotationEuler{k,1} = mot_in.rotationEuler{k,1}(:,range);
    end
end
for k = 1:size(mot_in.rotationQuat,1)
    if (~isempty(mot_in.rotationQuat{k,1}))
        mot_out.rotationQuat{k,1} = mot_in.rotationQuat{k,1}(:,range);
    end
end
    
mot_out.filename = mot_in.filename;
mot_out.documentation = vertcat(mot_in.documentation,{['The original file has been cropped to the range ' num2str(range) '!']});
mot_out.angleUnit = mot_in.angleUnit;
mot_out.boundingBox = computeBoundingBox(mot_out);