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

function mot_out = resampleMot_new(skel_in, mot_in,target_frame_rate)

target_frame_rate = round(target_frame_rate);

mot_out = emptyMotion;

mot_out.njoints = mot_in.njoints;
mot_out.frameTime = 1/target_frame_rate;
mot_out.samplingRate = target_frame_rate;
mot_out.jointNames = mot_in.jointNames;
mot_out.boneNames = mot_in.boneNames;
mot_out.nameMap = mot_in.nameMap;
mot_out.animated = mot_in.animated;
mot_out.unanimated = mot_in.unanimated;

if ~isempty(mot_in.rootTranslation)
    %mot_out.rootTranslation = resample(mot_in.rootTranslation',target_frame_rate,round(mot_in.samplingRate))';
    mot_out.rootTranslation = resampleTSData(mot_in.rootTranslation,round(mot_in.samplingRate),target_frame_rate);
end

if ~isempty(mot_in.rotationQuat)
    
    for k = 1:size(mot_in.rotationQuat,1)
        if (~isempty(mot_in.rotationQuat{k,1}))
            %mot_out.rotationQuat{k,1} = resample(mot_in.rotationQuat{k,1}',target_frame_rate,round(mot_in.samplingRate))';
            mot_out.rotationQuat{k,1} = resampleTSData(mot_in.rotationQuat{k,1},round(mot_in.samplingRate),target_frame_rate);
            mot_out.rotationQuat{k,1} = quatnormalize(mot_out.rotationQuat{k,1});
        end
    end
    
    mot_out.nframes= size(mot_out.rotationQuat{1}, 2);
    
    mot_out.jointTrajectories = forwardKinematicsQuat(skel_in, mot_out);
    for k = 1:size(mot_in.rotationEuler,1)
        if (~isempty(mot_in.rotationEuler{k}))
            mot_out.rotationEuler{k,1} = quat2euler(mot_in.rotationQuat{k, 1});
        end
    end
    
else
    error('This motion does not contain the field "rotationQuat" or the field is empty!\nIf you want to resize motions containing C3D marker positions, better use "resampleMotion"');
end



mot_out.nframes = max( size(mot_out.rootTranslation,2),size(mot_out.jointTrajectories{1},2));

mot_out.filename = mot_in.filename;
mot_out.documentation = vertcat(mot_in.documentation,{['The original file has been resampled to new frame rate ' num2str(target_frame_rate) ' Hz!']});
mot_out.angleUnit = mot_in.angleUnit;
mot_out.boundingBox = mot_in.boundingBox;