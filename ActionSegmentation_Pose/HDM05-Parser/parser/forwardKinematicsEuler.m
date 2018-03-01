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

function jointTrajectories = forwardKinematicsEuler(skel,mot)
switch lower(mot.angleUnit)
    case 'deg'
        conversion_factor = pi/180;
    case 'rad'
        conversion_factor = 1;
    otherwise
        error(['Unknown angle unit: ' mot.angleUnit]);
end

% root node involves special case for determination of rotation order (node.rotationOrder only concerns the global rotational offset in this case)
rootTransformationOrder = char(skel.nodes(1).DOF);
rootTransformationOrder = rootTransformationOrder(:,2)';
[x,IA,IB] = intersect({'rx','ry','rz'}, lower(skel.nodes(1).DOF));
rootRotationOrder = rootTransformationOrder(sort(IB));
mot.rotationQuat{node.ID,1} = euler2quat(flipud(completeEulers),fliplr(node.rotationOrder)); % ASF specs use opposite multiplication order as we do, hence fliplr() and flipud()!

jointTrajectories = recursive_forwardKinematicsEuler(skel,...
                                                mot,...
                                                1,...
                                                mot.rootTranslation + repmat(skel.nodes(1).offset,1,mot.nframes),...
                                                quatmult(repmat(skel.rootRotationalOffsetQuat,1,mot.nframes),euler2quat(conversion_factor*flipud(mot.rotationEuler{1}),rootRotationOrder)),...
                                                mot.jointTrajectories);