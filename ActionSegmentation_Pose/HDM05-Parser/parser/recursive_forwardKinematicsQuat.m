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

function trajectories = recursive_forwardKinematicsQuat(skel, mot, node_id, current_position, current_rotation, trajectories)
% trajectories = recursive_forwardKinematicsQuat(skel, mot, node_id, current_position, current_rotation, trajectories)
%
% skel:                 skeleton
% mot:                  motion
% node_id:              index of current node in node array
% current_position:     for all frames: current local coordinate offsets from world origin (3xnframes sequence of vectors)
% current_rotation:     for all frames: current local coordinate rotations against world coordinate frame (4xnframes sequence of quaternions)
% trajectories:         input & output cell array containing the trajectories that have been computed so far

trajectories{node_id,1} = current_position;

for child_id = skel.nodes(node_id).children'
    child = skel.nodes(child_id);

    if (~isempty(mot.rotationQuat{child_id}))
        child_rotation = quatmult(current_rotation,mot.rotationQuat{child_id});
    else
        child_rotation = current_rotation;
    end

    child_position = current_position + quatrot(repmat(child.offset,1,mot.nframes),child_rotation);

    trajectories = recursive_forwardKinematicsQuat(skel, mot, child_id, child_position, child_rotation, trajectories);
end