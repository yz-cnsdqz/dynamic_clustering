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

function trajectories = recursive_forwardKinematicsEuler(skel, mot, node_id, current_position, current_rotation, trajectories)
% trajectories = recursive_forwardKinematicsEuler(skel, mot, node_id, current_position, current_rotation, trajectories)
%
% skel:                 skeleton
% mot:                  motion
% node_id:              index of current node in node array
% current_position:     for all frames: current local coordinate offsets from world origin (3xnframes sequence of vectors)
% current_rotation:     for all frames: current local coordinate rotations against world coordinate frame (4xnframes sequence of quaternions)
% trajectories:         input & output cell array containing the trajectories that have been computed so far

switch lower(mot.angleUnit)
    case 'deg'
        conversion_factor = pi/180;
    case 'rad'
        conversion_factor = 1;
    otherwise
        error(['Unknown angle unit: ' mot.angleUnit]);
end

trajectories{node_id,1} = current_position;

for child_id = skel.nodes(node_id).children'
    child = skel.nodes(child_id);

    if (~isempty(mot.rotationEuler{child_id}))
        nRotDOF = size(mot.rotationEuler{child.ID,1},1); % number of rotational DOFs
        
        % zero-pad Euler array at the appropriate places, according to rotation order and presence/absence of respective DOFs
        if nRotDOF > 0
            completeEulers = zeros(3,mot.nframes);
            d = 1; % index for DOFs present in mot.rotationEuler
            for r = 1:3 % go through rotationOrder. 
            	idx = strmatch(['r' lower(child.rotationOrder(r))], lower(child.DOF), 'exact');
                if (~isempty(idx))
                    completeEulers(r,:) = mot.rotationEuler{child.ID,1}(d,:)*conversion_factor;
                    d = d+1;
                end
            end
            axis_quat = repmat(euler2quat(flipud(child.axis)*conversion_factor,'zyx'),1,mot.nframes); % According to ASF specs, rotation order for "axis" should be XYZ. However, they use the opposite multiplication order as we do!
            rotationQuat = quatmult(axis_quat,quatmult(euler2quat(flipud(completeEulers),fliplr(child.rotationOrder)),quatinv(axis_quat))); % ASF specs use opposite multiplication order as we do, hence fliplr() and flipud()!
        end
        child_rotation = quatmult(current_rotation,rotationQuat);
    else
        child_rotation = current_rotation;
    end

    child_position = current_position + quatrot(repmat(child.offset,1,mot.nframes),child_rotation);

    trajectories = recursive_forwardKinematicsEuler(skel, mot, child_id, child_position, child_rotation, trajectories);
end