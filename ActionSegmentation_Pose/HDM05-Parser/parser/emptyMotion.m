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

function mot = emptyMotion

        mot = struct('njoints',0,...                            % number of joints
                     'nframes',0,...                            % number of frames
                     'frameTime',1/120,...                      % inverse sampling rate: time per frame (in seconds)
                     'samplingRate',120,...                     % sampling rate (in Hertz) (120 Hertz is Carnegie-Mellon Mocap-DB standard)
                     'jointTrajectories',cell(1,1),...          % 3D joint trajectories
                     'rootTranslation',[],...                   % global translation data stream of the root
                     'rotationEuler',cell(1,1),...              % rotational data streams for all joints, including absolute root rotation at pos. 1, Euler angles
                     'rotationQuat',cell(1,1),...               % rotational data streams for all joints, including absolute root rotation at pos. 1, quaternions
                     'jointNames',cell(1,1),...                 % cell array of joint names: maps node ID to joint name
                     'boneNames',cell(1,1),...                  % cell array of bone names: maps bone ID to node name. ID 1 is the root.
                     'nameMap',cell(1,1),...                    % cell array mapping standard joint names to DOF IDs and trajectory IDs
                     'animated',[],...                          % vector of IDs for animated joints/bones
                     'unanimated',[],...                        % vector of IDs for unanimated joints/bones
                     'boundingBox',[],...                       % bounding box (given a specific skeleton)
                     'filename','',...                          % source filename
                     'documentation','',...                     % documentation from source file   
                     'angleUnit','deg');                        % angle unit, either deg or rad
