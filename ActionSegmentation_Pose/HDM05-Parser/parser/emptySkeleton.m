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

function skel = emptySkeleton
skel = struct('njoints',0,...                           % number of joints
              'rootRotationalOffsetEuler',[0;0;0],...   % global (constant) rotation of root in world system, Euler angles.
              'rootRotationalOffsetQuat',[1;0;0;0],...  % global (constant) rotation of root in world system, quaternion.
              'nodes',struct([]),...                    % struct array containing nodes of skeleton tree data structure
              'paths',cell(1,1),...                     % contains a set of edge-disjoint paths the union of which represents the whole tree; represented as cell array of joint ID vectors
              'jointNames',cell(1,1),...                % cell array of joint names: maps joint ID to joint name
              'boneNames',cell(1,1),...                 % cell array of bone names: maps bone ID to node name. ID 1 is the root.
              'nameMap',cell(1,1),...                   % cell array mapping standard joint names to DOF IDs and trajectory IDs
              'animated',[],...                         % vector of IDs for animated joints/bones
              'unanimated',[],...                       % vector of IDs for unanimated joints/bones
              'filename','',...                         % source filename
              'version','',...                          % file format version
              'name','',...                             % name for this skeleton
              'massUnit',1,...                          % unit divisor for mass
              'lengthUnit',1,...                        % unit divisor for lengths
              'angleUnit','deg',...                     % angle unit (deg or rad)
              'documentation','',...                    % documentation from source file   
              'fileType','',...                         % file type (BVH / ASF)
              'skin',cell(1,1));                        % cell array of skin filenames for this skeleton