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

function node = emptySkeletonNode

   node = struct('children',[],...            % struct array containing this node's child nodes' indices within "nodes" array of skeleton data structure
                 'jointName','',...           % name of this joint (if it is a joint (BVH))
                 'boneName','',...            % name of this bone (if it is a bone (ASF))
                 'ID',0,...                   % ID number of this joint/bone
                 'parentID',0,...             % parent node ID
                 'offset',[0;0;0],...         % offset vector from the parent of this node to this node
                 'direction',[0;0;0],...      % direction vector for this bone, given in global coordinates
                 'length',0,...               % length of this bone
                 'axis',[0;0;0],...           % axis of rotation for this node
                 'DOF',cell(1,1),...          % degrees of freedom and rotation order for this node, in the form {'tx','ty','tz','rx','ry','rz'} (e.g.)
                 'rotationOrder','',...       % rotation order in string representation (e.g., 'xyz')
                 'limits',[]);                % kx2, kx2 or kx2 matrix of limits for this node's DOFs
