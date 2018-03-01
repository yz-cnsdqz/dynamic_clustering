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

function clearTracePoses
global VARS_GLOBAL_ANIM
if (isempty(VARS_GLOBAL_ANIM)||isempty(VARS_GLOBAL_ANIM.graphics_data))
    return;
end

f = find(cell2mat({VARS_GLOBAL_ANIM.graphics_data.figure}) == gcf);
if (isempty(f))
    return;
end
VARS_GLOBAL_ANIM.graphics_data_index = f;

num_trace_poses = size(VARS_GLOBAL_ANIM.graphics_data(f).trace_poses,2);
num_lines_per_pose = size(VARS_GLOBAL_ANIM.graphics_data(f).trace_poses,1);
for k=1:num_trace_poses
    for j=1:num_lines_per_pose
        delete(VARS_GLOBAL_ANIM.graphics_data(f).trace_poses{j,k}(ishandle(VARS_GLOBAL_ANIM.graphics_data(f).trace_poses{j,k})));
    end
end
VARS_GLOBAL_ANIM.graphics_data(f).trace_poses = {};