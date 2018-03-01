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

function showTracePoses(skel,mot,frames)

global VARS_GLOBAL_ANIM

if (isempty(VARS_GLOBAL_ANIM)||isempty(VARS_GLOBAL_ANIM.graphics_data))
    VARS_GLOBAL_ANIM = emptyVarsGlobalAnimStruct;
    animate_initGraphics;
end

f = find(cell2mat({VARS_GLOBAL_ANIM.graphics_data.figure}) == gcf);
if (isempty(f))
    animate_initGraphics;
    f = gcf;
end

VARS_GLOBAL_ANIM.graphics_data_index = f;
for k=1:length(frames)
    VARS_GLOBAL_ANIM.graphics_data(f).trace_poses=...
        [VARS_GLOBAL_ANIM.graphics_data(f).trace_poses; ...
            createSkeletonLines({skel,mot},frames(k),...
                                                    VARS_GLOBAL_ANIM.trace_pose_Color,...
                                                    VARS_GLOBAL_ANIM.trace_pose_LineWidth,...
                                                    VARS_GLOBAL_ANIM.trace_pose_LineStyle,...
                                                    VARS_GLOBAL_ANIM.trace_pose_Marker,...
                                                    VARS_GLOBAL_ANIM.trace_pose_MarkerSize,...
                                                    VARS_GLOBAL_ANIM.trace_pose_MarkerEdgeColor,...
                                                    VARS_GLOBAL_ANIM.trace_pose_MarkerFaceColor)];
% {skel, mot}, frames, color, linewidth, linestyle, marker, markersize, markeredgecolor, markerfacecolor
end