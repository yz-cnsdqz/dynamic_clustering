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

function skel_lines = createSkeletonLines(skel_ind,current_frame,varargin)

global VARS_GLOBAL_ANIM

if (~iscell(skel_ind))
    skel = VARS_GLOBAL_ANIM.skel(skel_ind);
    mot = VARS_GLOBAL_ANIM.mot(skel_ind);
else
    skel = skel_ind{1};
    mot = skel_ind{2};
end

color = VARS_GLOBAL_ANIM.animated_skeleton_Color;
linewidth = VARS_GLOBAL_ANIM.animated_skeleton_LineWidth;
linestyle = VARS_GLOBAL_ANIM.animated_skeleton_LineStyle;
marker = VARS_GLOBAL_ANIM.animated_skeleton_Marker;
markersize = VARS_GLOBAL_ANIM.animated_skeleton_MarkerSize;
markeredgecolor = VARS_GLOBAL_ANIM.animated_skeleton_MarkerEdgeColor;
markerfacecolor = VARS_GLOBAL_ANIM.animated_skeleton_MarkerFaceColor;
% color, linewidth, linestyle, marker, markersize, markeredgecolor, markerfacecolor
if (nargin>8)
    markerfacecolor = varargin{9-2};
end
if (nargin>7)
    markeredgecolor = varargin{8-2};
end
if (nargin>6)
    markersize = varargin{7-2};
end
if (nargin>5)
    marker = varargin{6-2};
end
if (nargin>4)
    linestyle = varargin{5-2};
end
if (nargin>3)
    linewidth = varargin{4-2};
end
if (nargin>2)
    color = varargin{3-2};
end

%%%%%%%%%%% clear skel_lines if necessary
npaths = size(skel.paths,1);

skel_lines = cell(npaths,1);

for k = 1:npaths
    path = skel.paths{k};
    nlines = length(path)-1;
    px = zeros(2,1); py = zeros(2,1); pz = zeros(2,1);
    px(1) = mot.jointTrajectories{path(1)}(1,current_frame); 
    py(1) = mot.jointTrajectories{path(1)}(2,current_frame); 
    pz(1) = mot.jointTrajectories{path(1)}(3,current_frame);
    for j = 2:nlines % path number
        px(2) = mot.jointTrajectories{path(j)}(1,current_frame); 
        py(2) = mot.jointTrajectories{path(j)}(2,current_frame); 
        pz(2) = mot.jointTrajectories{path(j)}(3,current_frame);
        skel_lines{k}(j-1) = line(px,py,pz,'Color',color,'LineWidth',linewidth,'linestyle',linestyle,'Parent',gca,'marker',marker,'markersize',markersize,'markeredgecolor',markeredgecolor,'markerfacecolor',markerfacecolor);
        px(1) = px(2);
        py(1) = py(2);
        pz(1) = pz(2);
    end
    px(2) = mot.jointTrajectories{path(nlines+1)}(1,current_frame); 
    py(2) = mot.jointTrajectories{path(nlines+1)}(2,current_frame); 
    pz(2) = mot.jointTrajectories{path(nlines+1)}(3,current_frame);
    skel_lines{k}(nlines) = line(px,py,pz,'Color',color,'LineWidth',linewidth,'linestyle',linestyle,'Parent',gca,'marker',marker,'markersize',markersize,'markeredgecolor',markeredgecolor,'markerfacecolor',markerfacecolor);
end
