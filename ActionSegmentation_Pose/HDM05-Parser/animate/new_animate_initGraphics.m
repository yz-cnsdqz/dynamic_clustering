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

function new_animate_initGraphics

global VARS_GLOBAL_ANIM
 
VARS_GLOBAL_ANIM.loop_playback = 1;

bb = [inf -inf inf -inf inf -inf]';
ground_tile_size = -inf;
for k = 1:length(VARS_GLOBAL_ANIM.mot)
    bb([2 4 6]) = max(bb([2 4 6]), VARS_GLOBAL_ANIM.mot(k).boundingBox([2 4 6]));
    bb([1 3 5]) = min(bb([1 3 5]), VARS_GLOBAL_ANIM.mot(k).boundingBox([1 3 5]));

    lclavicleIdx = strMatch('lshoulder', {VARS_GLOBAL_ANIM.skel(1).nameMap{:,1}}, 'exact');
    if(isempty(lclavicleIdx))
        lclavicleIdx = strMatch('LSHO', {VARS_GLOBAL_ANIM.skel(1).nameMap{:,1}}, 'exact');
        if(isempty(lclavicleIdx))
            lclavicleIdx = strMatch('lsho', {VARS_GLOBAL_ANIM.skel(1).nameMap{:,1}}, 'exact');
        end
    end
    
    
    if isempty(lclavicleIdx) 
        warning('Could not find the clavicle in the nameMap. Using default size 10 [?] of the ground plane grid.');
        x=10;
    else
        lclavicleIdx = VARS_GLOBAL_ANIM.skel(1).nameMap{lclavicleIdx, 3};
            
        lelbowIdx = strMatch('lelbow', {VARS_GLOBAL_ANIM.skel(1).nameMap{:,1}});
        if(isempty(lelbowIdx))
            lelbowIdx = strMatch('LELB', {VARS_GLOBAL_ANIM.skel(1).nameMap{:,1}}, 'exact');
            if(isempty(lelbowIdx))
                lelbowIdx = strMatch('lelb', {VARS_GLOBAL_ANIM.skel(1).nameMap{:,1}}, 'exact');
            end
        end
        lelbowIdx = VARS_GLOBAL_ANIM.skel(1).nameMap{lelbowIdx, 3};
        x = VARS_GLOBAL_ANIM.ground_tile_size_factor*sqrt(sum((VARS_GLOBAL_ANIM.mot(1).jointTrajectories{lclavicleIdx}(:,1)-VARS_GLOBAL_ANIM.mot(1).jointTrajectories{lelbowIdx}(:,1)).^2));
    end

%     x = VARS_GLOBAL_ANIM.ground_tile_size_factor*VARS_GLOBAL_ANIM.skel(1).nodes(VARS_GLOBAL_ANIM.skel(1).nameMap{16,2}).length; % ground tile size will be a multiple of the length of longest lhumerus (DOF id 16)
    
    if (x>ground_tile_size)
        ground_tile_size = x;
    end
end
if (sum(isinf(bb)>0))
    ground_tile_size = 1;
    bb_extension = zeros(6,1);
    bb = [-1 1 0 1 -1 1]';
else
    bb_extension = [-1 1 0 1 -1 1]'*VARS_GLOBAL_ANIM.bounding_box_border_extension*ground_tile_size;
    bb = bb + bb_extension;
end
%nground_tiles = 10;
%ground_tile_size = min(bb(2) - bb(1), bb(6) - bb(5)) / nground_tiles;

create_new = false;
if (~isempty(VARS_GLOBAL_ANIM.graphics_data))
    % eliminate invalid figure handles
    valid_figure_handle_indices = find(ishandle(cell2mat({VARS_GLOBAL_ANIM.graphics_data.figure})));
    VARS_GLOBAL_ANIM.graphics_data = VARS_GLOBAL_ANIM.graphics_data(valid_figure_handle_indices);

    view_curr = view(gca);
%    newplot;
    current_figure = gcf;
    
    % find current figure in graphics_data
    VARS_GLOBAL_ANIM.graphics_data_index = find(cell2mat({VARS_GLOBAL_ANIM.graphics_data.figure}) == current_figure);
    if (isempty(VARS_GLOBAL_ANIM.graphics_data_index)) % current figure hasn't been used for our skeleton graphics output previously
        create_new = true;
        VARS_GLOBAL_ANIM.graphics_data_index = length(VARS_GLOBAL_ANIM.graphics_data)+1; % create new entry in graphics_data array
        VARS_GLOBAL_ANIM.graphics_data(VARS_GLOBAL_ANIM.graphics_data_index) = struct('figure',current_figure,...
                                                                                      'skel_lines',cell(1,1),...
                                                                                      'animated_patches',cell(1,1),...
                                                                                      'ground_plane',0,...
                                                                                      'text_handle',0,...
                                                                                      'view',[],...
                                                                                      'trace_poses',[]);
                                                                                      
    end
else % graphics_data IS empty!
    create_new = true;
    
    view_curr = view(gca);
%    newplot;
    current_figure = gcf;

    VARS_GLOBAL_ANIM.graphics_data = struct('figure',current_figure,...
                                            'skel_lines',cell(1,1),...
                                            'animated_patches',cell(1,1),...
                                            'ground_plane',0,...
                                            'text_handle',0,...
                                            'view',[],...
                                            'trace_poses',[]);
    VARS_GLOBAL_ANIM.graphics_data_index = 1;
end

i = VARS_GLOBAL_ANIM.graphics_data_index;

if (~create_new)
    delete(VARS_GLOBAL_ANIM.graphics_data(i).ground_plane(find(ishandle(VARS_GLOBAL_ANIM.graphics_data(i).ground_plane))));
    VARS_GLOBAL_ANIM.graphics_data(i).ground_plane = [];
    delete(VARS_GLOBAL_ANIM.graphics_data(i).text_handle(find(ishandle(VARS_GLOBAL_ANIM.graphics_data(i).text_handle) & (VARS_GLOBAL_ANIM.graphics_data(i).text_handle > 0))));
    VARS_GLOBAL_ANIM.graphics_data(i).text_handle = [];
    for k=1:length(VARS_GLOBAL_ANIM.graphics_data(i).skel_lines)
        if (~isempty(VARS_GLOBAL_ANIM.graphics_data(i).skel_lines{k}))
            for j = 1:size(VARS_GLOBAL_ANIM.graphics_data(i).skel_lines{k},1)
                delete(VARS_GLOBAL_ANIM.graphics_data(i).skel_lines{k}{j}(find(ishandle(VARS_GLOBAL_ANIM.graphics_data(i).skel_lines{k}{j}))));
                VARS_GLOBAL_ANIM.graphics_data(i).skel_lines{k}{j} = [];
            end
    	end
    end
    for k=1:length(VARS_GLOBAL_ANIM.graphics_data(i).animated_patches)
        if (~isempty(VARS_GLOBAL_ANIM.graphics_data(i).animated_patches{k}))
            delete(VARS_GLOBAL_ANIM.graphics_data(i).animated_patches{k}(find(ishandle(VARS_GLOBAL_ANIM.graphics_data(i).animated_patches{k}))));
            VARS_GLOBAL_ANIM.graphics_data(i).animated_patches{k} = [];
    	end
    end
end
reset(gca);
if (isempty(VARS_GLOBAL_ANIM.graphics_data(i).view))
    VARS_GLOBAL_ANIM.graphics_data(i).view = [-0.9129    0.0000    1.7180   -0.4026;...
                                               0.0461    1.6227    0.4335   -1.0511;...
                                               0.4056   -0.1843    3.8170   21.0978;...
                                                    0         0         0    1.0000];
else
    VARS_GLOBAL_ANIM.graphics_data(i).view = view_curr;
end

hld=ishold;
set(gcf,'NextPlot','add');
set(gca,'NextPlot','add');
%%%%%%%%%%% create ground plane
%bb([1 2 5 6],:) = bb([1 2 5 6],:)*10;
%ground_tile_size=ground_tile_size*10;
VARS_GLOBAL_ANIM.graphics_data(i).ground_plane = createGroundPlane(bb,bb_extension,ground_tile_size);

for k=1:length(VARS_GLOBAL_ANIM.skel)
%%%%%%%%%%% create text fields for frame number
    if (~isempty(VARS_GLOBAL_ANIM.draw_labels) && VARS_GLOBAL_ANIM.draw_labels(k))
        VARS_GLOBAL_ANIM.graphics_data(i).text_handle(k) = text(0,0,0,'');
        set(VARS_GLOBAL_ANIM.graphics_data(i).text_handle(k),'interpreter','none');
    end
%%%%%%%%%%% create lines for skeleton visualization
    range = VARS_GLOBAL_ANIM.range{k};
    current_frame = range(1);
    VARS_GLOBAL_ANIM.graphics_data(i).skel_lines{k,1} = createSkeletonLines(k,current_frame);
    
%%%%%%%%%%% create patches for animated plane/cylinder/sphere/etc visualization
    if (isfield(VARS_GLOBAL_ANIM.mot(k),'animated_patch_data'))
        VARS_GLOBAL_ANIM.graphics_data(i).animated_patches{k,1} = createAnimatedPatches(k,current_frame);
    end
end
if (~hld)
    %set(gcf,'NextPlot','replace');
    set(gca,'NextPlot','replace');
end

figure(VARS_GLOBAL_ANIM.graphics_data(i).figure);
set(gcf,'renderer','OpenGL','DoubleBuffer','on');
%set(gcf,'Color','w');
set(gcf,'Color',VARS_GLOBAL_ANIM.figure_color);
% set(gca,'Color','w',...
%         'XColor','w',...
%         'YColor','w',...
%         'ZColor','w');
set(gca,'visible','off');

axis equal;
axis(bb+[-1 1 -1 1 -1 1]'*ground_tile_size);
camup([0 1 0]);
%view(VARS_GLOBAL_ANIM.graphics_data(i).view);
camproj('perspective');
camlight headlight;
light;

%set(gcf,'position',[5 5 512 284]);
%set(gca,'Units','centimeter')
%camva('auto');
%camva('manual');
