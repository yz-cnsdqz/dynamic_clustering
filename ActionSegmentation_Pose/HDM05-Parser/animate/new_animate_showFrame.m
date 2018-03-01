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

function new_animate_showFrame(obj,event,varargin)

global VARS_GLOBAL_ANIM

tic;

if (VARS_GLOBAL_ANIM.animation_done)
    toc;
    return;
end


if (nargin>3) % direct frame draw mode!
    current_frame = varargin{2};
    VARS_GLOBAL_ANIM.current_frame = current_frame;
    frames_total = inf;
else
    current_frame = VARS_GLOBAL_ANIM.current_frame;
    frames_total = VARS_GLOBAL_ANIM.frames_total;
end

if (nargin>2)
    observer_fcn = varargin{1};
else
    observer_fcn = {};
end

if (current_frame>=frames_total)
    t = timerfind('Name','AnimationTimer');
    if (~isempty(t))
        if(VARS_GLOBAL_ANIM.kill_timer)
            delete(t);
            VARS_GLOBAL_ANIM.animation_done = true;
            avg_frame_time = VARS_GLOBAL_ANIM.frame_draw_time / VARS_GLOBAL_ANIM.frames_drawn;
            actual_frame_rate = 1/avg_frame_time;
            disp(['Frame rate: ' num2str(actual_frame_rate) ' fps.']);
            disp(['Total frame drawing time: ' num2str(VARS_GLOBAL_ANIM.frame_draw_time) ' sec.']);
        else
            stop(t);
            VARS_GLOBAL_ANIM.animation_paused = true;
        end
    end
    if(VARS_GLOBAL_ANIM.loop_playback)
        current_frame = 1;
        VARS_GLOBAL_ANIM.current_frame = 1;
        VARS_GLOBAL_ANIM.animation_paused = false;
		set(t,'TasksToExecute', VARS_GLOBAL_ANIM.frames_total+1);
        start(t);
    end
end

%if (current_frame>frames_total)
%	disp('Error: Timer sent too many calls to animate_showFrame!');
%	return;
%end

if(~isempty(observer_fcn))
    feval(observer_fcn{1},VARS_GLOBAL_ANIM.current_frame,observer_fcn{2:end});
end


% drop frames if more than 30 fps
elapsed = etime(clock,VARS_GLOBAL_ANIM.previous_call);
if (elapsed < 1/30 && current_frame ~= frames_total)
    VARS_GLOBAL_ANIM.current_frame = current_frame + 1;
    t2 = toc;
    VARS_GLOBAL_ANIM.frame_draw_time = VARS_GLOBAL_ANIM.frame_draw_time + t2;
    return;
end


%%%%%%%%%%%%% draw "current_frame" for all skeletons
for i = 1:length(VARS_GLOBAL_ANIM.skel)
    if (current_frame > length(VARS_GLOBAL_ANIM.range{i})) % if motion number i need not be animated anymore, skip it!
        continue;
    end
    
    % delay beginning of motion for a short period of time
    %current_frame = max(1,current_frame - 19);
    
    
    current_frame_in_mot = VARS_GLOBAL_ANIM.range{i}(current_frame);
    
    set(VARS_GLOBAL_ANIM.graphics_data.frameLabel, 'String', ['frame ' int2str(current_frame_in_mot)]);
    
    npaths = size(VARS_GLOBAL_ANIM.skel(i).paths,1);
    for k = 1:npaths
        path = VARS_GLOBAL_ANIM.skel(i).paths{k};
        nlines = length(path)-1;
        px = zeros(2,1); py = zeros(2,1); pz = zeros(2,1);
        px(1) = VARS_GLOBAL_ANIM.mot(i).jointTrajectories{path(1)}(1,current_frame_in_mot); 
        py(1) = VARS_GLOBAL_ANIM.mot(i).jointTrajectories{path(1)}(2,current_frame_in_mot); 
        pz(1) = VARS_GLOBAL_ANIM.mot(i).jointTrajectories{path(1)}(3,current_frame_in_mot);
        for j = 2:nlines % path number
            px(2) = VARS_GLOBAL_ANIM.mot(i).jointTrajectories{path(j)}(1,current_frame_in_mot); 
            py(2) = VARS_GLOBAL_ANIM.mot(i).jointTrajectories{path(j)}(2,current_frame_in_mot); 
            pz(2) = VARS_GLOBAL_ANIM.mot(i).jointTrajectories{path(j)}(3,current_frame_in_mot);
            set(VARS_GLOBAL_ANIM.graphics_data(VARS_GLOBAL_ANIM.graphics_data_index).skel_lines{i}{k}(j-1),'XData',px,'YData',py,'ZData',pz);
            
            px(1) = px(2);
            py(1) = py(2);
            pz(1) = pz(2);
        end
        px(2) = VARS_GLOBAL_ANIM.mot(i).jointTrajectories{path(nlines+1)}(1,current_frame_in_mot); 
        py(2) = VARS_GLOBAL_ANIM.mot(i).jointTrajectories{path(nlines+1)}(2,current_frame_in_mot); 
        pz(2) = VARS_GLOBAL_ANIM.mot(i).jointTrajectories{path(nlines+1)}(3,current_frame_in_mot);
        set(VARS_GLOBAL_ANIM.graphics_data(VARS_GLOBAL_ANIM.graphics_data_index).skel_lines{i}{k}(nlines),'XData',px,'YData',py,'ZData',pz);
    end
    
    if (isfield(VARS_GLOBAL_ANIM.mot(i),'animated_patch_data'))
        npatches = length(VARS_GLOBAL_ANIM.mot(i).animated_patch_data);
        for k=1:npatches
            X = VARS_GLOBAL_ANIM.mot(i).animated_patch_data(k).X(:,current_frame_in_mot);
            Y = VARS_GLOBAL_ANIM.mot(i).animated_patch_data(k).Y(:,current_frame_in_mot);
            Z = VARS_GLOBAL_ANIM.mot(i).animated_patch_data(k).Z(:,current_frame_in_mot);
            if (size(VARS_GLOBAL_ANIM.mot(i).animated_patch_data(k).color,1)==1)
                %C = repmat(VARS_GLOBAL_ANIM.mot(i).animated_patch_data(k).color,size(Z,1),1);
				C = VARS_GLOBAL_ANIM.mot(i).animated_patch_data(k).color(ones(size(Z,1),1),:);
            else
                C = repmat(VARS_GLOBAL_ANIM.mot(i).animated_patch_data(k).color(current_frame_in_mot,:),size(Z,1),1);
            end
            
            switch (lower(VARS_GLOBAL_ANIM.mot(i).animated_patch_data(k).type))
                case 'disc'
                    set(VARS_GLOBAL_ANIM.graphics_data(VARS_GLOBAL_ANIM.graphics_data_index).animated_patches{i}(k),...
                        'Vertices',[X Y Z],...
                        'FaceVertexCData',C);
                case 'polygondisc'
                    set(VARS_GLOBAL_ANIM.graphics_data(VARS_GLOBAL_ANIM.graphics_data_index).animated_patches{i}(k),...
                        'Vertices',[X Y Z],...
                        'FaceVertexCData',C);
                case 'griddisc'
                    set(VARS_GLOBAL_ANIM.graphics_data(VARS_GLOBAL_ANIM.graphics_data_index).animated_patches{i}(k),...
                        'Vertices',[X Y Z],...
                        'FaceVertexCData',C);
                case 'point'
                    if isempty(VARS_GLOBAL_ANIM.animated_point_MarkerEdgeColor)
                        markeredgecolor = C;
                    else
                        markeredgecolor = repmat(VARS_GLOBAL_ANIM.animated_point_MarkerEdgeColor,size(Z,1),1);
                    end
                    set(VARS_GLOBAL_ANIM.graphics_data(VARS_GLOBAL_ANIM.graphics_data_index).animated_patches{i}(k),...
                        'Vertices',[X Y Z],...
                        'MarkerEdgeColor',markeredgecolor,...
                        'MarkerFaceColor',C);
                case 'cappedcylinder'
                    set(VARS_GLOBAL_ANIM.graphics_data(VARS_GLOBAL_ANIM.graphics_data_index).animated_patches{i}(k),...
                        'Vertices',[X Y Z],...
                        'FaceVertexCData',C);
            end
            
        end
    end
    
    if (~isempty(VARS_GLOBAL_ANIM.draw_labels) && VARS_GLOBAL_ANIM.draw_labels(i))
        set(VARS_GLOBAL_ANIM.graphics_data(VARS_GLOBAL_ANIM.graphics_data_index).text_handle(i),'string',sprintf('              Frame %d',current_frame),'position',VARS_GLOBAL_ANIM.mot(i).jointTrajectories{1}(:,current_frame_in_mot));
    end
end

drawnow;
t2 = toc;

VARS_GLOBAL_ANIM.frames_drawn = VARS_GLOBAL_ANIM.frames_drawn + 1;
VARS_GLOBAL_ANIM.frame_draw_time = VARS_GLOBAL_ANIM.frame_draw_time + t2;

VARS_GLOBAL_ANIM.current_frame = current_frame + 1;

VARS_GLOBAL_ANIM.previous_call = clock;