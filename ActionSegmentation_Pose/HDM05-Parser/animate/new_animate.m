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

function new_animate(skel,mot,varargin)
% new_animate(skel,mot,num_repeats,time_stretch_factor,range,draw_labels,observer_fcn,start_animation)
% INPUT:  
%   - skel and mot are assumed to be struct arrays of identical length containing 
%     skeleton/motion pairs to be animated simultaneously. All mots are supposed to have the same samplingRate.
%   - range is a cell array which is supposed to be of same length as skel/mot. 
%     Empty cell array entries in range indicate that the entire frame range is to be played back.
%   - draw_labels is a logical vector the same length as skel indicating whether the corresponding skeleton 
%     is to be drawn with a frame counter label. Empty draw_labels means "draw no labels at all".
%	- observer_fcn function called at each frame animation; used for progress monitoring
%	- start_animation determines whether the actual animation shall start immediately or pause until resumeAnimation is called.  

global VARS_GLOBAL_ANIM
if (isempty(VARS_GLOBAL_ANIM))
    VARS_GLOBAL_ANIM = emptyVarsGlobalAnimStruct;
end

num_repeats = 1;
time_stretch_factor = 1;
VARS_GLOBAL_ANIM.range = cell(length(mot),1);
VARS_GLOBAL_ANIM.draw_labels = ones(length(mot),1);
start_animation = true;

timer_fcn = {'new_animate_showFrame'};

if (nargin>7)
	start_animation = varargin{6};
end

if (nargin>6)
	timer_fcn = {'new_animate_showFrame',varargin{5}};
end

if (nargin > 5)
    VARS_GLOBAL_ANIM.draw_labels = varargin{4};
end
if (nargin > 4)
    VARS_GLOBAL_ANIM.range = varargin{3};
    if (~iscell(VARS_GLOBAL_ANIM.range))
        VARS_GLOBAL_ANIM.range = {VARS_GLOBAL_ANIM.range};
    end
end
if (nargin > 3)
    time_stretch_factor = varargin{2};
end
if (nargin > 2)
    num_repeats = varargin{1};
end

VARS_GLOBAL_ANIM.skel = skel;
VARS_GLOBAL_ANIM.mot = mot;

num_frames = -inf;
if (isempty(VARS_GLOBAL_ANIM.range))
    VARS_GLOBAL_ANIM.range = cell(1,length(mot));
end
for k=1:length(VARS_GLOBAL_ANIM.mot)
    if (isempty(VARS_GLOBAL_ANIM.range{k}))
        VARS_GLOBAL_ANIM.range{k} = [1:VARS_GLOBAL_ANIM.mot(k).nframes];
    end
    if (length(VARS_GLOBAL_ANIM.range{k})>num_frames) % determine num_frames as maximum of range lengths
        num_frames = length(VARS_GLOBAL_ANIM.range{k});
    end
end

new_animate_initGraphics;

desired_frame_time = VARS_GLOBAL_ANIM.mot(1).frameTime;

if (~isempty(VARS_GLOBAL_ANIM.figure_camera_file))
    h = str2func(VARS_GLOBAL_ANIM.figure_camera_file);
    feval(h, gca);
end

if (~isempty(VARS_GLOBAL_ANIM.figure_position))
    set(gcf,'position',VARS_GLOBAL_ANIM.figure_position);
end
for (i=1:num_repeats)
    VARS_GLOBAL_ANIM.frame_draw_time = 0;
    VARS_GLOBAL_ANIM.frames_drawn = 0;
    VARS_GLOBAL_ANIM.animation_done = false;
    t=timerfind('Name','AnimationTimer');
    if (~isempty(t))
        delete(t);
    end
	t = timer('Name','AnimationTimer');
    try
        VARS_GLOBAL_ANIM.previous_call = clock;
		VARS_GLOBAL_ANIM.current_frame = 1;
        VARS_GLOBAL_ANIM.frames_total = num_frames;
		
        set(t,'TimerFcn',timer_fcn);
        set(t,'ExecutionMode','fixedRate');
		set(t,'TasksToExecute',num_frames);
        period = round(1000*desired_frame_time/time_stretch_factor)/1000;
        VARS_GLOBAL_ANIM.timer_period = period;
        if (period == 0)
            warning(['Requested timer period of ' num2str(desired_frame_time/time_stretch_factor) ' s is too short for Matlab! Setting period to minimum of 1 ms :-(']);
            period = 0.001;
        end
		set(t,'Period',period);
		set(t,'BusyMode','queue');
        if(start_animation)
		    start(t);
        end
    catch
		delete(t);
        return
    end
    
end
