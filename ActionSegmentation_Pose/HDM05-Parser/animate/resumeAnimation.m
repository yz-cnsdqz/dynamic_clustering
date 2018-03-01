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

function resumeAnimation(varargin)
global VARS_GLOBAL_ANIM

if (nargin >1)
    new_time_stretch_factor = varargin{2};
else
    new_time_stretch_factor = -1;
end

if (nargin >0)
    range = varargin{1};
    current_frame = range(1);
    if(length(range) == 2)
        last_frame = range(2);
    else
        last_frame = VARS_GLOBAL_ANIM.frames_total;
    end
else
    current_frame = VARS_GLOBAL_ANIM.current_frame;
    last_frame = VARS_GLOBAL_ANIM.frames_total;
end

t = timerfind('Name','AnimationTimer');
try
    if(~isempty(t))
        t = t(1);
        num_frames = last_frame - current_frame +1;
		set(t,'TasksToExecute',num_frames);
		VARS_GLOBAL_ANIM.current_frame = current_frame;
		
        if(new_time_stretch_factor ~= -1)
            desired_frame_time = VARS_GLOBAL_ANIM.mot(1).frameTime;
            period = round(1000*desired_frame_time/new_time_stretch_factor)/1000;
            VARS_GLOBAL_ANIM.timer_period = period;
            if (period == 0)
                warning(['Requested timer period of ' num2str(desired_frame_time/time_stretch_factor) ' s is too short for Matlab! Setting period to minimum of 1 ms :-(']);
                period = 0.001;
            end
		    set(t,'Period',period);      
        end
        
		start(t);
	end
catch
	delete(t);
    return
end

VARS_GLOBAL_ANIM.animation_paused = false;
