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

function varargout = animateGUI(varargin)
% animateGUI M-file for animateGUI.fig
%      animateGUI, by itself, creates a new animateGUI or raises the existing
%      singleton*.
%
%      H = animateGUI returns the handle to a new animateGUI or the handle to
%      the existing singleton*.
%
%      animateGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in animateGUI.M with the given input arguments.
%
%      animateGUI('Property','Value',...) creates a new animateGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before animateGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to animateGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help animateGUI

% Last Modified by GUIDE v2.5 18-Jul-2006 17:52:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @animateGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @animateGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes just before animateGUI is made visible.
function animateGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to animateGUI (see VARARGIN)

% Choose default command line output for animateGUI
handles.output = hObject;

handles.skel = varargin{1};
handles.mot = varargin{2};

mot_num = length(handles.mot);

frames_total = max([handles.mot.nframes]);

set(handles.slider_animate, 'Min',0.9,'Max', frames_total +.1,'Value', 1,'SliderStep', [1/frames_total 1/frames_total]);

% reset animation controls
set(handles.slider_animate, 'Enable','off');
set(handles.slider_speed,'Enable','off');
set(handles.pushbutton_pause,'Enable','off');
set(handles.pushbutton_play,'Enable','off'); 

% Update handles structure
guidata(hObject, handles);

cameratoolbar;
cameratoolbar('SetCoordSys','y');

handles = startAnimation(hObject,handles);

% frame counter
global VARS_GLOBAL_ANIM;
VARS_GLOBAL_ANIM.graphics_data.frameLabel = handles.text_frameCounter;

% UIWAIT makes animateGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Outputs from this function are returned to the command line.
function varargout = animateGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Start animating
function handles = startAnimation(hObject,handles)
% hObject   handle to UI-control calling this function
% handles   structure with handles and user data (see GUIDATA)
global VARS_GLOBAL_ANIM
% pause the animation that might currently be running
pauseAnimation;

% prepare animation controls
frames_total = max([handles.mot.nframes]);
set(handles.slider_animate, 'Enable','on','Max', frames_total+.001,'Value', 1,'SliderStep', [1/frames_total 1/frames_total]);
set(handles.slider_speed, 'Enable','on');
set(handles.pushbutton_play, 'Enable','on');
set(handles.pushbutton_pause, 'Enable','on');
speed = get(handles.slider_speed,'Value');


% set some properties for the animation
VARS_GLOBAL_ANIM = emptyVarsGlobalAnimStruct;
VARS_GLOBAL_ANIM.figure_color = get(gcf,'Color');
VARS_GLOBAL_ANIM.animated_skeleton_Color = [0 0 0]/255;
VARS_GLOBAL_ANIM.animated_skeleton_Marker = '.';
VARS_GLOBAL_ANIM.animated_skeleton_MarkerEdgeColor = [1 0 0];
VARS_GLOBAL_ANIM.animated_skeleton_MarkerFaceColor = [1 0 0];
VARS_GLOBAL_ANIM.kill_timer = false;
VARS_GLOBAL_ANIM.animated_skeleton_MarkerSize = 8;
VARS_GLOBAL_ANIM.animated_point_MarkerSize = 12;
VARS_GLOBAL_ANIM.animated_skeleton_LineWidth = 4;
VARS_GLOBAL_ANIM.ground_tile_size_factor = 1;
VARS_GLOBAL_ANIM.bounding_box_border_extension = 0.01;

%VARS_GLOBAL_ANIM.figure_position = [5 5 512 384];
%rehash path

set(gcf,'CurrentAxes',handles.axes_animate);
cla reset;
% start animation
VARS_GLOBAL_ANIM.animation_paused = false;
new_animate(handles.skel,handles.mot,1,speed,{},{},{'updateAnimationSlider',handles.slider_animate},true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton_pause.
function pushbutton_pause_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pauseAnimation;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton_play.
function pushbutton_play_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global VARS_GLOBAL_ANIM;

speed = get(handles.slider_speed,'Value');
frame = round(get(handles.slider_animate,'Value'));
max_frame = round(get(handles.slider_animate,'Max'));

if (frame == max_frame)
    resumeAnimation(1,speed);
elseif (VARS_GLOBAL_ANIM.animation_paused) 
    resumeAnimation(frame,speed);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function slider_animate_Callback(hObject, eventdata, handles)
% hObject    handle to slider_animate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global VARS_GLOBAL_ANIM

was_already_paused = VARS_GLOBAL_ANIM.animation_paused;
pauseAnimation;
current_frame = round(get(hObject,'Value'));
VARS_GLOBAL_ANIM.current_frame = current_frame;

new_animate_showFrame([],[],{},current_frame);

if (~was_already_paused)
    resumeAnimation;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function slider_animate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_animate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pauseAnimation;
t = timerfind('Name','AnimationTimer');
if (~isempty(t))
    delete(t);
end
delete(gcf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function slider_speed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function slider_speed_Callback(hObject, eventdata, handles)
% hObject    handle to slider_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global VARS_GLOBAL_ANIM
% adjust animation speed
if(~VARS_GLOBAL_ANIM.animation_paused)
    pauseAnimation;
    current_frame = VARS_GLOBAL_ANIM.current_frame;
    speed = get(hObject,'Value');
    resumeAnimation(current_frame,speed);
else
    speed = get(hObject,'Value');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on state change of checkbox_loop.
function checkbox_loop_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_loop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_loop
global VARS_GLOBAL_ANIM
VARS_GLOBAL_ANIM.loop_playback = get(hObject,'Value');


% --- Executes during object creation, after setting all properties.
function edit_gotoFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gotoFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit_gotoFrame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gotoFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gotoFrame as text
%        str2double(get(hObject,'String')) returns contents of edit_gotoFrame as a double

global VARS_GLOBAL_ANIM

was_already_paused = VARS_GLOBAL_ANIM.animation_paused;
pauseAnimation;
current_frame = round(str2num(char(get(hObject,'String'))));
VARS_GLOBAL_ANIM.current_frame = current_frame;

new_animate_showFrame([],[],{},current_frame);

if (~was_already_paused)
    resumeAnimation;
end
