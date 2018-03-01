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

function animateD(doc_id,varargin)

global VARS_GLOBAL;

repeat_num = 1;
speed = 1;
frames = 0;
if nargin > 3
    repeat_num = varargin{3};
end
if nargin > 2
    speed = varargin{2};
end
if nargin > 1
    frames = varargin{1};
end
    
[skel,mot] = readMocapD(doc_id);
if (frames == 0)
    frames = 1:mot.nframes;
end
%animate_initGraphics;
figure(100);
set(gcf, 'name', [num2str(doc_id), ': ' mot.filename]);
animate(skel,mot,repeat_num,speed,frames);