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

function handle = showTrajectory(varargin)
% showTrajectory(mot,trajname,linestyle,downsampling_fac)

switch (nargin)
    case 2
        mot = varargin{1};
        trajname = varargin{2};
        linestyle = 'red o';
        downsampling_fac = 1;
    case 3
        mot = varargin{1};
        trajname = varargin{2};
        linestyle = varargin{3};
        downsampling_fac = 1;
    case 4
        mot = varargin{1};
        trajname = varargin{2};
        linestyle = varargin{3};
        downsampling_fac = varargin{4};
    otherwise
        error('Wrong number of arguments!');
end
if (ischar(trajname))
    ID = trajectoryID(mot,trajname);
elseif (isnum(trajname))
    ID = trajname;
else
    error('Expected trajectory name or numeric trajectory ID!');
end

if (~ishold)
    hold;
end

handle = plot3(mot.jointTrajectories{ID}(1,1:downsampling_fac:end),mot.jointTrajectories{ID}(2,1:downsampling_fac:end),mot.jointTrajectories{ID}(3,1:downsampling_fac:end),linestyle);