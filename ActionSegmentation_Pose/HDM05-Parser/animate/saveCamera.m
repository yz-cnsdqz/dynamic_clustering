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

function saveCamera(varargin)

ax = gca;
file = 'camera.m';
if (nargin>1)
    file = varargin{2};
end
if (nargin>0)
    ax = varargin{1};
end

[fid,msg]=fopen(file,'w');
if (fid < 0)
    disp(msg);
    return;
end

fprintf(fid,'function camera(axs)\n\n');
fprintf(fid,'set(axs,''CameraPosition'', [%s],...\n',num2str(get(ax,'CameraPosition')));
fprintf(fid,'        ''CameraPositionMode'',''manual'',...\n');
fprintf(fid,'        ''CameraTarget'', [%s],...\n',num2str(get(ax,'CameraTarget')));
fprintf(fid,'        ''CameraTargetMode'',''manual'',...\n');
fprintf(fid,'        ''CameraUpVector'', [%s],...\n',num2str(get(ax,'CameraUpVector')));
fprintf(fid,'        ''CameraUpVectorMode'',''manual'',...\n');
fprintf(fid,'        ''CameraViewAngle'', [%s],...\n',num2str(get(ax,'CameraViewAngle')));
fprintf(fid,'        ''CameraViewAngleMode'',''manual'');');
fclose(fid);