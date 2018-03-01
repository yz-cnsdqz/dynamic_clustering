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

function [skel,mot] = readMocapD(file_num, range)

global VARS_GLOBAL

if ~(exist('DB_concat'))
    DB_concat = DB_index_load('HDM05_amc',{'AK_lower'},4);
    files_name = DB_concat.files_name;
end



amcfullpath = fullfile(VARS_GLOBAL.dir_root, files_name{file_num});
[info,OK] = filename2info(amcfullpath);
if (OK)
    [skel,mot] = readMocap([info.amcpath info.asfname], amcfullpath, [1 inf]);
    if (nargin > 2)
        mot = cropMot(mot, range);
    end
elseif strcmpi(amcfullpath(end-4:end),'.bvh')
    [skel,mot] = readMocap(amcfullpath);
else % assume an AMC file with an ASF that has the same name
    [skel,mot] = readMocap([amcfullpath(1:end-4) '.asf'],amcfullpath);
end
