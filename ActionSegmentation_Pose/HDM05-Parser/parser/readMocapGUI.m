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

function [skel, mot] = readMocapGUI(noCaching);
% [skel, mot] = readMocapGUI(noCaching)

if nargin < 1
    noCaching = false;
end

global VARS_GLOBAL;

oldPath = cd;

try
    defaultPath = VARS_GLOBAL.readMocapGUILastDir;
catch
    try
        defaultPath = VARS_GLOBAL.dir_root;
    catch
        defaultPath = cd;
    end
end

cd(defaultPath);

[datafile,datapath] = uigetfile('*.c3d; *.amc', 'Choose data file', 40, 40);
if datafile ~= 0    
    VARS_GLOBAL.readMocapGUILastDir = datapath;
    
    idx = findstr(datafile,'.');
    ext = upper(datafile(idx(end):end));
    
    if strcmp(ext, '.C3D')
        [skel, mot] = readMocap([datapath,datafile], [], noCaching);
    elseif strcmp(ext, '.AMC')
        asfFiles = dir([datapath '\*.ASF']);
        if isempty(asfFiles)
            cd(datapath);
            [asfFile, asfPath] = uigetfile('*.asf', 'Choose ASF file', 40, 40);
        else
            asfPath = datapath;
            asfFile = [datafile(1:6) '.ASF'];
        end
        [skel, mot] = readMocap([datapath,asfFile], [datapath,datafile], [], true, true, true, noCaching);
    end
end

cd(oldPath);