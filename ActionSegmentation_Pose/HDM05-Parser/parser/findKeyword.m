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

function [b, line, pos, line_count] = findKeyword(varargin)
% stops at first occurence of any of the keywords. NOT case sensitive!
% args: fid,keywords
if nargin < 2
    error('Not enough arguments!');
end
fid = varargin{1};

pos = 0;
line_count = 0;
line = [];
while ~feof(fid)
    l = eatWhitespace(fgetl(fid));
    line_count = line_count + 1;
    for i = 2:nargin
        k = strfind(upper(l),upper(varargin{i}));
        if size(k) > 0
            line = l(k:size(l,2));
            b = true;
            pos = ftell(fid);
            return;
        end
    end
end
b = false;
