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

function mot = addAnimatedPatches(mot,function_name,function_params,varargin)
% mot,function_name,function_params,type,color,alpha,overwrite)

type = cell(length(function_name),1);
color = cell(length(function_name),1);
alpha = zeros(length(function_name),1);
overwrite = false;
for k=1:length(function_name)
    type{k} = 'disc';
    color{k} = 'blue';
    alpha(k) = 1;
end
if (nargin>6)
    overwrite = varargin{4};
end
if (nargin>5)
    alpha = varargin{3};
end
if (nargin>4)
    color = varargin{2};
end
if (nargin>3)
    type = varargin{1};
end

if (overwrite | ~isfield(mot,'animated_patch_data'))
    mot.animated_patch_data = animatedPatchStruct(length(function_name));
    k0 = 0;
else
    k0 = length(mot.animated_patch_data);
end

for k = k0+1:k0+length(function_name)
    mot.animated_patch_data(k).function_name = function_name{k-k0};
    mot.animated_patch_data(k).function_params = function_params{k-k0};
    mot.animated_patch_data(k).type = type{k-k0};
    mot.animated_patch_data(k).color = color{k-k0};
    mot.animated_patch_data(k).alpha = alpha(k-k0);
end
