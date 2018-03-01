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

function [X,Y,Z] = createPlaneNormal(n, x0, sides_length)

X = [-sides_length(1)/2 -sides_length(1)/2;sides_length(1)/2 sides_length(1)/2];
Y = [sides_length(2)/2 -sides_length(2)/2;sides_length(2)/2 -sides_length(2)/2];
Z = [0 0;0 0];
n0 = [0;0;1];
x = cross(n,n0);
cphi = dot(n/norm(n),n0/norm(n0));


rotate(h,x,phi*180/pi);

set(h,'facecolor','black');
alpha(h,0.25);
