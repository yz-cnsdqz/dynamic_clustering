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

cs = emptySkeleton;
cm = emptyMotion;
doc_IDs = [20:33];
for k=1:length(doc_IDs);
    [cs(k),cm(k)] = readMocapD(doc_IDs(k));
    for j=1:length(cm(k).jointTrajectories)
%        offset_xz = (k-1)*[5;0;5]; % skeletons placed along diagonal
        chkr_width = 4;
        x_ofs = 10; z_ofs = 10;
        offset_xz = [mod((k-1),chkr_width)*x_ofs;0;fix((k-1)/chkr_width)*z_ofs]; %  skeletons placed in checkerboard formation
        cm(k).jointTrajectories{j} = cm(k).jointTrajectories{j} + repmat(offset_xz,1,size(cm(k).jointTrajectories{j},2));
        cm(k).boundingBox = computeBoundingBox(cm(k));
    end
end