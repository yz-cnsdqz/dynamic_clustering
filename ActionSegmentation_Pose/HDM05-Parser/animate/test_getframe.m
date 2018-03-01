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

figure('units','normalized','outerposition',[0 0 1 1]); 


% outerposition is an undocumented feature. It allows you to 
% define the size of the outerposition of the figure window. 
% The above sets the figure to full screen or you could also 
% try uncommenting the statements after the surf(peaks) 


plot([0:99],sin(1/100*[0:99]*2*pi));  
%disp ('Maximize the screen, then press any key to continue'); 
%pause; 


mov=avifile('test','compression','indeo5'); 
%makes test.avi w/ Indeo 5 


%f2=getframe(gcf); % gets the gcf  
mov=addframe(mov,gcf); % adds the frame into mov  
mov=close(mov); % closes the mov  


%In order to run the avi from MATLAB command window:  
%!test.avi&  