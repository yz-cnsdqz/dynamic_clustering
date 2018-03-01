function [ skel_new, mot_new ] = generateSkel( skel, mot, callbackFunctionNames)

if nargin < 3
    callbackFunctionNames = generateDefaultFunctionNames;
end

skel_new = skel;
skel_new.nameMap = [];
skel_new.njoints = 24;
mot_new = mot;
mot_new.nameMap = [];
mot_new.jointTrajectories = [];
mot_new.njoints = 24;

for i = 1:size(callbackFunctionNames,2)
    traj = feval(char(callbackFunctionNames(2,i)), mot);
    mot_new.jointTrajectories{i,1} = traj;
    mot_new.nameMap{i,1} = char(callbackFunctionNames(1,i));
    mot_new.nameMap{i,2} = 0;
    mot_new.nameMap{i,3} = i;
end
    
skel_new.nameMap = mot_new.nameMap;
skel_new.paths = buildSkelPathsFromMot(mot_new);

