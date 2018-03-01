function [y_left_arm, y_right_arm, y_torso] = ...
    readTUMKitchenAnnotation(anno)    


%% read annotation from the file
left_arm_states = {};
right_arm_states = {};
torso_states = {};
frame_idx = [];
for ii = 2:length(anno)
    gt_info = strsplit(anno{ii},',');
    frame_idx = [frame_idx; str2num(gt_info{1})];

% 
%         if ~isempty(left_arm_states)
%             if ~strcmp(gt_info{2}, left_arm_states(end))
%                 timestamp_left_arm_change = [timestamp_left_arm_change; str2num(gt_info{1})];
%             end
%         end
% 
%         if ~isempty(right_arm_states)
%             if ~strcmp(gt_info{3}, right_arm_states(end))
%                 timestamp_right_arm_change = [timestamp_right_arm_change; str2num(gt_info{1})];
%             end
%         end
% 
% 
%         if ~isempty(torso_states)
%             if ~strcmp(gt_info{4}, torso_states(end))
%                 timestamp_torso_change = [timestamp_torso_change; str2num(gt_info{1})];
%             end
%         end

    left_arm_states{end+1} = (gt_info{2});
    right_arm_states{end+1} =(gt_info{3});
    torso_states{end+1} = (gt_info{4});
end

%% encode string labels to numbers
torso_labels = {'StandingStill','HumanWalkingProcess'};
hand_labels = {'Reaching','TakingSomething','LoweringAnObject','ReleasingGraspOfSomething',...
    'OpeningADoor','ClosingADoor','OpeningADrawer','ClosingADrawer','CarryingWhileLocomoting'};

y_left_arm = zeros(length(anno)-1,1);
y_right_arm = zeros(length(anno)-1,1);
y_torso = zeros(length(anno)-1,1);

for ii = 1:length(torso_labels)
    idx = find(contains(torso_states,torso_labels{ii}));
    y_torso(idx) = ii;
end

for ii = 1:length(hand_labels)
    idx = find(contains(left_arm_states,hand_labels{ii}));
    y_left_arm(idx) = ii;
    
    idx = find(contains(right_arm_states,hand_labels{ii}));
    y_right_arm(idx) = ii;
    
end
    





end