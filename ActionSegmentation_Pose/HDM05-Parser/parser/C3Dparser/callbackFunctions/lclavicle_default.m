function trajectories = lclavicle_default( mot )

viconNames = {'C7', 'LSHO'};

for i = 1:length(viconNames)
    idx(i) = strMatch(viconNames(i), mot.nameMap(:,1));
end
    
% can be adjusted to slide over to shoulder if needed: alpha \in [0:1]
alpha = 1.0;
trajectories = alpha * mot.jointTrajectories{idx(1)} + (1-alpha) * mot.jointTrajectories{idx(2)};
