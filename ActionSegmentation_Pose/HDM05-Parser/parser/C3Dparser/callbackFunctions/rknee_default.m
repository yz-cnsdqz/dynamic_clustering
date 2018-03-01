function trajectories = rknee_default( mot )

viconNames = {'RKNE', 'RTOE', 'RMT5'};

for i = 1:length(viconNames)
    idx(i) = strMatch(viconNames(i), mot.nameMap(:,1));
end

% Approximation: knee joint is located 2/5 foot-width from surface
dist = mot.jointTrajectories{idx(2)} - mot.jointTrajectories{idx(3)};

trajectories = mot.jointTrajectories{idx(1)} + 0.4 * dist;

