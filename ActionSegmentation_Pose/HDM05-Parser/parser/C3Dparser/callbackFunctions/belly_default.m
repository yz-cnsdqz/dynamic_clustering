function trajectories = belly_default( mot )

viconNames = {'STRN', 'T10', 'RFWT', 'LFWT', 'LBWT', 'RBWT'};

for i = 1:length(viconNames)
    idx(i) = strMatch(viconNames(i), mot.nameMap(:,1));
end
    
trajectories = 0.4 * mot.jointTrajectories{idx(1)} + 0.4 * mot.jointTrajectories{idx(2)} + 0.05 * mot.jointTrajectories{idx(3)}  ...
             + 0.05 * mot.jointTrajectories{idx(4)} + 0.05 * mot.jointTrajectories{idx(5)} + 0.05*mot.jointTrajectories{idx(6)};
