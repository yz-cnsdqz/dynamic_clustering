function mot = convert2quat(skel,mot)

switch lower(mot.angleUnit)
    case 'deg'
        conversion_factor = pi/180;
    case 'rad'
        conversion_factor = 1;
    otherwise
        error(['Unknown angle unit: ' mot.angleUnit]);
end
nTraj = size(mot.rotationEuler,1);
mot.rotationQuat = cell(nTraj, 1);
for k = 1:nTraj
    node = skel.nodes(k);
    nRotDOF = size(mot.rotationEuler{node.ID,1},1); % number of rotational DOFs
    
    % zero-pad Euler array at the appropriate places, according to rotation order and presence/absence of respective DOFs
    if nRotDOF > 0
        completeEulers = zeros(3,mot.nframes);
        d = 1; % index for DOFs present in mot.rotationEuler
        for r = 1:3 % go through rotationOrder. 
        	idx = strmatch(['r' lower(node.rotationOrder(r))], lower(node.DOF), 'exact');
            if (~isempty(idx))
                completeEulers(r,:) = mot.rotationEuler{node.ID,1}(d,:)*conversion_factor;
                d = d+1;
            end
        end
        axis_quat = repmat(euler2quat(flipud(node.axis)*conversion_factor,'zyx'),1,mot.nframes); % According to ASF specs, rotation order for "axis" should be XYZ. However, they use the opposite multiplication order as we do!
        if (node.ID == 1) % root node? => special case for determination of rotation order (node.rotationOrder only concerns the global rotational offset in this case)
            rootTransformationOrder = char(node.DOF);
            rootTransformationOrder = rootTransformationOrder(:,2)';
            [x,IA,IB] = intersect({'rx','ry','rz'}, lower(node.DOF));
            rootRotationOrder = rootTransformationOrder(sort(IB));
            mot.rotationQuat{node.ID,1} = euler2quat(flipud(completeEulers),fliplr(rootRotationOrder)); % ASF specs use opposite multiplication order as we do, hence fliplr() and flipud()!
        else            
            mot.rotationQuat{node.ID,1} = quatmult(axis_quat,quatmult(euler2quat(flipud(completeEulers),fliplr(node.rotationOrder)),quatinv(axis_quat))); % ASF specs use opposite multiplication order as we do, hence fliplr() and flipud()!
        end
    end
end
