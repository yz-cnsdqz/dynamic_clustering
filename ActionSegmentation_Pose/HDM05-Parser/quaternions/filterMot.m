function motOut = filterMot(skel,mot,w,method)
motOut = mot;
for k=1:length(mot.animated)
    switch (lower(method))
        case 'r4'
        motOut.rotationQuat{mot.animated(k)} = filterR4(w,mot.rotationQuat{mot.animated(k)},1,'sym');
        case 's3'
        b = orientationFilter(w);
        motOut.rotationQuat{mot.animated(k)} = filterS3(b,mot.rotationQuat{mot.animated(k)},'sym');
        case 'sphericalaverage'
        motOut.rotationQuat{mot.animated(k)} = filterSphericalAverageA1(w,mot.rotationQuat{mot.animated(k)},1,'sym');
    end        
end
motOut.jointTrajectories = forwardKinematicsQuat(skel,motOut);
motOut.boundingBox = computeBoundingBox(motOut);