function [F, Fs1, Fs2] = setInterpFunc(Element, pathDir, ReversePath)
    %Natural interpolation method is used to form a stress function to then
    %compute the paths.
    surr_elements = RunLibrary_surrounding_elemnts(Element, Element);
    nodes = unique([surr_elements(:).nodes]);
    coordx = [nodes(:).xCoordinate]';
    coordy = [nodes(:).yCoordinate]';
    coordz = [nodes(:).zCoordinate]';

    stress
    coordinates
    stressTensor = [stress(:,1), stress(:,2), stress(:,3)];

    switch pathDir
        case 'X'
            stress_tensor = [[nodes(:).xStress]', [nodes(:).xyStress]', [nodes(:).xzStress]'];
        case 'Y'
            stress_tensor = [[nodes(:).yStress]', [nodes(:).xyStress]', [nodes(:).yzStress]'];
        case 'Z'
            stress_tensor = [[nodes(:).zStress]', [nodes(:).xzStress]', [nodes(:).yzStress]'];
    end
    F = scatteredInterpolant(coordx, coordy, coordz, stress_tensor(:,1), 'natural');
    Fs1 = scatteredInterpolant(coordx, coordy, coordz, stress_tensor(:,2), 'natural');
    Fs2 = scatteredInterpolant(coordx, coordy, coordz, stress_tensor(:,3), 'natural');
    %Flips the stress function when the backwards path is being computed.
    if ReversePath
        F=@(x,y, z) -F(x,y, z);
        Fs1=@(x,y,z) -Fs1(x,y,z);
        Fs2=@(x,y,z) -Fs2(x,y,z);
    end
end