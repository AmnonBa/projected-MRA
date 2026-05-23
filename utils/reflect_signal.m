function thetaRef = reflect_signal(theta)
    p = numel(theta);
    thetaRef = zeros(p,1);

    thetaRef(1) = theta(1);

    for j = 1:p-1
        reflectedIndex = mod(-j, p);
        thetaRef(j+1) = theta(reflectedIndex + 1);
    end
end
