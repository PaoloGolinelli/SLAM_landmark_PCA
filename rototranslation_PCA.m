function [R1to2, G1, G2] = rototranslation_PCA(P1, P2)

    % remove the baricenter from each set of points
    G1 = mean(P1);
    G2 = mean(P2);
    P1c = P1 - G1;
    P2c = P2 - G2;

    % find the rotation matrix through SVD
    [~, ~, V1] = svd(P1c);
    [~, ~, V2] = svd(P2c);

    R1to2 = V1 * V2';
end