close all
clear 
clc

%% Note
% In this file the matching of two scan only is performed
% Used for debugging purposes

%% Main RANSAC: extract the data

load lidarScans.mat

ReferenceScan = lidarScans(201);
CurrentScan = lidarScans(211);

% Scan analysis data
Res = 4;            % span of point window for landmark extractions algorithm
Nrand = 1000;       % repetitions of randsac
randsac_trhs = 0.1; % randsac threshold

pplot = 0;
nplot = 3;

%% RANDSAC implementation

% EXTRACT LANDMARKS from reference and current scans
cell_array_reference = keyP2(ReferenceScan, Res, pplot, nplot);
cell_array_current = keyP2(CurrentScan, Res, pplot, nplot);
ReferencePoints = ReferenceScan.Cartesian;

% matrix containing the landmark coordinates:
ScanCart_reference_Landmarks = cell_array_reference{2};
ScanCart_current_Landmarks = cell_array_current{2};

NConsMax = 0;

% perform Nrand iterations of RANSAC algorithm
for i = 1:Nrand
   
    terminating_cond = true;
    while( terminating_cond )
        % extract 2 random landmarks from each scan
        landmark_sample_reference = ceil(rand(2,1) * length(ScanCart_reference_Landmarks));
        landmark_sample_current = ceil(rand(2,1) * length(ScanCart_current_Landmarks));

        p11 = ScanCart_reference_Landmarks(landmark_sample_reference(1),:);
        p12 = ScanCart_reference_Landmarks(landmark_sample_reference(2),:);

        p21 = ScanCart_current_Landmarks(landmark_sample_current(1),:);
        p22 = ScanCart_current_Landmarks(landmark_sample_current(2),:);

        % ensure that the indeces are different between each others
        terminating_cond = (landmark_sample_reference(1) == landmark_sample_reference(2));
        terminating_cond = or(terminating_cond, (landmark_sample_current(1) == landmark_sample_current(2)));
        % ensure the distance between the two landmarks is closer than threshold
        terminating_cond = or(terminating_cond, abs(norm(p12-p11) - norm(p22-p21)) >= randsac_trhs);
    end

    % find the corresponding rototranslation between landmark_sample_reference and landmark_sample_current
    ps1 = zeros(2);
    ps1(1,:) = p11;
    ps1(2,:) = p12;

    ps2 = zeros(2);
    ps2(1,:) = p21;
    ps2(2,:) = p22;

    [Rmat, G1, G2] = rototranslation_PCA(ps1, ps2);

    % rototranslate ALL the landmarks of currentScan on the referenceScan
    % CurrentPoints_rotated = (CurrentScan.Cartesian-G2)*Rmat' + G1;
    % ps2Rot = (ps2 - G2)*Rmat' + G1;
    ScanCurrentLandRot = (ScanCart_current_Landmarks-G2) * Rmat' + G1;

    % find the consensus set, set a threshold of 0.2 meters 
    % and save the corresponding cardinality on a variable in order to find
    % at the end the solution with the higher cardinality

    Nconsensus = 0;
    ConsensusSet1 = [];
    ConsensusSet2 = [];

    % iterate over all landmarks of reference scan
    L1s = length(ScanCart_reference_Landmarks);
    for j = 1:L1s
        Pj = ScanCart_reference_Landmarks(j,:);

        % calculate distance to each landmark of current scan
        dist_to_land2 = vecnorm(ScanCurrentLandRot - Pj, 2, 2);
        % extract the minimum distance 
        [Dmin,IdxMin] = min(dist_to_land2);
        % if min distance is smaller than threshold then belong to consensus
        if (Dmin < randsac_trhs)
            Nconsensus = Nconsensus + 1;
            ConsensusSet1 = [ConsensusSet1; Pj];
            ConsensusSet2 = [ConsensusSet2; ScanCart_current_Landmarks(IdxMin,:)];
            ScanCurrentLandRot(IdxMin, :) = [NaN,NaN];  % remove the matched landmark
        end
    end

    % save the best, yet found, consensus
    if (Nconsensus > NConsMax)
        NConsMax = Nconsensus;
        ConsensusSet1Max = ConsensusSet1;
        ConsensusSet2Max = ConsensusSet2;
        Rbest = Rmat; G1best = G1; G2best = G2;
    end
end 

% find the corresponding rototranslation between
% landmark_consensusSet_reference and landmark_consensusSet_current
% with SVD
[Rmat, G1, G2] = rototranslation_PCA(ConsensusSet1Max, ConsensusSet2Max);

angleR = atan2(Rmat(1,2),Rmat(1,1));
angleRbest = atan2(Rbest(1,2),Rbest(1,1));
if (abs(angleRbest-angleR) > 3)
    fprintf("need a rotation %d\n", angleR-angleRbest)
    Rmat = [-1 0; 0 -1]*Rmat;
elseif (abs(angleRbest-angleR) > 0.5) 
    fprintf("don't use refined transformation %d\n", angleR-angleRbest)
    Rmat = Rbest; G1 = G1best; G2 = G2best;
end
%Rmat = Rbest; G1 = G1best; G2 = G2best; 
sprintf('Matching based on %d landmarks\n', NConsMax)

% rotate the current scan to match the reference scan
CurrentPoints_rotated = (CurrentScan.Cartesian-G2)*Rmat' + G1;
ScanCurrentLandRot = (ScanCart_current_Landmarks-G2)*Rmat' + G1;
ConsensusSet2MaxRot = (ConsensusSet2Max-G2)*Rmat' + G1;

% plot the two matched scans and the landmark in the consensus
figure(30);
plot(ReferencePoints(:,1),ReferencePoints(:,2),'r.');
hold on
plot(CurrentPoints_rotated(:,1),CurrentPoints_rotated(:,2),'g.');
plot(ConsensusSet1Max(:,1), ConsensusSet1Max(:,2),'ro', 'MarkerSize',10);
plot(ConsensusSet2MaxRot(:,1), ConsensusSet2MaxRot(:,2),'go', 'MarkerSize',10);

