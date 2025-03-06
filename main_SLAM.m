close all
clear 
clc

%% Note
% In this file the matching of consecutive lidar scans is performed in
% order to estimate the pose of the robot

% if it starts from different scan (RefId != 0) it might not work
% if the threshold is changed (!= 0.1) it might not work 
% there are too few landmarks to work in reliable and consistent way
% also the SVD might find the opposite direction to the correct one when
% performed on multiple landmarks
% SLAM performed until 570 scan, because after finds only one landmark

%% Initialization and parameters

load lidarScans.mat

% Scan analysis data
Res = 4;                % span of point window for landmark extractions algorithm
randsac_trhs = 0.1;     % threshold for randsac algorithm
Nrand = 2000;           % repetitions of randsac algorithm
RefId = 1;              % starting scan index
Dscan = 15;             % how many scans it skips while performing SLAM
EndId = 570;            % ending scan index

pplot = 0;
nplot = 3;

% initialize pose of the robot
r_yaw = 0;      % heading angle of the robot
r_pos = [0,0];  % xy coordinates of the robot
rob_pose = [r_pos(1), r_pos(2), r_yaw]; % list containing all robot poses

% extract the reference (first) scan
ReferenceScan = lidarScans(RefId);
ReferencePoints = ReferenceScan.Cartesian;
cell_array_reference = keyP2(ReferenceScan, Res, pplot, nplot); % extract landmarks
ScanCart_reference_Landmarks = cell_array_reference{2}; % extract landmarks coords

% initialize the plot
figure(30)
plot(ReferencePoints(:,1),ReferencePoints(:,2),'k.');
hold on

% initialize the occupancy map
mapR = occupancyMap(20, 20, 10);
mapR.GridLocationInWorld = [-10, -10]; % [x,y] world cooridnates of the bottom-left corner
insertRay(mapR, [r_pos(1),r_pos(2),r_yaw], ReferenceScan, 10, [0.1 0.9]); % add first scan

%% SLAM
% cycle over all the scans from RefId to EndId by steps of Dscan
for ScanId = RefId+Dscan:Dscan:EndId

    % RANDSAC implementation:
    % EXTRACT LANDMARKS from current scans
    CurrentScan = lidarScans(ScanId);
    CurrentPoints = CurrentScan.Cartesian;
    cell_array_current = keyP2(CurrentScan, Res, pplot, nplot);
    ScanCart_current_Landmarks = cell_array_current{2};

    NConsMax = 0; % initialize the max number of landmarks of consensus
    
    % perform Nrand iterations of RANSAC algorithm
    for i = 1:Nrand

        % randomly sample two landmarks from each scan
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

        % rotate all landmarks  of current scan
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
                ScanCurrentLandRot(IdxMin, :) = [NaN, NaN];  % remove the matched landmark
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

    fprintf("scan %d, matched %d landmarks\n", ScanId, NConsMax)

    % refine the rototranslation by calculating it over all the consensus
    % set (not always working) 
    %[Rmat, G1, G2] = rototranslation_PCA(ConsensusSet1Max, ConsensusSet2Max);

    % try to fix the problem described in the first section
    % angleR = atan2(Rmat(1,2),Rmat(1,1));
    % angleRbest = atan2(Rbest(1,2),Rbest(1,1));
    % if (abs(angleRbest-angleR) > 3)
    %     fprintf("need a rotation %d\n", angleR-angleRbest)
    %     Rmat = [-1 0; 0 -1]*Rmat;
    % elseif (abs(angleRbest-angleR) > 0.5) 
    %     fprintf("don't use refined transformation %d\n", angleR-angleRbest)
    %     Rmat = Rbest; G1 = G1best; G2 = G2best;
    % end
    
    % better to use the rototransformation calculated over the two landmarks only
    Rmat = Rbest; G1 = G1best; G2 = G2best;

    % rotate the current scan to match the reference scan
    CurrentPoints_rotated = (CurrentScan.Cartesian-G2)*Rmat' + G1;
    ScanCurrentLandRot = (ScanCart_current_Landmarks-G2)*Rmat' + G1;
    ConsensusSet2MaxRot = (ConsensusSet2Max-G2)*Rmat' + G1;

    % plot the two matched scans and the landmark in the consensus
    figure(10)
    clf(10)
    plot(ReferencePoints(:,1), ReferencePoints(:,2),'r.', 'MarkerSize',10);
    hold on
    plot(CurrentPoints_rotated(:,1), CurrentPoints_rotated(:,2),'g.', 'MarkerSize',10);
    plot(ConsensusSet1Max(:,1), ConsensusSet1Max(:,2),'ro', 'MarkerSize',10);
    plot(ConsensusSet2MaxRot(:,1), ConsensusSet2MaxRot(:,2),'go', 'MarkerSize',10);
    legend(sprintf('Scan %g', ScanId))

    % update robot's pose according to the rototranslation found
    r_yaw = atan2(Rmat(2,1),Rmat(1,1)); % tan(yaw) = sin(yaw)/cos(yaw)
    if (r_yaw > pi/2) r_yaw = r_yaw - 2*pi; end % mirror (plotting purposes)
    
    r_pos = ([0,0]-G2)*Rmat' + G1;  % robot is at the center of the scan ([0,0])
    rob_pose = [rob_pose; [r_pos(1),r_pos(2), r_yaw]]; % update dynamically list of poses

    % add new rays in the occupancy map
    insertRay(mapR, [r_pos(1),r_pos(2),r_yaw], CurrentScan, 10, [0.1 0.9]);

    % plot the occupancy map
    figure(20)
    show(mapR);
    title('ReferenceScan Occupancy grid map');

    % add the new rotated point to the plot
    figure(30)
    plot(CurrentPoints_rotated(:,1), CurrentPoints_rotated(:,2),'k.');
    if (ScanId > RefId + Dscan)
        line([rob_pose(end-1,1),r_pos(1)], [rob_pose(end-1,2),r_pos(2)], 'Color','b', 'LineWidth',2)
    end
    legend(sprintf('Scan %g', ScanId))

    pause(0.1)  % allow to update plots
    % pause(1)    % slow computation to debug

    % current scan rotated becomes reference of next iteration
    ReferencePoints = CurrentPoints_rotated;
    ScanCart_reference_Landmarks = ScanCurrentLandRot;
end

% plot the trajectory of the robot and the angle over time
figure(40)
subplot(2,1,1);
plot(rob_pose(:,1),rob_pose(:,2), 'LineWidth',2)
title('Robot position')
xlabel('x [m]')
ylabel('y [m]')
subplot(2,1,2); 
plot(rob_pose(:,3), 'LineWidth',2)
title('Robot heading')
xlabel('iteration')
ylabel('angle [°]')