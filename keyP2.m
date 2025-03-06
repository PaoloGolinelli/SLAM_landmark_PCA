function [keyPoints] = keyP2(Scan, Res, pplot, nplot)
    % Function that compute the curvature
    %%%% INPUT
    % Scan           Scan
    %   lidarScan with properties:
    %        Ranges: [N—1 double]
    %        Angles: [N—1 double]
    %     Cartesian: [N—2 double]
    %         Count: N
    %
    %   Res          +- range of points to compute indicators
    %                  
    %   pplot        == 1 plot the output
    %
    %   nplot        number of figure
    %
    %   ind          index of the scan point to plot
    %                  
    %%%% OUTPUT
    % keyPoints 1×3 cell array:
    % keyPoints{1} = indices n
    % keyPoints{2} = cartesian coordinates n x 2
    % keyPoints{3} = properties n x 4 [curvature, PCA angle, 1st comp, 2nd
    % comp]
    
    %%%% DEBUG SECTION:
    % close all
    % clear 
    % clc
    % load lidarScans.mat
    % Scan = lidarScans(90);
    % Res = 4;
    % pplot = 1;
    % nplot = 2;
    % figure(nplot-1)
    % plotScan([Scan], nplot-1)
    
    %%%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    %%%% Descriptors (Curvature, PCA, ...) profiles
    %%%%,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    ScanCart = Scan.Cartesian;
    ScanRange = Scan.Ranges;
    ll = length(ScanCart);
    % D_profile is a matrix containing the descriptors for each point:
    % - Column 1: Curvature at each point.
    % - Column 2: Angle between PCA components.
    % - Column 3: Eigenvalue of the first PCA component.
    % - Column 4: Eigenvalue of the second PCA component.
    D_profile = zeros(ll, 4);
    
    for i = 1+Res : ll-Res
        
        % ...... Curvature ......
        G = zeros(1,3);
        F = zeros(1,3);
        H = zeros(1,3);
        G(1:2) = Scan.Cartesian(i-Res,:);
        F(1:2) = Scan.Cartesian(i,:);
        H(1:2) = Scan.Cartesian(i+Res,:);
    
        A = norm(cross(H-F,G-F))/2;
    
        curvature = 4*A/(norm(G-H)*norm(F-H)*norm(G-F));
        
        D_profile(i, 1) = curvature;
    
        % ...... PCA ......
        % Given a data matrix X [n-by-p]: 
        % rows correspond to observations columns to variables
        % [coeff, score, latent] = pca(X)  returns:
        % - The coefficient matrix [p-by-p]. Each column of coeff contains the versor
        % for one principal component, the columns are in descending order 
        % of component variance. 
        % - the principal component scores in score (points reprojected in the
        % main reference system)
        % - the principal component variances in vector latent [p-by-1]
        % NOTE: By default, pca centers the data and uses the singular value decomposition (SVD) algorithm.
        % 
        X = Scan.Cartesian(i-Res:i+Res,:);
        [coeff, score, latent] = pca(X);
    
        % extract the angle of the first component
        princ_dir = coeff(:,1);
        theta = atan2(princ_dir(2),princ_dir(1));
    
        % extract the eigenvalues of the two principal components
        PCAcomp1 = latent(1);
        PCAcomp2 = latent(2);
    
        D_profile(i, 2) = theta;
        D_profile(i, 3) = PCAcomp1;
        D_profile(i, 4) = PCAcomp2;         
    end
    
    % Plot the results
    if (pplot)
        Screen = get(0, 'ScreenSize');
        figFeatures = figure(nplot+1);
        figFeatures.Position = [Screen(1) Screen(2) Screen(3)/3-20 Screen(4)];
        
        subplot(4,1,1), plot(D_profile(:,1), 'r');
        subplot(4,1,2), plot(D_profile(:,2), 'b');
        subplot(4,1,3), plot(D_profile(:,3), 'g');
        subplot(4,1,4), plot(D_profile(:,4), 'c');
    end
    
    %%%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    %%%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    %%%% Compute Keypoints/Landmark
    %%%%,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    %%%%,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
   
    keyPointsRange = zeros(ll, 1);
    keyPointsSvd = zeros(ll, 1);
    
    %%%%.................................
    %
    % 1. ENDPOINTS from distance jumps
    %
    %%%%.................................
    jump_thrs = 0.5;
    for i = 2:ll-1
        % exctract current range the one before and the next
        rp = Scan.Ranges(i-1);
        rc = Scan.Ranges(i);
        rn = Scan.Ranges(i+1);
    
        % if the first jump is big exor the next is big
        if xor(abs(rc-rp) > jump_thrs, abs(rc-rn) > jump_thrs) 
            keyPointsRange(i) = 1;
        end
    end
    
    
    %%%%...........................................
    %
    % 2. CORNERS from second singular value maxima
    %
    %%%%...........................................
    
    svd_var_thrs = 0.001;
    [~,peaks_idx] = findpeaks(D_profile(:,4),"MinPeakHeight",svd_var_thrs);
    keyPointsSvd(peaks_idx) = 1;
    
    % find indices of the keypoints
    indKey = find(or(keyPointsRange > 0, keyPointsSvd > 0));
    
    % OUTPUT: keypoints index, cartesian coordinates, properties
    % ind_coordinates_properties = ...
    %     [indKey ScanCart(indKey, :) D_profile(indKey, :)];
    keyPoints = {indKey, ScanCart(indKey, :), D_profile(indKey, :)};
    

    % PLOT:
    if pplot == 1
        plotScan(Scan, nplot);
        hold on
        colori = ['k' 'b' 'r' 'm' 'g'];
    
        for i = 1:length(indKey)
            
            ind = indKey(i);
            
            % Plot scan and points in the neigborhood:
            figScan = figure(nplot); plot(ScanCart(ind-Res:ind+Res, 1), ...
                ScanCart(ind-Res:ind+Res, 2), [colori(1 + mod(i,5)) '.'])
            Screen = get(0, 'ScreenSize');
            figScan.Position = [Screen(3)/3 Screen(2) 2*Screen(3)/3 Screen(4)];
            
            % Plot Keypoints:
            plot(ScanCart(ind, 1), ...
                ScanCart(ind, 2), [colori(1 + mod(i,5)) 'x'], ...
                'MarkerSize',20)
    
            figFeatures = figure(nplot+1);
            figFeatures.Position = [Screen(1) Screen(2) Screen(3)/3-20 Screen(4)];
    
            subplot(5,1,1), plot(D_profile(:,1), 'r'), hold on, ...
            if (keyPointsRange(ind) > 0)
                plot(ind, D_profile(ind,1), 'xr'), legend('Curvature')
            else
                plot(ind, D_profile(ind,1), 'xg'), legend('Curvature')
            end
            
            subplot(5,1,2), plot(D_profile(:,2), 'b'), hold on, ...
            if (keyPointsRange(ind) > 0)
                plot(ind, D_profile(ind,2), 'xr'), legend('Curvature')
            else
                plot(ind, D_profile(ind,2), 'xg'), legend('Curvature')
            end
            
            subplot(5,1,3), plot(D_profile(:,3), 'g'), hold on, ...
            if (keyPointsRange(ind) > 0)
                plot(ind, D_profile(ind,3), 'xr'), legend('Curvature')
            else
                plot(ind, D_profile(ind,3), 'xg'), legend('Curvature')
            end
            
            subplot(5,1,4), plot(D_profile(:,4), 'c'), hold on, ...
            if (keyPointsRange(ind) > 0)
                plot(ind, D_profile(ind,4), 'xr'), legend('Curvature')
            else
                plot(ind, D_profile(ind,4), 'xg'), legend('Curvature')
            end
            
            subplot(5,1,5), plot(ScanRange, 'k'), hold on, ...
            if (keyPointsRange(ind) > 0)
                plot(ind, ScanRange(ind), 'xr'), legend('Curvature')
            else
                plot(ind, ScanRange(ind), 'xg'), legend('Curvature')
            end
        end
    end
end
