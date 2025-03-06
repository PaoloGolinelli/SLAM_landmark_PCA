%% Plot Lidar Scan superimposed

function plotScan(lidarScans, numFigure, color)

colori = ['k' 'g' 'r' 'b' 'm'];
figure(numFigure)

for i = 0 : length(lidarScans)-1
    ScanCart = lidarScans(i+1).Cartesian;
    
    plot(ScanCart(:,1),ScanCart(:,2),[colori(1 + mod(i,5)) '.']);
    hold on
end
figure(numFigure)
grid
axis equal
hold off

end
