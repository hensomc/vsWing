function draw_panel_fib(plot_flag, panXY, lamDATA, numFibPts, refPath, panFibAxis)

% panXY = panel corner point coordinates

if(plot_flag > 0)
    %% Fiber Paths
    % Draw a panel in xy & st system
    figure;ISO_PlotQ4_Sub(panXY)
    
    % Draw ply courses in xy and st system
    if sum(abs(lamDATA.theta0-lamDATA.theta1))>0
        draw_ply_courses(panXY, lamDATA, 2, numFibPts, refPath, panFibAxis);
    end
    
    % Save to image file
    print('fiber_paths','-dpng');
       
end