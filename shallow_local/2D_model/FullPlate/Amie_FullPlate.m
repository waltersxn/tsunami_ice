%%% modified code from Nurbek Tazhimbetov's PhD thesis, 2022 %%%
%%% changes made by Nestor Walters, advised by Eric Dunham %%%

% The original FullPlate.m is the primary code that time steps the
% simulation. Modified here to allow input of the matrix file name
% with variable mat_file_name


function [q] = Amie_FullPlate(f_1, f_1_t, f_2, ...
                         bc_ice1, bc_ice2, bc_water, bc_water_t, bc_water_NR, ...
                         w0, phi0, order, m, k, iter, ...
                         g, rho_ice, rho_water, gravity, ...
                         nu, BendS, bathy, ice_thick, ...
                         anime, w_exact, phi_exact, ...
                         symmetry_check, time_solver, ...
                         MMS, save_matrix, non_reflective, ...
                         g_ice, time_width, in_angle, video_mode, ...
                         boundary_data, tsunami_data_construct, ...
                         water_block_N, x0, y0, iLU_tol, ...
                         plot_details, q_file_name, mat_file_name, mat_only)

    
    n_ice = g_ice.nPoints;
    n_water = g.nPoints;
    switch MMS
        case 'no'
            nuBendS = nu .* BendS;
            ice1 = multiblock.LaplaceSquared(g_ice, order, 1, nuBendS(1:n_ice));
            ice2 = multiblock.DiFourth(g_ice, order, nu(1:n_ice), BendS(1:n_ice));
            Dice1 = ice1.D;
            Dice2 = ice2.D;
        case 'yes'
            nuBendS = @(x, y) nu(x, y) .* BendS(x, y);
    end    
    
    ice_1 = multiblock.LaplaceSquared(g, order, 1, nuBendS);
    ice_2 = multiblock.DiFourth(g, order, nu, BendS);
    water = multiblock.DiffOp(@scheme.LaplaceCurvilinearMin, g, order, {1, bathy});

    Dice_1 = ice_1.D;
    Dice_2 = ice_2.D;
    Dwater = water.D;
    [Size, ~] = size(Dice_1);

    % ENFORCE FORCING TERMS
    if ~isempty(f_1)
        F1 = @(t) grid.evalOn(g, @(x, y) f_1(x, y, t));
    else
        F1 = @(t) grid.evalOn(g, @(x, y) 0 * x);
    end
    
    if ~isempty(f_1_t)
        F1t = @(t) grid.evalOn(g, @(x, y) - rho_ice * spdiag(ice_thick(x, y)) * f_1_t(x, y, t));
    else
        F1t = @(t) grid.evalOn(g, @(x, y) 0 * x);
    end
    
    if ~isempty(f_2)
        F2 = @(t) grid.evalOn(g, @(x, y) f_2(x, y, t));
    else
        F2 = @(t) grid.evalOn(g, @(x, y) 0 * x);
    end

    % OBTAIN CLOSURE AND PENALTY TERMS
    switch MMS
        case 'yes'
            [closure_ice_1, penalty_ice_1] = scheme.bcSetup(ice_1, bc_ice1);
            [closure_ice_2, penalty_ice_2] = scheme.bcSetup(ice_2, bc_ice2);
            [closure_water, penalty_water] = scheme.bcSetup(water, bc_water);
            N = - (Dice_1 + closure_ice_1 + Dice_2 + closure_ice_2);
        case 'no'
            [closure_ice_1, ~] = scheme.bcSetup(ice1, bc_ice1);
            [closure_ice_2, ~] = scheme.bcSetup(ice2, bc_ice2);
            [closure_water, ~] = scheme.bcSetup(water, bc_water);
            N = - (Dice1 + closure_ice_1 + Dice2 + closure_ice_2);
            N = [N, sparse(n_ice, n_water-n_ice); sparse(n_water-n_ice, n_ice), sparse(n_water-n_ice, n_water-n_ice)];
    end
    
    I_m = speye(Size);
    N_m = 0 * I_m;
    M = - (Dwater + closure_water);

    % PERFOMING SYMMETRY CHECK AND THROWING AN EXCEPTION IN CASE OF
    % NEGATIVE JACOBIAN WITH THE LOCATION OF THE BLOCK IN THE GRID
    switch symmetry_check
        case 'yes'
            block_number = 0;
            for i = 1:water_block_N
                jacobian = diag(water.diffOps{i}.J);
                if nnz(jacobian(jacobian < 0)) ~= 0
                    jacobian
                    block_number = i
                end
            end
    
            if block_number ~= 0
                return
            end            
            
            NN = ice_1.H * N;
            max(max(abs(NN - NN'))) / max(max(abs(NN)))
            switch MMS
                case 'yes'
                    eigs(NN, 3, 'smallestabs')
            end
            MM = ice_1.H * M;
            max(max(abs(MM - MM'))) / max(max(abs(MM)))
    end
    
    switch MMS
        case 'yes'
            ice_thick = grid.evalOn(g, @(x, y) ice_thick(x, y));
        case 'no'
            ice_thick = grid.evalOn(g, @(x, y) ice_thick + x * 0);
    end

    % NON-REFLECTIVE BOUNDARY CONDITIONS USING A GAUSSIAN OR A DATASET
    switch non_reflective
        case 'yes'
            Phi0 = .05;
            time_delay = 4;
            
            switch boundary_data
                case 'gaussian'
                    coord = g.getBoundary(bc_water_NR{1}.boundary);
                    n1     = diag(ice_2.getBoundaryNormal('n1', bc_water_NR{1}.boundary));
                    n2     = diag(ice_2.getBoundaryNormal('n2', bc_water_NR{1}.boundary));
                    inv_velocity = double((gravity * bathy(coord(:, 1), coord(:, 2))) .^ (-.5));
            
                    sx = inv_velocity * cos(in_angle);
                    sy = inv_velocity * sin(in_angle);
                    
                    GaussBorder = @(Z) exp(-.5 * (Z - time_delay * time_width).^2 / time_width^2);
                    dGaussBorder_dZ = @(Z) exp(-.5 * (Z - time_delay * time_width).^2 / time_width^2) .* (-(Z - time_delay * time_width)/time_width^2);
            
                    non_refl.boundary = bc_water_NR{1}.boundary;
                    non_refl.type     = 'N';
                    non_refl.data     = @(t, x, y)     Phi0 * dGaussBorder_dZ(t - sx .* (x - x0) ...
                                                                                - sy .* (y - y0)) .* inv_velocity ...
                                                     - Phi0 * dGaussBorder_dZ(t - sx .* (x - x0) - ...
                                                                                  sy .* (y - y0)) .* sx .* n1 ...
                                                     - Phi0 * dGaussBorder_dZ(t - sx .* (x - x0) - ...
                                                                                  sy .* (y - y0)) .* sy .* n2;
                    non_refl = {non_refl};
                    [~, penalty_nr] = scheme.bcSetup(water, non_refl);
                    penalty_nr_1 = @(t) penalty_nr(t);
                    S_nr = penalty_nr_1;
                    
                case 'tsunami_data'
                    if (exist(tsunami_data_construct.interpolated, 'file') == 2)
                        load(tsunami_data_construct.interpolated)
                    else
                        run(tsunami_data_construct.file)
                    end
            end
            
            [~, penalty_2] = water.boundary_condition(bc_water_NR{1}.boundary, 'N');
            e = water.getBoundaryOperator('e', bc_water_NR{1}.boundary);
            G = grid.evalOn(g, bathy);
            vel_inv = (gravity * abs(G)) .^ (-.5);
            E = vel_inv .* (penalty_2 * e');
            

        case 'no'
            E = N_m;
    end

    A = [I_m, -E;
         N_m, rho_ice * spdiag(ice_thick) * M + rho_water * I_m];
     
    B = [N_m, M;
         N - rho_water * gravity * I_m, N_m];

    A = sparse(A);
    B = sparse(B);

    % FOR LARGE MATRICES SAVE THE MATRIX AND USE A NON-MATLAB SOLVER
    switch save_matrix
        case 'yes'
            q0 = [w0; phi0];
            
            %%% changed for experiment on constancy effects on eigenmodes
            
            if (exist(mat_file_name,'file') == 2)
                delete(mat_file_name);
            end
            save(mat_file_name,'A', 'B', 'q0');
            if mat_only == 1
                return
            end
            %%% previous version:
            % matname = 'mat_ABq_bvar_tvar.mat';
            % if (exist(matname,'file') == 2)
            %     delete(matname);
            % end
            % save(matname,'A', 'B', 'q0');

            %%% original save matrix cases
%              if (exist('mat_A_B_q.mat', 'file') == 2)
%                  delete('mat_A_B_q.mat');
%              end
%              save('mat_A_B_q.mat', 'A', 'B', 'q0');
            
            q = q0;
            % return
    end

    switch MMS
        case 'yes'
            [~, penalty_water_t] = scheme.bcSetup(water, bc_water_t);
    
            penalty_ice_1 = @(t) - penalty_ice_1(t);
            penalty_ice_2 = @(t) - penalty_ice_2(t);
            penalty_water_t = @(t) rho_ice * spdiag(ice_thick) * penalty_water_t(t);
            data_funcs_ice = cell(5, 1);
            data_funcs_ice{1} = F2;
            data_funcs_ice{2} = penalty_ice_1;
            data_funcs_ice{3} = penalty_ice_2;
            data_funcs_ice{4} = F1t;
            data_funcs_ice{5} = penalty_water_t;
            S_ice = [];
    
            for i = 1:numel(data_funcs_ice)
                data = data_funcs_ice{i};
                if ~isempty(data)
                    if isempty(S_ice)
                        S_ice = data;
                    else
                        S_ice = @(t) S_ice(t) + data(t);
                    end
                end
            end
    
            penalty_water = @(t) - penalty_water(t);
            data_funcs_water = cell(2, 1);
            data_funcs_water{1} = F1;
            data_funcs_water{2} = penalty_water;
            S_water = [];
            for i = 1:numel(data_funcs_water)
                data = data_funcs_water{i};
                if ~isempty(data)
                    if isempty(S_water)
                        S_water = data;
                    else
                        S_water = @(t) S_water(t) + data(t);
                    end
                end
            end
    end
    
    % TIME SOLVER; IMPLICIT SOLVER USES CRANK-NICOLSON AND EXPLICIT SOLVER
    % USER 4TH ORDER RUNGE-KUTTA
    switch time_solver
        case 'implicit'
            switch MMS
                case 'yes'
                    S = @(t) [S_water(t); S_ice(t)];
            end
            
            % FORCING IN THE TIME SOLVER CHANGES BASED ON THE TYPE OF
            % PROBLEM WE ARE SOLVING
            switch non_reflective
                case 'yes'
                    switch boundary_data
                        case 'gaussian'
                            S = @(t) [S_nr(t); 0 * S_nr(t)];
                        case 'tsunami_data'
                            S = [S_nr; 0 * S_nr];
                        case 'none'
                            S = 0 * [w0; phi0];
                    end
            end
            q = [w0; phi0];
            
            
%             fileID = fopen(q_file_name, 'w');
%             fprintf(fileID, '%.7f ', q);
%             fprintf(fileID, '\n');

            fileID = fopen(q_file_name, 'w');
            fwrite(fileID, q, 'single');
            fclose(fileID);
            
            
            w_abs = abs(w0);
            LHS = (A - k/2 * B);
            
            % SETTING THE OPTIONS FOR INCOMPLETE LU DECOMPOSITION
            options.type = 'crout';
            options.milu = 'row';
            options.droptol = iLU_tol; 
            [iL, iU] = ilu(LHS, options);
            RHS = (A + k/2 * B);
            
            switch anime
                case 'yes'
                    switch video_mode
                        case 'yes'
                            myVideo = VideoWriter('myVideoTsunamiData', 'MPEG-4');
                            myVideo.FrameRate = 15;
                            open(myVideo)
                    end
                    figure()
                    hold on
                    cmap % leave this unchanged
                    for i = 1:length(g_ice.boundaryGroups.FR)
                        free_coord = g.getBoundary(g_ice.boundaryGroups.FR(i));
                        p1 = plot3(free_coord(:, 1), free_coord(:, 2), 1000 + 0*free_coord(:, 1), 'black-', 'LineWidth', 1.5);
                        uistack(p1, 'top')                        
                        hold on
                    end
                    for i = 1:length(g.boundaryGroups.R)
                        reflectN = g.getBoundary(g.boundaryGroups.R(i));
                        p2 = plot3(reflectN(:, 1), reflectN(:, 2), 1000 + 0*reflectN(:, 1), 'black-', 'LineWidth', 0.5);
                        uistack(p2, 'top')
                        hold on
                    end

                    xlabel('$x$ (km)', 'Interpreter', 'latex')
                    ylabel('$y$ (km)', 'Interpreter', 'latex')
                    xlim([plot_details.left, plot_details.right])
                    ylim([plot_details.bottom, plot_details.top])

                    ph = multiblock.Surface(g, real(q(1:Size)));
                    caxis([plot_details.caxis_low, plot_details.caxis_up])

                    c = colorbar;
                    c.Label.String = '$w$ in meters';
                    c.Label.Interpreter = 'latex';
                    c.Label.FontSize = 18;
                    ph.set('LineStyle', 'none');
                    ph.set('FaceColor', 'interp');
                    count = 1;
                    while count < iter + 1
                        count = count + 1
                        switch MMS
                            case 'yes'
                                w_val = grid.evalOn(g, @(x, y) w_exact((count-1) * k, x, y));
                                phi_val = grid.evalOn(g, @(x, y) phi_exact((count-1) * k, x, y));
                                b = RHS * q + k / 2 * (S((count - 2) * k) + S((count - 1) * k));
                            case 'no'
                                switch non_reflective
                                    case 'yes'
                                        switch boundary_data
                                            case 'gaussian'
                                                b = RHS * q + k / 2 * (S((count - 2) * k) + S((count - 1) * k));
                                            case 'tsunami_data'
                                                if count < iter
                                                    b = RHS * q + k / 2 * (S(:, count) + S(:, count - 1));
                                                else
                                                    b = RHS * q;
                                                end
                                            case 'none'
                                                b = RHS * q;
                                        end
                                        (count - 1) * k
                                    case 'no'
                                        b = RHS * q;
                                end
                        end
                        q = gmres(LHS, b, [], 1e-9, 300, iL, iU, q);
                        
%                         fileID = fopen(q_file_name, 'a');
%                         fprintf(fileID, '%.7f ', q);
%                         fprintf(fileID, '\n');

                        fileID = fopen(q_file_name, 'a');
                        fwrite(fileID, q, 'single');
                        fclose(fileID);

                        
                        ph.CData = q(1:Size);
                        ph.ZData = q(1:Size);
                        switch video_mode
                            case 'no'
                                % drawnow;
                                pause(0.01)
                        end
                        legend([p1, p2], {'Ice shelf front', 'Wall'}, 'Interpreter', 'latex', 'Location', 'southeast')
                        
                        switch video_mode
                            case 'yes'
                                pause(0.01)
                                frame = getframe(gcf);
                                writeVideo(myVideo, frame);
                        end
                    end
                    switch video_mode
                        case 'yes'
                            close(myVideo)
                    end
                case 'no'
                    count = 1;
                    while count < iter + 1
                        count = count + 1;
                        switch MMS
                            case 'yes'
                                b = RHS * q + k / 2 * (S((count - 1) * k) + S((count - 2) * k));
                            case 'no'
                                switch non_reflective
                                    case 'yes'
                                        switch boundary_data
                                            case 'gaussian'
                                                b = RHS * q + k / 2 * (S((count - 2) * k) + S((count - 1) * k));
                                            case 'tsunami_data'
                                                if count < iter
                                                    b = RHS * q + k / 2 * (S(:, count) + S(:, count - 1));
                                                else
                                                    b = RHS * q;
                                                end
                                            case 'none'
                                                b = RHS * q;
                                        end
                                        (count - 1) * k
                                    case 'no'
                                        b = RHS * q;
                                end
                        end
                        q = gmres(LHS, b, [], 1e-9, 300, iL, iU, q);
%                         fileID = fopen(q_file_name, 'a');
%                         fprintf(fileID, '%.7f ', q);
%                         fprintf(fileID, '\n');
                        fileID = fopen(q_file_name, 'a');
                        fwrite(fileID, q, 'single');
                        fclose(fileID);
                    end
            end

        case 'explicit'
            D = A \ B;
            S = @(t) A \ [S_water(t); S_ice(t)];
            q0 = [w0; phi0];
            ts = time.Rungekutta4(D, S, k, 0, q0);
            switch anime
                case 'yes'
                    figure()
                    ph = multiblock.Surface(g, real(q0(1:Size)));
                    ph.set('LineStyle', 'none');
                    ph.set('FaceColor', 'interp');
                    colorbar
                    count = 1;
                    while count < iter + 1
                        count = count + 1
                        [q, t] = ts.getV();
                        w_val = grid.evalOn(g, @(x, y) w_exact((count-1) * k, x, y));
                        ph.CData = abs(w_val - q(1:Size));
                        ph.ZData = abs(w_val - q(1:Size));
                        ts.step();
                        drawnow;
                    end
                case 'no'
                    ts.stepN(iter);
                    q = ts.getV();
            end
    end

    q = q(1:Size);
end