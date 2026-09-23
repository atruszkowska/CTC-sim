function compute_area_ratio

    % Fold change of cell area
    
    % Cell radius 
    R_cell = 1.01e-5;

    % Number of cells
    N_prt = 10;
    % Number of collected time steps
    N_times = 60;
    % Number of repetitions
    N_rep = 1;

    % Results (Repetition # | Ai/A0 for each i)
    area_change = zeros(N_rep, N_times);

    % Dataset
    data = load('no_bonds_results.txt');

    for ir = 1:N_rep
        % Datast portion for this repetition
        idx0 = (ir-1)*N_prt*N_times + 1;
        stride = N_prt;
        A0 = get_area(data(idx0:idx0 + stride - 1,1), data(idx0:idx0 + stride - 1,2), R_cell);
        area_change(ir, 1) = 1;

        % Compute all remaining areas, store each ratio
        for jt = 2:N_times
            idx = (ir-1)*N_prt*N_times + 1 + (jt-1)*N_prt;
            % Compute the current area
            A = get_area(data(idx:idx+stride-1,1), data(idx:idx+stride-1,2), R_cell);
            area_change(ir, jt) = A/A0;
        end
    end
   
    % Write area change to file
    save('area_change_mat.txt', 'area_change', '-ascii');

    % Compute the mean area 
    mean_area = mean(area_change, 1);
    save('mean_area_change_mat.txt', 'mean_area', '-ascii');
    for jt = 1:N_times
        fprintf('Time step %d, Mean Area Change: %f\n', jt, mean_area(jt));
    end
    
    save('area_res')
end

function A = get_area(x, y, R_cell)
    % Compute the area spanned by cells using convex hull - if the cells
    % are aligned, perturb their positions slightly to allow the
    % computation - this does not change the actual positions
    
    while true
        try 
            [k, A] = convhull(x, y);
            break;
        catch ME
            if x > 0 & y > 0
                x = x.*(1+rand(size(x))/100);
                y = y.*(1+rand(size(x))/100);
                continue;
            else
                x = rand(size(x))*R_cell;
                y = rand(size(x))*R_cell;
                continue;
            end
        end
    end

end
