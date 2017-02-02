classdef PoseEstimator < handle
    
     properties (GetAccess = public, SetAccess = public)
        
        num_parts = 4
        num_x_buckets = 20
        num_y_buckets = 20
        num_theta_buckets = 10
        num_scale_buckets = 5
        model_len = [160, 95, 95, 65, 65, 60];
        min_scale = 0.5
        max_scale = 1.5
        min_theta = -pi / 2
        max_theta = pi / 2
        
        %[x, y, theta, scale], scale: [0.5, 1.5]
        ideal_parameters
        step_size
        search_step
        momentum = [-1 * ones(4, 1), ones(4, 1)]
        momentum_step = [-4, -4, -3, -2]
        
        %define the order to calculate energy
        table_set_order
        child_relation
        parent_relation
        energy_map
        all_combos
        match_cost_cache
        change_base_cache
        last_optimal
        
        %need to be tuned
        %[variable X partNum X partNum]
        deform_cost_weights
        random_radius = 0
        match_cost_weights = 5e-2
        
        %define energy functions
        match_cost = @match_energy_cost
        
        %image info
        image_dir = '../buffy_s5e2_original/';
        all_names
        seq
        img_height
        img_width
        lF
        
     end
     
     methods (Access = public)
         function obj = PoseEstimator(ideal_parameters, table_set_order, child_relation, deform_cost_weights)
            
            all_names = dir(obj.image_dir);
            obj.all_names = containers.Map();
            obj.lF = ReadStickmenAnnotationTxt('../data/buffy_s5e2_sticks.txt');
            counter = 0;
            for n = 1: numel(all_names)
                if length(all_names(n).name) ~= 10
                    continue
                else
                    counter = counter + 1;
                    obj.all_names(all_names(n).name) = counter;
                end
            end
            obj.table_set_order = table_set_order;
            obj.child_relation = child_relation;
            obj.parent_relation = cell(numel(obj.child_relation), 1);
            %get parent order             
            for i = 1: numel(obj.child_relation)
                for c = 1: numel(obj.child_relation{i})
                    obj.parent_relation{obj.child_relation{i}(c)} = ...
                        [obj.parent_relation{obj.child_relation{i}(c)}, i];
                end
            end
            obj.num_parts = numel(ideal_parameters);
            obj.deform_cost_weights = zeros(4, obj.num_parts, obj.num_parts);
            
            obj.deform_cost_weights = deform_cost_weights;
            
            obj.ideal_parameters = ideal_parameters;
            for j  = 1: 4
                obj.ideal_parameters{j} = zeros(length(ideal_parameters{j}));
                if j ~= 3
                    for p = 1: numel(child_relation)
                        for c = 1: numel(child_relation{p})
                            if j == 4 %scale
                                obj.ideal_parameters{j}(p, child_relation{p}(c)) = ...
                                    ideal_parameters{child_relation{p}(c)}(j) / ideal_parameters{p}(j);
                            else
                                obj.ideal_parameters{j}(p, child_relation{p}(c)) = ...
                                    ideal_parameters{child_relation{p}(c)}(j) - ideal_parameters{p}(j);
                            end
                        end
                    end
                end
            end
            
            %cell array, one cell for a part, 
            %a matrix :[min_energy, optimal_location_idx (in obj.all_combos)]
            obj.energy_map = cell(numel(obj.ideal_parameters), 1);
            %initialize match cost cache
            obj.match_cost_cache = cell(numel(obj.ideal_parameters),1);
            
            search_step = [1, obj.num_x_buckets, ...
                            obj.num_y_buckets * obj.num_x_buckets, ...
                            obj.num_theta_buckets * obj.num_y_buckets * obj.num_x_buckets];
            obj.search_step = [-1 * search_step, search_step];
         end
         
         function cost = deformCost(obj, part_p, part_c, lp_idx, lc_idx)
            lp = obj.all_combos(lp_idx, :);
            lc = obj.all_combos(lc_idx, :);
            coor_p = obj.change_base_cache{part_p}(lp_idx, :);
            coor_c = obj.change_base_cache{part_c}(lc_idx, :);
            if isnan(coor_p(1))
                coor_p = obj.changeBase(lp_idx, part_p);
                obj.change_base_cache{part_p}(lp_idx, :) = coor_p;
            end
            if isnan(coor_c(1))
                coor_c = obj.changeBase(lc_idx, part_c);
                obj.change_base_cache{part_c}(lc_idx, :) = coor_c;
            end
            
            coor_C = [coor_c; [coor_c(3: 4), coor_c(1: 2)]];
            dists = coor_C - repmat(coor_p, [2, 1]);
            dists = [dists(:, [1, 2]); dists(:, [3, 4])];
            [~, I] = min(sum(dists .^ 2, 2));
            diff_junct = dists(I, :);
            
            x_diff = obj.deform_cost_weights(part_p, part_c, 1) * ...
                abs(diff_junct(1));
            y_diff = obj.deform_cost_weights(part_p, part_c, 2) * ...
                abs(diff_junct(2));
            theta_diff = obj.deform_cost_weights(part_p, part_c, 3) * ...
                abs(lc(3) - lp(3));
            scale_diff = obj.deform_cost_weights(part_p, part_c, 4) * ...
                abs(log(lc(4)) - log(lp(4)));
            cost = x_diff + y_diff + theta_diff + scale_diff;
         end
         
         function energy = calcEnergy(obj, self_part_idx, l_self_idx, ...
                                      parent_part_idx, l_parent_idx)
             %check invalid location
             if l_self_idx <= 0 || l_self_idx > size(obj.all_combos, 1)
                energy = inf;
                return;
             end
             
             %pairwise
             if parent_part_idx
                 pair_wise_energy = obj.deformCost(parent_part_idx, ...
                                                  self_part_idx, ...
                                                  l_parent_idx, ...
                                                  l_self_idx);
             else
                 pair_wise_energy = 0;
             end
             
             %children
             children_energy = 0;
             for c = 1: numel(obj.child_relation{self_part_idx})
                assert(sum(isnan(obj.energy_map{obj.child_relation{self_part_idx}(c)}(l_self_idx, :))) == 0, ...
                    'Havent done children parts');
                temp_energy = obj.energy_map{obj.child_relation{self_part_idx}(c)}(l_self_idx, 1);
                children_energy = children_energy + temp_energy;
             end
             
             %match
             if(isnan(obj.match_cost_cache{self_part_idx}(l_self_idx)))
                 match_energy = obj.match_cost_weights * ...
                     match_energy_cost(obj.all_combos(l_self_idx, :), self_part_idx, obj.seq, obj.lF);
                 obj.match_cost_cache{self_part_idx}(l_self_idx) = match_energy;
             else
                match_energy = obj.match_cost_cache{self_part_idx}(l_self_idx);
             end
             %total
             %fprintf('match cost: %f\tdeform cost: %f\n', match_energy, pair_wise_energy);
             energy = pair_wise_energy + match_energy + children_energy;
         end
         
         function [current_min_energy, current_min_idx] = localMin(obj, self_part_idx, ...
                 parent_part_idx, l_parent_idx)            
            
            if isnan(obj.last_optimal)
                current_min_idx = randi(size(obj.all_combos, 1), 1);
            else
                current_min_idx = obj.last_optimal;
            end
            
            current_min_energy = obj.calcEnergy(self_part_idx, current_min_idx, ...
                                                parent_part_idx, l_parent_idx);
            obj.momentum = ones(1, 8);
            while true
                all_neighbors_idx = current_min_idx * ones(1, 8) + obj.search_step .* obj.momentum;
                energies = zeros(9, 1);
                for i = 2: 9
                    energies(i) = obj.calcEnergy(...
                        self_part_idx, all_neighbors_idx(i - 1), parent_part_idx, l_parent_idx);
                end
                energies(1) = current_min_energy; %min return first element when equal
                [current_min_energy, best_idx] = min(energies);
                
                if best_idx > 1
                    current_min_idx = all_neighbors_idx(best_idx - 1);
                    obj.momentum = ones(1, 8);
                    wild_try = mod(best_idx + 4, 8);
                    wild_try = wild_try + (wild_try == 0) * 8;
                    m_id = mod(best_idx, 4);
                    m_id = m_id + (m_id == 0) * 4;
                    obj.momentum(wild_try) = obj.momentum_step(m_id);
                else
                    obj.last_optimal = current_min_idx;
                    return;
                end
                
            end            
         end
         
         function updateEnergymap(obj, part_idx)
            obj.last_optimal = nan;
            
            if part_idx == obj.table_set_order(end)
                for j = 1: size(obj.all_combos, 1)
                    if mod(j, 100) == 0
                        fprintf('Part: %d, possiblility %d/%d\n', part_idx, j, ...
                            obj.num_x_buckets * obj.num_y_buckets * obj.num_theta_buckets * obj.num_scale_buckets);
                    end
                    total = obj.calcEnergy(part_idx, j, [], []);
                    obj.energy_map{part_idx}(j, :) = [total, nan];
                end
            else
                for i = 1: size(obj.all_combos, 1)
                    if mod(i, 100) == 0
                        fprintf('Part: %d, possiblility %d/%d\n', part_idx, i, ...
                            obj.num_x_buckets * obj.num_y_buckets * obj.num_theta_buckets * obj.num_scale_buckets);                
                    end
                    [temp_min_energy, temp_min_idx] = ...
                        obj.localMin(part_idx, obj.parent_relation{part_idx}, i);                
                    obj.energy_map{part_idx}(i, :) = ...
                        [temp_min_energy, temp_min_idx];
                end
            end
            
         end
         
         function parts = estimate(obj, seq)
            obj.seq = obj.all_names(seq);
            img = imread(fullfile(obj.image_dir, seq));
            [obj.img_height, obj.img_width, ~] = size(img);
            obj.step_size = [floor(obj.img_width / obj.num_x_buckets), ...
                             floor(obj.img_height / obj.num_y_buckets), ...
                             pi / obj.num_theta_buckets, 1 / obj.num_scale_buckets];
            
            %set up all_combo
            xs = (obj.step_size(1) / 2): obj.step_size(1): obj.img_width;
            ys = (obj.step_size(2) / 2): obj.step_size(2): obj.img_height;
            thetas = (-pi / 2 + (obj.step_size(3) / 2)): obj.step_size(3): pi / 2;
            scales = 0.5 + obj.step_size(4) / 2: obj.step_size(4): 1.5;
            obj.all_combos = combvec(xs(1: obj.num_x_buckets), ...
                                 ys(1: obj.num_y_buckets), ...
                                 thetas(1: obj.num_theta_buckets), ...
                                 scales(1: obj.num_scale_buckets)).';
            
            for i = 1: numel(obj.ideal_parameters)
                obj.energy_map{i} = nan(size(obj.all_combos, 1), 2);
                obj.match_cost_cache{i} = nan(size(obj.all_combos, 1), 1);
                obj.change_base_cache{i} = nan(size(obj.all_combos, 1), 4);
            end
            
            %forward calculate energies
            for i = 1: numel(obj.table_set_order)
                obj.updateEnergymap(obj.table_set_order(i));
            end
            %backward return optimal values for parts
            parts = zeros(obj.num_parts, 4);
            parts_Lidx = zeros(obj.num_parts, 1);
            for j = numel(obj.table_set_order): -1: 1
                % not backtrace yet
                if sum(parts(obj.table_set_order(j), :) == 0) == 4
                    %not root
                    if obj.parent_relation{obj.table_set_order(j)}
                        temp = obj.energy_map{obj.table_set_order(j)} ...
                            (parts_Lidx(obj.parent_relation{obj.table_set_order(j)}), :);
                        parts_Lidx(obj.table_set_order(j)) = temp(2);
                        parts(obj.table_set_order(j), :) = ...
                            obj.changeBase(temp(2), obj.table_set_order(j));
                    else%root
                        vals = obj.energy_map{obj.table_set_order(j)};
                        [~, idx] = min(vals(:, 1));
                        parts_Lidx(obj.table_set_order(j)) = idx;
                        parts(obj.table_set_order(j), :) = ...
                            obj.changeBase(idx, obj.table_set_order(j));
                    end
                end
                %find the child of this part
                for c = 1: numel(obj.child_relation{obj.table_set_order(j)})
                    child_part_idx = obj.child_relation{obj.table_set_order(j)}(c);
                    temp = obj.energy_map{child_part_idx} ...
                            (parts_Lidx(obj.table_set_order(j)), :);
                    parts_Lidx(child_part_idx) = temp(2);
                    parts(child_part_idx, :) = ...
                        obj.changeBase(temp(2), child_part_idx);
                end
            end
            parts = parts.';
         end
         
         %coor is in format [x1, y1, x2, y2]
         function coor = changeBase(obj, location_idx, part_idx) 
            location = obj.all_combos(location_idx, :);
            stick_len = location(4) * obj.model_len(part_idx);
            if part_idx == 2 || part_idx == 3
                dx = 0.5 * stick_len * cos(location(3));
                dy = 0.5 * stick_len * sin(location(3));
            else
                dx = 0.5 * stick_len * sin(location(3));
                dy = 0.5 * stick_len * cos(location(3));
            end
            coor = [location(1), location(2), location(1), location(2)] ...
                 + [dx, dy, -dx, -dy];
         end
         
     end
end