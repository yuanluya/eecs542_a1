classdef PoseEstimator < handle
    
     properties (GetAccess = public, SetAccess = public)
        
        num_parts = 4
        num_x_buckets = 25
        num_y_buckets = 25
        num_theta_buckets = 15
        num_scale_buckets = 5
        %model_len = [160, 95, 95, 65, 65, 60];
        model_len = [160, 95, 95, 60];
        min_scale = 0.5
        max_scale = 1.5
        min_theta = -pi / 2
        max_theta = pi / 2
        
        %[x, y, theta, scale], scale: [0.5, 1.5]
        ideal_parameters
        step_size
        
        %define the order to calculate energy
        table_set_order
        child_relation
        parent_relation
        energy_map
        match_cost_cache
        last_optimal
        
        %need to be tuned
        %[variable X partNum X partNum]
        deform_cost_weights
        random_init_radius = [-0, 0]
        match_cost_weights = 1e-1
        
        
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
            
            %cell array, one cell for a part, a long map, 500000
            %Mapping: mat2str(l_parent)->[min_energy, best_location]
            obj.energy_map = cell(numel(obj.ideal_parameters), 1);
            for i = 1: numel(obj.ideal_parameters)
                obj.energy_map{i} = containers.Map();
            end
            
            %initialize match cost cache
            obj.match_cost_cache = cell(numel(obj.ideal_parameters),1);
            for j = 1:numel(obj.match_cost_cache)
                obj.match_cost_cache{j} = containers.Map();
            end
         end
         
         function cost = deformCost(obj, part_p, part_c, lp, lc)
            coor_p = obj.changeBase(lp, part_p);
            coor_c = obj.changeBase(lc, part_c); %x1,x2,y1,y2
            
            coor_C = [coor_c; [coor_c(3: 4), coor_c(1: 2)]];
            dists = bsxfun(@minus, coor_C, coor_p);
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
         
         %return true if not in image
         function in_or_not = checkInPicture(obj, l_self)
             in_or_not = ~(sum(l_self([1, 2]) > 0) ~= 2 ...
                        || l_self(3) < obj.min_theta ...
                        || l_self(1) > obj.img_width ...
                        || l_self(2) > obj.img_height ...
                        || l_self(3) > obj.max_theta ...
                        || l_self(4) > obj.max_scale...
                        || l_self(4) < obj.min_scale);
         end
         
         function energy = calcEnergy(obj, self_part_idx, l_self, ...
                                      parent_part_idx, l_parent)
             %check invalid location
             if ~obj.checkInPicture(l_self)
                energy = inf;
                return;
             end
             %pairwise
             if parent_part_idx
                 pair_wise_energy = obj.deformCost(parent_part_idx, ...
                                                  self_part_idx, ...
                                                  l_parent, ...
                                                  l_self);
             else
                 pair_wise_energy = 0;
             end
             
             %children
             children_energy = 0;
             for c = 1: numel(obj.child_relation{self_part_idx})
                assert(isKey(obj.energy_map{obj.child_relation{self_part_idx}(c)}, ...
                    mat2str(l_self)), 'Havent done children parts');
                energy_pair = obj.energy_map{obj.child_relation{self_part_idx}(c)}(mat2str(l_self));
                children_energy = children_energy + energy_pair(1);
             end
             
             %match
             if(~isKey(obj.match_cost_cache{self_part_idx},mat2str(l_self)))
                 fixed_part_id = self_part_idx;
                 if fixed_part_id == 4
                     fixed_part_id = 6;
                 end
                 match_energy = obj.match_cost_weights * match_energy_cost(l_self, fixed_part_id, obj.seq, obj.lF);
                 obj.match_cost_cache{self_part_idx}(mat2str(l_self)) = match_energy;
             else
                match_energy = obj.match_cost_cache{self_part_idx}(mat2str(l_self));
             end
             %total
             %fprintf('match cost: %f\tdeform cost: %f\n', match_energy, pair_wise_energy);
             energy = pair_wise_energy + match_energy + children_energy;
         end
         
         
         function real_child = sampleFromParent(obj, self_part_idx, parent_part_idx, l_parent)
             %random sample from optimal_child
             real_child = optimal_child;
             return
             %... + [randi(obj.random_init_radius), randi(obj.random_init_radius), 0, 0];
             if ~obj.checkInPicture(real_child)
                real_child = l_parent;
             end
         end
         
         function [current_min_energy, current_min] = localMin(obj, self_part_idx, parent_part_idx, l_parent)            
            
            if isnan(obj.last_optimal)
                init_idx = [randi([1, obj.num_x_buckets]), ...
                            randi([1, obj.num_y_buckets]), ...
                            floor(obj.num_theta_buckets / 2), floor(obj.num_scale_buckets / 2)];
                current_min = [0, 0, obj.min_theta, obj.min_scale] ...
                    + 0.5 * obj.step_size + (init_idx - 1) .* obj.step_size;
            else
                current_min = obj.last_optimal;
            end
            
            %current_min = obj.sampleFromParent(self_part_idx, parent_part_idx, l_parent);
            %current_min = l_parent;
            current_min_energy = obj.calcEnergy(self_part_idx, current_min, parent_part_idx, l_parent);    
            while true
                neighbors1 = repmat(current_min, [4, 1]) - eye(4) .* diag(obj.step_size);
                neighbors2 = repmat(current_min, [4, 1]) + eye(4) .* diag(obj.step_size);
                all_neighbors = [neighbors1; neighbors2];
                energies = zeros(9, 1);
                for i = 2: 9
                    energies(i) = obj.calcEnergy(...
                        self_part_idx, all_neighbors(i - 1, :), parent_part_idx, l_parent);
                end
                energies(1) = current_min_energy; %min return first element when equal
                [current_min_energy, best_idx] = min(energies);
                
                if best_idx > 1
                    current_min = all_neighbors(best_idx - 1, :);
                else
                    obj.last_optimal = current_min;
                    %l_parent
                    return;
                end
                
            end            
         end
         
         function updateEnergymap(obj, part_idx)
            obj.last_optimal = nan;
            xs = (obj.step_size(1) / 2): obj.step_size(1): obj.img_width;
            ys = (obj.step_size(2) / 2): obj.step_size(2): obj.img_height;
            thetas = (-pi / 2 + (obj.step_size(3) / 2)): obj.step_size(3): pi / 2;
            scales = 0.5: obj.step_size(4): 1.5;
            all_combos = combvec(xs(1: obj.num_x_buckets), ...
                                 ys(1: obj.num_y_buckets), ...
                                 thetas(1: obj.num_theta_buckets), ...
                                 scales(1: obj.num_scale_buckets)).';
            if part_idx == obj.table_set_order(end)
                for j = 1: size(all_combos, 1)
                    if mod(j, 50) == 0
                        fprintf('Part: %d, possiblility %d/%d\n', part_idx, j, ...
                            obj.num_x_buckets * obj.num_y_buckets * obj.num_theta_buckets * obj.num_scale_buckets);
                    end
                    total = obj.calcEnergy(part_idx, all_combos(j, :), [], []);
                    obj.energy_map{part_idx}(mat2str(all_combos(j, :))) = total;
                end
            else
                for i = 1: size(all_combos, 1)
                    if mod(i, 50) == 0
                        fprintf('Part: %d, possiblility %d/%d\n', part_idx, i, ...
                            obj.num_x_buckets * obj.num_y_buckets * obj.num_theta_buckets * obj.num_scale_buckets);                
                    end
                    [temp_min_energy, temp_min] = ...
                        obj.localMin(part_idx, obj.parent_relation{part_idx}, all_combos(i, :));                
                    obj.energy_map{part_idx}(mat2str(all_combos(i, :))) = ...
                        [temp_min_energy, temp_min];
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
            
            %forward calculate energies
            for i = 1: numel(obj.table_set_order)
                obj.updateEnergymap(obj.table_set_order(i));
            end
            %backward return optimal values for parts
            parts = zeros(obj.num_parts, 4);
            for j = numel(obj.table_set_order): -1: 1
                % not backtrace yet
                if sum(parts(obj.table_set_order(j), :) == 0) == 4
                    %not root
                    if obj.parent_relation{obj.table_set_order(j)}
                        temp = obj.energy_map{obj.table_set_order(j)} ...
                            (mat2str(parts(obj.parent_relation{obj.table_set_order(j)}, :)));
                        parts(obj.table_set_order(j), :) = temp(2: 5);
                    else%root
                        vals = cell2mat(values(obj.energy_map{obj.table_set_order(j)}));
                        [~, idx] = min(vals);
                        all_keys = keys(obj.energy_map{obj.table_set_order(j)});
                        parts(obj.table_set_order(j), :) = eval(all_keys{idx});
                    end
                end
                %find the child of this part
                for c = 1: numel(obj.child_relation{obj.table_set_order(j)})
                    child_part_idx = obj.child_relation{obj.table_set_order(j)}(c);
                    temp = obj.energy_map{child_part_idx} ...
                            (mat2str(parts(obj.table_set_order(j), :)));
                    parts(child_part_idx, :) = temp(2: 5);
                end
            end
         end
         
         %coor is in format [x1,x2,y1,y2]
         function coor = changeBase(obj, location, part_idx) 
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
         
         function reset(obj)
            obj.energy_map = cell(numel(obj.ideal_parameters), 1);
            for i = 1: numel(obj.ideal_parameters)
                obj.energy_map{i} = containers.Map('ValueType', 'any');
                obj.match_cost_cache{i} = containers.Map('ValueType', 'any');
            end
         end
     end
end