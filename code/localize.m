function localize(img_dir)
    startup;
    child_relation = cell(4, 1);
    child_relation{1} = [2, 3, 6];              %torso
    child_relation{2} = 4;                      %left upper arm
    child_relation{3} = 5;                      %right upper arm
    child_relation{4} = [];                     %left lower arm
    child_relation{5} = [];                     %right lower arm
    child_relation{6} = [];                     %head

    %A_ij means relation between parent i, child j
    deform_param = ...
      [0, 0, 0, 0, 0, 1; 
       0, 0, 0, 1, 0, 0; 
       0, 0, 0, 0, 1, 0; 
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; %x
    deform_param_y = ...
      [0, 0.3, 0.3, 0, 0, 1;
       0, 0, 0, 1, 0, 0; 
       0, 0, 0, 0, 1, 0; 
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; %y
    deform_param_theta = ...
      [0, 0, 0, 0, 0, 0; 
       0, 0, 0, 0, 0, 0; 
       0, 0, 0, 0, 0, 0; 
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; %theta
    deform_param_scale = ...
      [0, 1, 1, 0, 0, 0.5; 
       0, 0, 0, 1, 0, 0; 
       0, 0, 0, 0, 1, 0; 
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; %scale

    deform_param(:, :, 2) = deform_param_y;
    deform_param(:, :, 3) = 100 * deform_param_theta;
    deform_param(:, :, 4) = 100 * deform_param_scale;

    a = PoseEstimator([4, 5, 2, 3, 6, 1], child_relation, [160, 95, 95, 65, 65, 60], deform_param);
    sticks = a.estimate(img_dir);
    DrawStickman(sticks, imread(fullfile('../buffy_s5e2_original', img_dir)));
end

