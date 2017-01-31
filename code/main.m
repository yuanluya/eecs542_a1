startup;
ideal_parameters = cell(4, 1);
ideal_parameters{1} = [0, 80, 0, 1];        %torso
ideal_parameters{2} = [-107.5, 20, 0, 1];    %left upper arm
ideal_parameters{3} = [107.5, 20, 0, 1];     %right upper arm
ideal_parameters{4} = [0, -45, 0, 1];       %head
%ideal_parameters{5} = [0, -45, 0, 1];       %head
%ideal_parameters{6} = [0, -45, 0, 1];       %head

child_relation = cell(4, 1);
child_relation{1} = [2, 3, 4];              %torso
child_relation{2} = [];                     %left upper arm
child_relation{3} = [];                     %right upper arm
child_relation{4} = [];                     %head
%child_relation{5} = [];                     %head
%child_relation{6} = [];                     %head

%A_ij means relation between parent i, child j
deform_param = ...
  [1, 1, 1, 1; 
   0, 0, 0, 0; 
   0, 0, 0, 0; 
   0, 0, 0, 0]; %x
deform_param_y = ...
  [1, 1, 1, 1; 
   0, 0, 0, 0; 
   0, 0, 0, 0; 
   0, 0, 0, 0]; %y
deform_param_theta = ...
  [0.5, 0.5, 0.5, 1.0;
   0, 0, 0, 0; 
   0, 0, 0, 0; 
   0, 0, 0, 0]; %theta
deform_param_scale = ...
  [1, 1, 1, 1;
   0, 0, 0, 0; 
   0, 0, 0, 0; 
   0, 0, 0, 0]; %scale

deform_param(:, :, 2) = deform_param_y;
deform_param(:, :, 3) = 100 * deform_param_theta;
deform_param(:, :, 4) = 100 * deform_param_scale;

a = PoseEstimator(ideal_parameters, [2, 3, 4, 1], child_relation, 0 * deform_param);
proposed_parts = a.estimate('058918.jpg')
[parts, ~] = size(proposed_parts);
sticks = zeros(4, parts);
for i = 1: parts
    sticks(:, i) = a.changeBase(proposed_parts(i,:), i);
end
sticks
stick_hdl = DrawStickman(sticks, imread('../buffy_s5e2_original/000063.jpg'));
