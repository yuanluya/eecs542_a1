ideal_parameters = cell(4, 1);
ideal_parameters{1} = [0, 80, 0, 1];        %torso
ideal_parameters{2} = [-47.5, 20, 0, 1];    %left upper arm
ideal_parameters{3} = [47.5, 20, 0, 1];     %right upper arm
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

deform_param = [0.7; 0.2; 0.6; 0.4];
deform_param = repmat(deform_param, [1, 4, 4]);
deform_param = permute(deform_param, [3, 2, 1]);

a = PoseEstimator(ideal_parameters, [2, 3, 4, 1], child_relation, deform_param);
a.estimate('000063.jpg');