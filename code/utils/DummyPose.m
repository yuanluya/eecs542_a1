function coor =  DummyPose(img,det_bb)
% dummy routine generating stickman (it always generate the same one just scaled and shifted according to the bounding box)
% Input 
%   img - input image
%   det_bb - detection bouding box [minx miny maxx maxy]
%
% Output
%   coor - stickman coordinates coor(:,n) = [x1 y1 x2 y2]' in the following order: 
%          torso, left upper arm, right upper arm, left lower arm, right lower arm, head



  normalizded_det_size = [90 100];
  x1 = rand(1,6)*normalizded_det_size(2);
  x2 = rand(1,6)*normalizded_det_size(2);
  y1 = rand(1,6)*normalizded_det_size(1);
  y2 = rand(1,6)*normalizded_det_size(1);
 coor = round([x1; y1; x2; y2]);
  % move and scale stickman accoridng to det_bb
  scale = normalizded_det_size(1)/(det_bb(4)-det_bb(2)+1);
  coor = coor/scale+repmat([det_bb(1) det_bb(2) det_bb(1) det_bb(2)]',1,size(coor,2));
end