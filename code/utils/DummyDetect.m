function out = DummyDetect(img)
% this is a dummy detection routine which produces from 0 up to 10 random not overlaping boxes with fixed aspect-ratio within the image
% Input:
%   img - input image
% Output:
%   out - array of structs (one entry per bouding box) with .det field containing bouding box coordinates formated as [minx miny maxx maxy]
  max_overlaping = 0.5;
  hwboxaspect = 0.9;
  maxboxes = 10;
  minwidth = 70;
  imsize = size(img);
  maxwidth = imsize(2)/2;

  nbox = round(rand(1)*maxboxes);

  if nbox == 0
    out = [];
    return
  end
  out(nbox) = struct('det',[]);
  for i=1:nbox
    found = true;
    while(found)
      boxW = (maxwidth-minwidth)*rand(1)+minwidth;
      boxsize = [boxW*hwboxaspect boxW];
      minyx = rand(1,2).*imsize([1 2]);
      if minyx < (imsize([1 2]) - boxsize)
        found = false;
      end
    end
    out(i).det = round([minyx(2) minyx(1) minyx(2)+boxsize(2) minyx(1)+boxsize(1)]);
  end
  
  % remove overlaping random boxes, to not break the evaluation routine
  overlap = zeros(length(out));
  for i=1:length(out)
    for j=(i+1):length(out)
      overlap(i,j) = all_pairs_bb_iou(out(i).det', out(j).det');
    end
  end
  [i1 i2]= ind2sub([length(out) length(out)], find(overlap > max_overlaping));
  for i=1:length(i1)
    if ~isempty(out(i1(i)).det)
      out(i2(i)).det = [];
    end
  end
  
  out(arrayfun(@(x)isempty(x.det), out)) = []; % keep only not empty entries
end
