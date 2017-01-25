function out = DummyBuffyPoseEstimationPipeline(buffydir,episodenr)
% function out = DummyBuffyPoseEstimation(buffydir,episodenr)
% this routine provides a dummy example of pose estimation pipeline
% for each image in the dataset it returns from 0-10 random detections with a dummy pose
% this routine was prepared to demonstrate the way to obtain data structure suitable for our evaluation routine.
% Input:
%   buffydir - relative or absolute path to buffy frames directory
%   episodenr - episode number (value stored in out(i).episode
% Output:
%   out - is a stucture as excpected by the BatchEvalBuffy evaluation routine, it is an array of structs (one entry per frame) with the following fields:
%           .episode - corresponding episode number
%           .frame - corresponding frame number
%           .stickmen - array of structs containing multiple stickmen to be evaluated with fields:
%             .coor - stickman end-points coordinates coor(:,nparts) = [x1 y1 x2 y2]'
%             .det - is the detection bounding box associated with the stickman in [minx miny maxx maxy]
%           

  Files = dir(buffydir);
  invalid = false(length(Files),1);
  for i=1:numel(Files)
    invalid(i) = isempty(regexpi(Files(i).name, '.jpg'));
  end
  Files(invalid) = [];
  N = length(Files);

  out(N) = struct('frame',[],'stickmen',[]);

  for i=1:N
    img = imread(fullfile(buffydir,Files(i).name));
    out(i).frame = str2double(Files(i).name(1:end-4));
    out(i).episode = episodenr;
    out(i).stickmen = DummyDetect(img);
    for j=1:length(out(i).stickmen)
      out(i).stickmen(j).coor = DummyPose(img,out(i).stickmen(j).det);
    end

  end

end




  