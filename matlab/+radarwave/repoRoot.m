function root = repoRoot()
%REPOROOT Folder that holds matlab/, configs/ and reports/.
here = fileparts(mfilename('fullpath'));   % .../matlab/+radarwave
root = fileparts(fileparts(here));
end
