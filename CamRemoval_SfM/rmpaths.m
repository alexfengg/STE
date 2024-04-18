function rmpaths
% add the current directory and all of its subdirectories to PATH

rmpath(genpath(['.' filesep]));

end