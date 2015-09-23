function s = ENT_TestHandle(H)

s = ~isempty(H);                            % Check that input isn;t empty                               
if s == 1                                   % If not empty...
    if exist('ishandle.m','file')           % Pre-MATLAB 2014b
        s = ishandle(H);                    % Check if input is a graphics handle
    elseif exist('isgraphics.m','file')     % MATLAB 2014b onward
        s = isgraphics(H);                	% Check if input is a graphics handle
    end
end