% 
% customPos By Morten Juelsgaard
% Last edit 16/02/2012    
%
% Custom positioning of figures on the screen
%
% Syntax:
%   customPos(nrCols, nrRows, plotNr)
%
% nrCols: Positive scalar integer 
%    The number of columns the screen should be divided into
% nrRows: Positive scalar integer 
%    The number of rows the screen should be divided into
% plotNr: Positive scalar integer less than or equal nrCols*nrRows 
%    The number of the plot
%
% Demonstration:
% The function call
%   customPos('example')
%
% Runs the following example code:
% 
% close all;
% nrCols = 2;
% nrRows = 3;
%
% x = -10:10;
% for i = 1:nrCols*nrRows
%   figure;
%   customPos(nrCols,nrRows,gcf);
%   plot(x,x.^i);
%   title(sprintf('Sample figure in position nr. %d',gcf));
% end

function customPos(varargin)

    if(length(varargin) == 1 && strcmp(varargin{1},'example'))
        close all;
        nrCols = 2;
        nrRows = 3;

        x = -10:10;
        for i = 1:nrCols*nrRows
            figure;
            customPos(nrCols,nrRows,gcf);
            plot(x,x.^i);
            title(sprintf('Sample figure in position nr. %d',gcf));
        end
    elseif(length(varargin) == 3 && ... %nr of inputs
            varargin{1} > 0 && varargin{2} > 0 && varargin{3} > 0 && ... % positive inputs
            varargin{3} <= varargin{1}*varargin{2} && ... % nr. of plot bounds
            rem(varargin{1},1)==0 && rem(varargin{2},1)==0 && rem(varargin{3},1)==0 &&... % integer check
            max(size(varargin{1}))==1 && max(size(varargin{2}))==1 && max(size(varargin{3})==1) &&... % scalar check
            min(size(varargin{1}))==1 && min(size(varargin{2}))==1 && min(size(varargin{3})==1)) % scalar check
        
        % Retreive inputs
        col = varargin{1};
        row = varargin{2};
        nr = varargin{3};
        
        % Screen size
        screen = get(0,'ScreenSize');
        
        % Standard sizes
        heightOverhead = 0;%30; % Pixel - Screen Estimate
        stdHeight = 420; % Pixel
        stdWidth = 600; % Pixel
        xMin = 0;%1; % Pixel
        yMin = 60;%25; % Pixel

        % Find width and height of figures when printed
        curHeight = max(stdHeight,(screen(4)-yMin -heightOverhead)/row);
        curWidth = min(stdWidth,(screen(3)-2*xMin)/col);
        
        % Calculate positions
        posMat = zeros(col*row,4);
        n = 1;
        for i = 1:col
            for j= 1:row
                % startX, startY, width, height
                posMat(n,:) = [xMin+curWidth*(i-1) yMin+(curHeight)*(j-1)...
                                curWidth, curHeight];
                n = n+1;
            end
        end
        pos = posMat(nr,:);
        
        % Set figure position
        % Use 'OuterPosition' to account for toolbar, menu etc.
        set(gcf,'OuterPosition',pos);
    else
        disp(sprintf('\nInvalid parameter input, try: \n\n  help customPos \n\nfor assistance'))
    end
end