function [params, dirs] = makeFigs(params, dirs, varargin)
%% Initialize variables based on varargin input

    if(nargin>=3)
        switch nargin
            case 3
                sMax = varargin{1};
                figR = 300;
                figW = params.width;
                figH = params.height;
            case 4
                sMax = varargin{1};
                figR = varargin{2};
                figW = params.width;
                figH = params.height;
            case 5
                sMax = varargin{1};
                figR = 300;
                figW = varargin{2};
                figH = varargin{3};
            case 6
                sMax = varargin{1};
                figR = varargin{2};
                figW = varargin{3};
                figH = varargin{4};
        end
    else
        sMax = params.sMax;
    end

%% The following generates a pseudocolored figure displaying velocity fields
% 

    load([dirs.pivDir 'masterPIV.mat']);
    load([dirs.pivDir 'masterVels.mat']);
    load([dirs.kymo 'kymograph.mat']);
    
    % The following is futuristic. Commenting out for future support.
    %
    %if(params.edges==2)
    %    figW = params.width/200;
    %else
        figW = figW/figR;
    %end
    figH = figH/figR;
    
    pivcbfile = [ dirs.kymo 'PIV Colorbar.tif'];
    opcbfile = [ dirs.kymo 'S Colorbar.tif'];
    vkymotif = [ dirs.kymo 'Velocity_kymograph.tif'];
    skymotif = [ dirs.kymo 'OrderParam_kymograph.tif'];
    %kymomat = [ dirs.kymo 'kymograph.mat'];
    %mastervelmat = [ dirs.pivDir 'masterVels.mat'];
    
    % Local versions of structure variables for parfor
    
    vfImages  = dirs.vfImages;
    SImages   = dirs.SImages;
    iframes   = dirs.iframes;
    overwrite = params.overwrite;
    expType   = params.expType;

    %gridX = masterPIV(1).x/params.pixelSize;
    %gridY = masterPIV(1).y/params.pixelSize;
    [gridX, gridY] = meshgrid(1:size(masterMags(:,:,1),2),1:size(masterMags(:,:,1),1));
    [finer_gridX, finer_gridY] = meshgrid(1:0.1:size(masterMags(:,:,1),2),1:0.1:size(masterMags(:,:,1),1));
    
    %{
    params.pivx = size(masterPIV(1).x,1);
    params.pivy = size(masterPIV(1).y,1);
    
    % Preallocating memory for kymographs and master PIV and Order
    % parameter matrices

    vkymograph = zeros(size(masterPIV(1).u,2),params.nFrames-1);
    skymograph = vkymograph;
    
    meanu       = NaN.*zeros(params.nFrames-1,1);
    meanv       = meanu;
        
    masterVelu = zeros(size(masterPIV(1).u,1),size(masterPIV(1).u,2),size(masterPIV,2));
    masterVelv = masterVelu;
    masterMags = masterVelu;
    masterTheta = masterVelu;
    masterCosT = masterVelu;
    
    [tmpX, tmpY] = meshgrid(1:params.width,1:params.height);

    for i=1:params.nFrames-1
        pivmaskfile       = [ dirs.maskmat filesep sprintf('%03d',i) '.mat'];
        pivmask           = ones(size(tmpX));
        if(exist(pivmaskfile,'file'))
            load(pivmaskfile);
            pivmask       = double(maske.msk);
        end
        
        interp_pivmask = interp2(tmpX.*params.pixelSize, tmpY.*params.pixelSize, pivmask, masterPIV(i).x, masterPIV(i).y);
 
        interp_pivmask(interp_pivmask==0)=NaN;

        masterVelu(:,:,i) = masterPIV(i).final_u.*interp_pivmask;
        masterVelv(:,:,i) = masterPIV(i).final_v.*interp_pivmask;
        masterMags(:,:,i) = sqrt((masterVelu(:,:,i).^2) + (masterVelv(:,:,i).^2)).*60;%.*interp_pivmask; % multiplication by 60 converts to um/hr
        masterTheta(:,:,i) = atan2(masterVelv(:,:,i),masterVelu(:,:,i));
        masterCosT(:,:,i) = cos(masterTheta(:,:,i));%.*interp_pivmask;
        
        meanu(i)          = nanmean(nanmean(masterVelu(:,:,i)));
        meanv(i)          = nanmean(nanmean(masterVelv(:,:,i)));
                
    end
    
    % Figuring out top speed. Not velocity in this case.
    vMax = nanmax(nanmax(nanmax(masterMags)));
    params.vMax = vMax;
    %}
    % One-time PIV colorbar printing
    if(~exist(pivcbfile,'file'))
        cbf=figure('visible','off');
        cbfgcf = gcf;
        cbfgcf.PaperPositionMode='manual';
        cbfgcf.PaperPosition = [0 0 0.8 figH];
        cbfgcf.InvertHardcopy='off';
        cbfgcf.Color = [1 1 1];
        axis off;
        colormap('jet');
        caxis([0 sMax]);
        cbf.Position=[100 100 80 figH];
        cb=colorbar;
        cb.Location='west';
        cb.Position=[ 0.05 0.05 0.5 0.9 ];
        cb.AxisLocation='in';
        print('-dtiff',pivcbfile,['-r' num2str(figR)]);
        close(cbf);
    end
    
    % One-time Order param colorbar printing
    if(~exist(opcbfile,'file'))
        cbf=figure('visible','off');
        cbfgcf = gcf;
        cbfgcf.PaperPositionMode='manual';
        cbfgcf.PaperPosition = [0 0 0.8 figH];
        cbfgcf.InvertHardcopy='off';
        cbfgcf.Color = [1 1 1];
        axis off;
        colormap('jet');
        caxis([-1 1]);
        cbf.Position=[100 100 80 figH];
        cb=colorbar;
        cb.Location='west';
        cb.Position=[ 0.05 0.05 0.5 0.9 ];
        cb.AxisLocation='in';
        print('-dtiff',opcbfile,['-r' num2str(figR)]);
        close(cbf);
    end
                
    parfor i=1:params.nFrames-1
        
        % Xavier Trepat reports median in Kymographs. I use mean
        
        %vkymograph(:,i)=nanmean(masterMags(:,:,i));
        %skymograph(:,i)=nanmean(masterCosT(:,:,i));
 
        % Define output file names
        %pivhmfile= [ vfImages sprintf('%03d',i) '.tif'];
        pivhmfile = [vfImages sprintf('%03d',i) '.eps'];
        %Sfile = [ SImages sprintf('%03d',i) '.tif'];
        Sfile = [ SImages sprintf('%03d',i) '.eps'];
 
        if ~(exist(pivhmfile, 'file') && exist(Sfile, 'file') && overwrite==0)
                                          
            pivfig=figure('visible','off');
            pivgcf = gcf;
            
            pivgcf.PaperPositionMode = 'manual';
            pivgcf.PaperPosition = [0 0 figW figH ];
            pivgcf.Color = 'black';
            pivgcf.InvertHardcopy = 'off';
            
            if (strcmpi(expType,'scratch'))
                %magMat = interp2(gridX,gridY,masterMags(:,:,i),finer_gridX,finer_gridY);
                %%smoothed_mat = smooth2a(magMat,params.smoothX,params.smoothY);
                %smoothed_mat = smooth2a(magMat,5,5);
            	colormap('jet');
            	image(gridX(1,:),gridY(:,1),masterMags(:,:,i),'CDataMapping','scaled');
                %imagesc(smoothed_mat);
            	caxis([0 sMax]);
            	
                hold on
                
            	q = quiver(gridX,gridY,masterVelu(:,:,i),masterVelv(:,:,i));
            	q.LineWidth=2;
            	q.Color='w';
                %q.MaxHeadSize=0.6;
                q.AutoScaleFactor=1.0;
            	hold off
                ax=gca;
            	ax.Position=[0 0 1 1];            
            	ax.XTick = [];
                ax.YTick = [];
                ax.Color = 'none';
                %print('-dtiff', pivhmfile, ['-r' num2str(figR)]);
                print('-depsc',pivhmfile);
            else
            	%image(tX(1,:),tY(:,1),masterMags(:,:,i),'CDataMapping','scaled');
            	origImg = imread([iframes sprintf('%03d',i) '.tif']);
				imagesc(origImg);
                colormap('gray');
				ax=gca;
            	ax.Position=[0 0 1 1];  
                ax.XTick=[];
                ax.YTick=[];
            	hold on
            	q = quiver(gridX,gridY,masterVelu(:,:,i),masterVelv(:,:,i));
            	q.LineWidth=1;
            	q.Color='r';
                q.MaxHeadSize=0.4;
            	hold off
                print('-dtiff', pivhmfile, ['-r' num2str(figR)]);
			end

            

            close(pivfig);

            oparfig=figure('visible','off');
            opargcf = gcf;
            
            opargcf.PaperPositionMode = 'manual';
            opargcf.PaperPosition = [0 0 figW figH ];
            opargcf.Color = 'black';
            opargcf.InvertHardcopy = 'off';
            
            colormap('jet');

            %pcolor(flipud(masterPIV(i).x),flipud(masterPIV(i).y),masterCosT(:,:,i)), axis off equal tight, shading flat%, colorbar 
            image(gridX(1,:),gridY(:,1),masterCosT(:,:,i),'CDataMapping','scaled');
            
            caxis([-1 1]);
            ax = gca;
            ax.Position = [0 0 1 1];
            ax.XTick = [];
            ax.YTick = [];
            ax.Color = 'none';
        
            hold on
            q = quiver(gridX,gridY,masterVelu(:,:,i),masterVelv(:,:,i));
            %q=quiver(flipud(masterPIV(i).x),flipud(masterPIV(i).y),masterPIV(i).final_u,masterPIV(i).final_v);
            q.LineWidth=2;
            q.Color='w';
            %q.MaxHeadSize=0.6;
            q.AutoScaleFactor=1.0;
            hold off

            %print('-dtiff', Sfile, ['-r' num2str(figR)]);
            print('-depsc', Sfile);

            close(oparfig);
        end          
    end
    
    %vkymograph(isnan(vkymograph))=0;
    %skymograph(isnan(skymograph))=0;
    
    %save(kymomat,'vkymograph','skymograph','-mat');
    %save(mastervelmat, 'masterVelu','masterVelv','masterMags','masterCosT', 'masterTheta', 'meanu','meanv','-mat');
    
    vkymofig = figure('visible','off');
    
    kymox=linspace(1,size(vkymograph,2),size(vkymograph,2)).*(params.timePerFrame/60); % Division by 60 converts to hrs
    
    kymoy=linspace(1,size(vkymograph,1),size(vkymograph,1)).*(params.min_win_size*params.pixelSize); 

    [kymogridx,kymogridy]=meshgrid(1:size(vkymograph,2),1:size(vkymograph,1));
    [kymofinex,kymofiney]=meshgrid(1:0.1:size(vkymograph,2),1:0.1:size(vkymograph,1));
    
    smoothed_kymo=interp2(kymogridx,kymogridy,vkymograph,kymofinex,kymofiney,'linear');
    smoothed_cost=interp2(kymogridx,kymogridy,skymograph,kymofinex,kymofiney,'linear');
    %pcolor(kymox, kymoy, vkymograph), axis tight, shading flat, colorbar
    %image(kymox,kymoy,vkymograph,'CDataMapping','scaled');
    
    imagesc(smoothed_kymo);
    set(gca,'YDir','normal');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    
    
    colorbar
    colormap('jet');
    caxis([0 max(max(vkymograph))]);
    saveas(vkymofig,vkymotif,'tiff');
    close(vkymofig);
    
    skymofig = figure('visible','off');
    %pcolor(kymox, kymoy, skymograph), axis tight, shading flat, colorbar
    %image(kymox,kymoy,skymograph,'CDataMapping','scaled');
    
    imagesc(smoothed_cost);
    set(gca,'YDir','normal');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    
    colorbar
    colormap('jet');
    caxis([-1 1]);
    saveas(skymofig, skymotif, 'tiff');
    close(skymofig);
    
end
    
            
