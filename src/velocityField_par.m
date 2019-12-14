function [params, dirs] = velocityField_par(params,dirs)

    %% Processing each frame to compute velocity field
    
    masterpivfile = [dirs.pivDir filesep 'masterPIV.mat'];
    kymomat       = [ dirs.kymo 'kymograph.mat'];
    mastervelmat  = [ dirs.pivDir 'masterVels.mat'];
    params.min_win_size = 16;   
    % Check to see if PIV has already been done on this image. The only way
    % I check is to see if masterPIV.mat exists in the PIV sub-directory.
    
    if(~exist(masterpivfile,'file'))
        local_ext = params.ext;
        local_frdir = dirs.iframes;
        local_nfrms = params.nFrames;
        local_width = params.width;
        local_height = params.height;
        local_tpf   = params.timePerFrame;
        local_olap  = params.overlap;
        local_woco  = params.woco;

        %local_vmax = params.vmax; %need to update original when parfor is done

        % Compute PIV matrix size. This is done in the PIV code. We need to
        % generate this here for initializing arrays to store PIV data from
        % the parfor loop.
        
        interrMatSize = [floor((local_height-params.min_win_size)/((1-local_olap)*params.min_win_size))+1,...
            floor((local_width-params.min_win_size)/((1-local_olap)*params.min_win_size))+1];
        
        %[params.interrWinX, params.interrWinY] = interrMatSize;

        for i = local_nfrms-1:-1:1
            masterPIV(i).u = NaN.*zeros(interrMatSize);
            masterPIV(i).v = NaN.*zeros(interrMatSize);
            masterPIV(i).x = NaN.*zeros(interrMatSize);
            masterPIV(i).y = NaN.*zeros(interrMatSize);
            masterPIV(i).snr_u = NaN.*zeros(interrMatSize);
            masterPIV(i).snr_v = NaN.*zeros(interrMatSize);
            masterPIV(i).final_u = NaN.*zeros(interrMatSize);
            masterPIV(i).final_v = NaN.*zeros(interrMatSize);
            masterPIV(i).snr = NaN.*zeros(interrMatSize);
            masterPIV(i).pkh = NaN.*zeros(interrMatSize);
        end
        
        parfor i = 1 : local_nfrms-1

            % define input file names -- two successive frames
            ifname1 = [ local_frdir filesep sprintf('%03d',i) local_ext ];
            ifname2 = [ local_frdir filesep sprintf('%03d',i+1) local_ext ];

            Img1 = imread(ifname1);
            Img2 = imread(ifname2);

            % Leaving the following here for the future. Currently only one
            % edge is supported
            
            %if(local_edges==2)
            %    Img1 = Img1(:,1:local_width/2);
            %    Img2 = Img2(:,1:local_width/2);
            %end

            [x,y,u,v,snr,pkh] = matpiv(Img1, Img2, [64 64;64 64;32 32;32 32;16 16;16 16],local_tpf,local_olap,'multin',local_woco);

            [su,sv] = snrfilt(x,y,u,v,snr,1.3);
            [pu,pv] = peakfilt(x,y,su,sv,pkh,0.2);
            [gu,gv] = globfilt(x,y,pu,pv,3);
            
            % Commented out the following because masking is done later
            % when figures are generated and kymographs are computed.
            
            %if(exist(pivmaskfile,'file'))
            %    [mu,mv] = localfilt(x,y,gu,gv,3,'median',3,pivmaskfile);       
            %    [fu,fv] = naninterp2mask(mu,mv,pivmaskfile,x,y);
            %else
                [mu,mv] = localfilt(x,y,gu,gv,3,'median',3);
                [fu,fv] = naninterp2(mu,mv);
            %end

            masterPIV(i).x = x;
            masterPIV(i).y = y;
            masterPIV(i).u = u;
            masterPIV(i).v = v;
            masterPIV(i).snr_u = su;
            masterPIV(i).snr_v = sv;
            masterPIV(i).final_u = fu;
            masterPIV(i).final_v = fv;

        end
        save(masterpivfile,'masterPIV','-mat');
    else
        load(masterpivfile);
    end
      
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

        masterVelu(:,:,i)  = masterPIV(i).final_u.*interp_pivmask;
        masterVelv(:,:,i)  = masterPIV(i).final_v.*interp_pivmask;
        masterMags(:,:,i)  = sqrt((masterVelu(:,:,i).^2) + (masterVelv(:,:,i).^2)).*60;%.*interp_pivmask; % multiplication by 60 converts to um/hr
        masterTheta(:,:,i) = atan2(masterVelv(:,:,i),masterVelu(:,:,i));
        masterCosT(:,:,i)  = cos(masterTheta(:,:,i));%.*interp_pivmask;
        
        meanu(i)           = nanmean(nanmean(masterVelu(:,:,i)));
        meanv(i)           = nanmean(nanmean(masterVelv(:,:,i)));
        
        vkymograph(:,i)    = nanmean(masterMags(:,:,i));
        skymograph(:,i)    = nanmean(masterCosT(:,:,i));
                
    end
    
    % Figuring out top speed. Not velocity in this case.
    params.sMax = nanmax(nanmax(nanmax(masterMags)));
    
    vkymograph(isnan(vkymograph))=0;
    skymograph(isnan(skymograph))=0;
    
    save(kymomat,'vkymograph','skymograph','-mat');
    save(mastervelmat, 'masterVelu','masterVelv','masterMags','masterCosT', 'masterTheta', 'meanu','meanv','-mat');
    
    
end
