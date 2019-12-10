function [params, dirs] = spatialCorrelation(params,dirs)
    %% Computes spatial velocity correlation function and bins it by distance.
    
    % Check to see if start and end frames have been set. If not,
    % initialize them to first and last frames of the movie.
    
    if(~isfield(params,'sframe'))
        params.sframe=1;
    end
    
    if(~isfield(params,'eframe'))
        params.eframe=params.nFrames-1;
    end
    
    % Sanity check to make sure no crazy start and end frames were
    % requested.
    params.sframe = max(params.sframe, 1);
    params.eframe = min(params.eframe, params.nFrames-1);
    
    % Setup file names to write output
    masterpivfile = [dirs.pivDir 'masterVels.mat'];
    spCorrMatFile = [dirs.spCorr filesep 'Spatial_Correlation.mat'];
    
    if(~exist(spCorrMatFile, 'file'))
        if(exist(masterpivfile,'file'))
            load(masterpivfile);
        else
            error(['Something went wrong! Cannot find masterVels.mat in the PIV folder. - ' sprintf('\n%s\n',masterpivfile) '- I am quitting!!']);
        end

        [velX,velY]=size(masterVelu(:,:,1));
        %velX=3;velY=3; % Temp matrix for troubleshooting
        distmat=zeros(velX,velY);

        % generate a distance matrix such that each cell contains euclidean
        % distance to the bottom right corner of the matrix. the size of the
        % distance matrix is the same as the velocity field. this matrix is later
        % expanded for easier masking (see below)


        for i = velX:-1:1
            for j = velY:-1:1
                dist = sqrt((velX-i)^2+(velY-j)^2);
                distmat(i,j)=dist;
            end
        end
        clear dist;
        all_dists = unique(distmat);
        all_dist_mats = NaN.*zeros(velX,velY,velX,velY);
        all_nums = NaN.*zeros(size(all_dists,1),params.eframe-params.sframe+1);

        % Reflect distmat vertically and horizontally such that the center element
        % is 0 and radial distances are approximated by pixels around the center
        % as an example, a 3x3 matrix is converted as follows
        %
        %                    NaN       NaN          NaN       NaN       NaN   
        %   0 0 0            NaN       NaN          NaN       NaN       NaN
        %   0 0 0    ===>    NaN       NaN            0    1.0000    2.0000
        %   0 0 0            NaN       NaN       1.0000    1.4142    2.2361
        %                    NaN       NaN       2.0000    2.2361    2.8284

        %{
        tmp_mat=fliplr(distmat);
        distmat=horzcat(distmat,tmp_mat(:,2:end));
        tmp_mat=flipud(distmat);
        distmat=vertcat(distmat,tmp_mat(2:end,:));
        clear tmp_mat;
        %}

        tmp_mat=NaN.*ones(size(distmat,1),size(distmat,2));

        distmat=horzcat(tmp_mat(:,2:end),rot90(distmat,2));
        tmp_mat=horzcat(tmp_mat,tmp_mat(:,2:end));
        distmat=vertcat(tmp_mat(2:end,:),distmat);

        % i and j are indices of the center of the matrix. there is only one zero
        % in the entire matrix.
        [i,j]=ind2sub(size(distmat),find(distmat==0));

        % now loop backwards from the center to the top left of the expanded
        % distance matrix. This way, outcome of each loop is the value of the
        % interrogation window in the forward direction

        a=1;
        for p=i:-1:1
            b=1;
            for q=j:-1:1           
                all_dist_mats(:,:,a,b) = distmat(p:(p+velX-1),q:(q+velY-1));
                b=b+1;
            end
            a=a+1;
        end
        clear a b

        %speed=sqrt(masterVelu.^2+masterVelv.^2);
        % Migration proceeds in the x-direction. So, velocity correlation
        % must be computed in the y-direction. i.e., the lateral component
        % of velocity, v_i,j
        
        speed=masterVelv;
        spCorr=zeros(numel(all_dists),params.nFrames-1);
        
        for currDist=1:numel(all_dists)

            d1=find(all_dist_mats==all_dists(currDist));

            % j1, j1 represent indices of the cell that is currDist units away
            % from the cell represented by j3, j4
            % In other words, j3,j4 represents the cell with 0 in it and j1,j2
            % has currDist in it.

            [j1, j2, j3, j4]=ind2sub(size(all_dist_mats),d1);

            distInds=sub2ind(size(speed),j1,j2);
            zeroInds=sub2ind(size(speed),j3,j4);

            for frame=1:params.nFrames-1
                currFrame=speed(:,:,frame)-meanv(frame);

                distVels=currFrame(distInds);
                zeroVels=currFrame(zeroInds);


                % Some cells could be NaN depending on PIV and/or masking
                % We need to remove the cells that have NaNs in them, from both
                % vectors to equalize number of elements in summation

                % Find NaNs in both vectors together
                allnans=double(~isnan(distVels.*zeroVels));
                allnans(allnans==0)=NaN;

                % Make distVels(NaNs) || zeroVels(NaNs) = NaN
                distVels=distVels.*allnans;
                zeroVels=zeroVels.*allnans;

                % Remove indices that are NaN from both vectors
                distVels(isnan(distVels))=[];
                zeroVels(isnan(zeroVels))=[];

                % Numerator in the correlation coefficient is
                %
                %       d
                %    \-----                       
                %     \     /                    \
                %      |   | v*(i) x v*(i+del_i) |
                %     /     \                    /
                %    /-----
                %      i=1
                % The summation over time come in later. Here we are processing
                % only the current frame

                num=nansum(distVels.*zeroVels); %may not need nansum here

                % Denominator in the correlation coefficient is
                %
                %                                             1/2 
                %  -----                                 -----
                %  |     d                 d                 |
                %  |  \-----            \-----               |
                %  |   \                 \                   |
                %  |    |   v*(i)^2  x    |    v*(i+del_i)^2 |
                %  |   /                 /                   | 
                %  |  /-----            /-----               |
                %  |    i=1               i=1                |
                %  -----                                 -----
                %

                den=sqrt(nansum(distVels.^2).*nansum(zeroVels.^2));

                % Correlation coefficient is then the ratio of numerator to
                % denominator

                spCorr(currDist,frame)= num ./ den;

            end

        end
        save([dirs.spCorr filesep 'Spatial_Correlation.mat'], 'spCorr', 'all_dists', '-mat');
        save([dirs.spCorr filesep 'Spatial_Correlation.txt'], 'spCorr', '-ascii', '-tabs');
    else
        % If correlation file is already generated, load it
        load(spCorrMatFile);
    end
    
    % Spatial velocity correlation coefficient is computed for the whole
    % length of the movie here. Alternatively, parts of the movie could be
    % analyzed.
    
    %avg_cc = mean(spCorr(:,48:end),2);
    avg_cc = mean(spCorr(:,:),2);
    %Bad programming practice. 16 is hard coded as smallest window size.
    real_dist = all_dists.*(16*params.pixelSize);
    
    bin_size = 30;
    d_bins = 0:bin_size:max(real_dist);
    
    binned_cc = zeros(numel(d_bins),1);
    for db=1:numel(d_bins)
        
        if(db==1)
            dinds = find(real_dist==0);
        else
            dinds = find((real_dist > d_bins(db-1)) & (real_dist <= d_bins(db)));
        end
        
        
        binned_cc(db) = sum(avg_cc(dinds))/numel(dinds);
        
        
    end
    
    binned_spCorr = [d_bins' binned_cc];
    save([dirs.spCorr filesep 'Spatial_Correlation.mat'], 'binned_spCorr', 'bin_size', 'spCorr', 'all_dists', '-mat');
    save([dirs.spCorr filesep 'Binned_Spatial_Correlation.txt'], 'binned_spCorr', '-ascii', '-tabs');

    %win_width=10;
    %sliding_window = conv(ones(1,size(avg_cc,1)), ones(1,win_width), 'same');
    %sliding_window = sliding_window';
    %mean_win=conv(avg_cc,ones(1,win_width),'same');%/win_width);
    %mean_win=mean_win./sliding_window;

    %plot(real_dist,avg_cc,'b-')
    %plot(d_bins,binned_cc)
    
    %avg_cc=mean(cc(:,1:48)');
    %plot(linspace(1,numel(all_dists),numel(all_dists)),cc(:,1));

    %plot(linspace(1,numel(avg_cc),numel(avg_cc')),avg_cc,'ob-')
    %hold on
    %avg_cc=mean(cc(:,48:96)');
    %plot(linspace(1,numel(avg_cc),numel(avg_cc')),avg_cc,'.r-')
    %avg_cc=mean(cc(:,96:144)');
    %plot(linspace(1,numel(avg_cc),numel(avg_cc')),avg_cc)
    %avg_cc=mean(cc(:,144:192)');
    %plot(linspace(1,numel(avg_cc),numel(avg_cc')),avg_cc)
    %avg_cc=mean(cc(:,192:240)');
    %plot(linspace(1,numel(avg_cc),numel(avg_cc')),avg_cc)
    %avg_cc=mean(cc(:,240:end)');
    %plot(linspace(1,numel(avg_cc),numel(avg_cc')),avg_cc)
    %hold off
    %ylim([0 1]);
    %xlim([0 100]);
     
end