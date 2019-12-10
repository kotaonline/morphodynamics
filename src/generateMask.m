function [params, dirs] = generateMask(Img,olayImg,fno,params,dirs)

    %% Find scratch in frame 'Img'

    % define output file names
        
    scratchtif  = [ dirs.overlay filesep sprintf('%03d',fno) '.tif' ];
    scrMaskFile = [ dirs.masktif filesep sprintf('%03d',fno) '.tif' ];
    pivmaskfile = [ dirs.maskmat filesep sprintf('%03d',fno) '.mat' ];
        
    % Check to see if all the files above exist in their respective
    % locations. Generate mask tif and overlay files if they are not
    % present. Load mask tif files if mat files are not present and
    % generate mat files based on the tif files as masks.

    if ~(exist(scratchtif, 'file') && exist(scrMaskFile, 'file') && params.overwrite==0)
            
        [r,c] = size(Img);
        
        Img_MF=medfilt2(Img, [5 5]);
        Img_EF=entropyfilt(Img_MF,true(9));
        Img_EF=mat2gray(Img_EF);
        Img_EF=im2bw(Img_EF,graythresh(Img_EF));

        erode_size=5000;
        disk_el = strel('disk',4);
        Egw = imdilate(Img_EF,disk_el);

        Egw1 = bwareaopen(Egw,erode_size);
        Egw1 = ~bwareaopen(~Egw1,erode_size);
        Egw1 = imerode(Egw1,disk_el);
        Egw1 = bwareaopen(Egw1,erode_size);
        Egw1 = ~bwareaopen(~Egw1,erode_size);
        if(strcmp(params.expType,'scatter'))
            Egw1 = ~Egw1;
        end
        Egw1 = ~bwfill(~Egw1,'holes');
        
        % Debugging
        %imshowpair(Img,Egw1,'blend');return;
        
        scrfig=figure('visible','off');
        scrgcf = gcf;
        scrgcf.PaperPosition = [ 0 0 c/100 r/100 ];
            
        
        imagesc(olayImg);
        axis image off;
        colormap('gray');
        ax = gca;
        ax.Position = [ 0 0 1 1 ];
        
        mskid = 1;
        
        if(strcmpi(params.expType,'scatter'))
            [B,L,n,A] = bwboundaries(Egw1);
            imwrite( ~Egw1, scrMaskFile,  'tif');
            maske(mskid).msk = ~Egw1;
        else
            [B,L,n,A] = bwboundaries(~Egw1);
            imwrite( Egw1, scrMaskFile, 'tif');
            maske(mskid).msk = Egw1;
        end
        
        maske(mskid).idx = [];
        maske(mskid).idy = [];
        hold on
        
        if(n>1)
            for k = 1:length(B)
                boundary = B{k};
                if(numel(boundary) > params.bound_cutoff)
                    maske(mskid).idx = vertcat(maske(mskid).idx,boundary(:,2));
                    maske(mskid).idy = vertcat(maske(mskid).idy,boundary(:,1));

                    plot(boundary(:,2), boundary(:,1), 'r', 'linewidth', 2);
                end
            end
            
        elseif(n==1)

            [maske(mskid).idx, maske(mskid).idy] = ind2sub(size(maske(mskid).msk),find(maske(mskid).msk==1));
            
            if(strcmp(params.expType,'scratch'))
 
                XY=zeros(size(Img,1),size(Img,2));
                XY(sub2ind(size(XY),B{1}(:,1),B{1}(:,2)))=1;
                Mas=zeros(size(Img,1),size(Img,2));
                Mas(1,:)=1;
                Mas(:,1)=1;
                Mas(end,:)=1;
                Mas(:,end)=1;

                [a b]=find(XY(1,:)==1,1);
                XY(a,b)=2;
                [~, b] = find(XY(end,:)==1,1);
                XY(end,b)=2;
                %M=vertcat(ones(1,9),horzcat(ones(7,1),zeros(7),ones(7,1)),ones(1,9));
                Y=XY-Mas;
                XY=XY&Y;
                [px py]=find(XY');
                B=bwboundaries(XY'); 
                            plot(B{1}(:,1),B{1}(:,2),'r','linewidth',2);
            end
            
            if(strcmp(params.expType,'scatter'))
                B=bwboundaries(~Egw1);
                            plot(B{1}(:,2),B{1}(:,1),'r','linewidth',2);
            end
        end
        hold off

        print('-dtiff',scratchtif,'-r100');
        close(scrfig);

        if(n>=1)
            D=dir(params.woco);
            if(size(D,1)==1)
                for mask_cor=1:length(maske)
                    mapp=load(params.woco);
                    [maske(mask_cor).idxw,maske(mask_cor).idyw]=pixel2world(double(maske(mask_cor).idx),double(maske(mask_cor).idy),...
                        double(maske(mask_cor).idx),double(maske(mask_cor).idy),mapp.comap(:,1),mapp.comap(:,2));
                end
            else
                disp('No such world coordinate file present');
            end

            save(pivmaskfile, 'maske');
        end
    else
        if(strcmpi(params.expType,'scatter'))
            maskim = imread(scrMaskFile);
            maske(1).msk=(maskim==255);

            [maske(1).idx, maske(1).idy] = ind2sub(size(maske(1).msk),find(maske(1).msk==1));

            D=dir(params.woco);
            if(size(D,1)==1)
                mapp=load(params.woco);
                [maske(1).idxw,maske(1).idyw]=pixel2world(double(maske(1).idx),double(maske(1).idy),...
                            double(maske(1).idx),double(maske(1).idy),mapp.comap(:,1),mapp.comap(:,2));

            else
                disp('No such world coordinate file present');
            end
        

            save(pivmaskfile, 'maske');
        end
    end
end
