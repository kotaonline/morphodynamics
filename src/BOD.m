%{
clc
clear all
close all

%masterpivfile='../10A-CPyo-Scratch-0001_shNeg-1-scaled-drifted/PIV/masterPIV.mat';

%load masterpivfile;

shnegMat = '../test/test-L/PIV/masterVels.mat';
shsurMat = '../shSur8-scatter-3/PIV/masterVels.mat';

[absN, recovN] = bodAnalysis(shnegMat);
[absS, recovS] = bodAnalysis(shsurMat);

figure;
plot(absN, recovN, 'ob')
hold on
plot(absS, recovS, 'or')
xlim([0 20]);
ylim([0 1]);
%}
function [xvals, recov_E] = BOD(pivmatfile, varargin)
    
    fileFlag = 0;
    
    switch nargin
        case 4
            if(ischar(varargin{1}))
                fileFlag = 1;
                outFile  = varargin{1};
            else
                error('Invalid input: Usage - bodAnalysis(pivmatfile, <outputFile>, <startFrame>, <endFrame>)');
            end
            
            startFrame = varargin{2};
            endFrame   = varargin{3};
        case 2
            if(ischar(varargin{1}))
                fileFlag = 1;
                outFile = varargin{1};
            else
                error('Invalid input: Usage - bodAnalysis(pivmatfile, <outputFile>, <startFrame>, <endFrame>)');
            end
    end
    
    pvf = pivmatfile;
    load(pvf);
    
    nFrm = size(masterVelu,3);
    
    if(nargin<4)
        startFrame = 1;
        endFrame   = nFrm;
    end
    
    endFrame = min(endFrame,nFrm);
    startFrame=1;
    endFrame = 168;
    uMat = masterVelu(:,:,startFrame:endFrame);
    vMat = masterVelv(:,:,startFrame:endFrame);
    [pivx, pivy, pivt] = size(uMat);

    tcorr_mat = NaN.*zeros(pivt);

    for nowT=1:pivt
        for nexT=nowT:pivt
            uframe1=uMat(:,:,nowT);
            uframe2=uMat(:,:,nexT);
            vframe1=vMat(:,:,nowT);
            vframe2=vMat(:,:,nexT);

            uframe12=uframe1.*uframe2;
            vframe12=vframe1.*vframe2;
            uvframe12 = uframe1.*vframe2;
            vuframe12 = vframe1.*uframe2;

            tcorr_frame = uframe12+vframe12+uvframe12+vuframe12;
            tcorr_coef = nansum(nansum(tcorr_frame));
            tcorr_mat(nowT,nexT)=tcorr_coef;
            tcorr_mat(nexT,nowT)=tcorr_coef;
        end
    end

    kE = sum(diag(tcorr_mat));
    tcorr_mat  = tcorr_mat./kE;

    [eigvec, eigval] = eig(tcorr_mat);
    
    [eigval, eigind] = sort(diag(eigval),'descend');

    %eigvec = flipud(eigvec);
    %eigval = flipud(eigval);
    eigvec = eigvec(:,eigind);
    eigval2 = cumsum(eigval);

    totalE = sum(eigval);

    recov_E = eigval2./totalE;
    xvals = linspace(1,pivt,pivt);

    if(fileFlag)
        fileH = fopen([outFile '.txt'],'w');
        tmpArr = [xvals; recov_E'];
        fprintf(fileH, '%f %f\n', tmpArr);
        fclose(fileH);
        clear tmpArr;
    
        [fdir, fname, fext] = fileparts(outFile);
        bodfig = figure('Visible','off');
    
        ax = axes('Position', [0.1 0.1 0.8 0.8]);
        plot(ax, xvals, recov_E, 'bo')
        
        title('Biorthogonal Decomposition');
        xlabel('Mode');
        ylabel('{\xi} (Fraction of energy recovered)');
        xlim([0 20]);
        ylim([0 1]);
        
        print('-dpsc', [fdir filesep fname '.ps']);
        close(bodfig);

    end
    %{
    [~, ~, evmat] = meshgrid(masterVelu(1,:,1)',masterVelu(:,1,1),eigvec(:,1));
    u = (1/sqrt(eigval(1,1))).*nansum(masterVelu(:,:,startFrame:endFrame).*evmat,3);
    v = (1/sqrt(eigval(1,1))).*nansum(masterVelv(:,:,startFrame:endFrame).*evmat,3);
    x = masterVelu(:,:,1);
    y = masterVelv(:,:,1);
    %}
    %u=masterVelu(:,:,48);
    %v=masterVelv(:,:,48);

    %[curlz,cav] = curl(u,v);
    %figure;
    %hist(reshape(cav,[numel(cav),1]));
    %xlim([-1 1]);
    %{
    figure;
    pcolor(cav), 
    shading interp,
    
    colorbar
    caxis([-1,1])
    hold on
    quiver(u,v);
    hold off
    colormap('jet');
    %}
    %figure;
    %hist(reshape(cav,[numel(cav),1]))
    %quiver(proj_frame_u,proj_frame_v);
    %{
    scorr_mat = zeros(39*2,2*pivy,pivt);
    for t=1:pivt
        scorr_uu_mat = masterVelu(:,:,t).*masterVelu(:,:,t);
        scorr_uv_mat = masterVelu(:,:,t).*masterVelv(:,:,t);
        scorr_vv_mat = masterVelv(:,:,t).*masterVelv(:,:,t);
        temp_mat = horzcat(scorr_uu_mat(1:39,:),scorr_uv_mat(1:39,:));
        temp_mat = vertcat(temp_mat,horzcat(scorr_uv_mat(1:39,:),scorr_vv_mat(1:39,:)));
        scorr_mat(:,:,t) = temp_mat;
    end
    scorr_mat = nansum(scorr_mat,3);
    
    %{
    for t=1:pivt
        uframe = masterVelu(:,:,t);
        vframe = masterVelv(:,:,t);
        
        for x=1:pivx
            for y=1:pivy
                scorr_uu_mat(x,y) = scorr_uu_mat(x,y)+(uframe(x,y)*uframe(x,y));
                scorr_uv_mat(x,y) = scorr_uv_mat(x,y)+(uframe(x,y)*vframe(x,y));
                %scorr_vu_mat(x,y) = scorr_vu_mat(x,y)+(vframe(x,y)*uframe(x,y));
                scorr_vv_mat(x,y) = scorr_vv_mat(x,y)+(vframe(x,y)*vframe(x,y));
            end
        end
    end
    %}
    %figure;
    %plot(xvals,recov_E,'ob')
    %hold on
    
    kE = sum(diag(scorr_mat));
    scorr_mat  = scorr_mat./kE;
    
    [eigvec, eigval] = eig(scorr_mat);
    
    eigval = diag(eigval);
    
    %eigvec = flipud(eigvec);
    eigval = flipud(eigval);

    eigval2 = cumsum(eigval);

    totalE = sum(eigval);

    recov_E = eigval2./totalE;
    
    xvals = linspace(1,78,78);
    
    %plot(xvals,recov_E,'or')
    %hold off
    %figure;
    %plot(xvals,recov_E,'or');
    %plot(xvals,eigvec(:,1)','-b')
    %ylim([-1,1])
    %hold on
    %plot(xvals,eigvec(:,2)','-r')
    %plot(xvals,eigvec(:,3)','-g')
    %plot(xvals,eigvec(:,4)','-k')
    %hold off
    %}
end