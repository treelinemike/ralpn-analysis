% Compare phantom anatomy across multiple scans
% Display/animate overlay in 3D
% Also compute centroid displacements
% Reference rigid registration via the SVD: https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
%
% Author: M. Kokko
% Updated: 18-Aug-2021

% restart
close all; clear; clc;

% options
base_path = 'S:/eit/mak/20210811-phantom-ct/'; % path to where the segmented STL files are located...
doAnimate = 1;   % animate on the screen?
doMakeVideo = 1; % set to 1 to produce a video file; requires imagemagick ('convert') and ffmpeg
videoFileName = 'phantom_align_';
videoFrameRate = 6; % [frames/sec]

% anatomy to load with color values
anatomy_list = {...
    'RK',[0 1 0]; ...
    'LK',[0 1 0]; ...
    'AA',[1 .33 .33];...
    'IVC',[.33 .33 1];...
    'RU',[1 1 .33];...
    'LU',[1 1 .33];...
    'RV',[.33 .33 1];...
    'LV',[.33 .33 1];...
    'RA',[1 .33 .33];...
    'LA',[1 .33 .33];...
    };

% positions of fiducial "donut" markers #3-8 as identified on CT for each scan
fiducial_positions = zeros(3,6,3);
fiducial_positions(:,:,1) = [...
    2.20, 4.13, 99.02, -133.59, -136.13, -135.98;
    -201.47, -201.72, -150.69, -57.59, -137.23, -65.76;
    222.40, 330.39, 380.97, 378.74, 341.19, 192.40;
    ];
fiducial_positions(:,:,2) = [...
    6.24, 8.46, 104.30, -128.91, -131.80, -131.51;
    2.37, 2.40, 53.12, 146.27, 66.40, 138.37;
    220.50, 328.49, 378.67, 378.04, 339.89, 191.10;
    ];
fiducial_positions(:,:,3) = [...
    -98.57, -96.06, 0.78, -231.84, -235.67, -236.27;
    -98.69, -97.06, -46.16, 47.85, -32.48, 37.03;
    224.1, 332.09, 381.29, 379.83, 342.89, 193.50;
    ];

% check distance between fiducals in each scan
dist_check = zeros(size(fiducial_positions,3),size(fiducial_positions,2)-1);
for i = 1:size(fiducial_positions,3)
   for j = 2:size(fiducial_positions,2) 
       dist_check(i,j-1) = norm( fiducial_positions(:,1,i) - fiducial_positions(:,j,i)  );
   end
end
dist_check

% storage for structure displacements
anat_orig = zeros(3,size(fiducial_positions,2));
anat_disp = zeros(size(fiducial_positions,3),size(fiducial_positions,2));

% perform registration
for moving_idx = 2:3
    
    % center fixed points
    q = fiducial_positions(:,:,1);
    q_bar = mean(fiducial_positions(:,:,1),2);
    Y = q - q_bar;
    
    % center moving points
    p = fiducial_positions(:,:,moving_idx);
    p_bar = mean(fiducial_positions(:,:,moving_idx),2);
    X = p - p_bar;
    
    % compute covariance matrix
    W = eye(size(Y,2)); % identity weighting matrix
    S = X*W*Y';
    
    % compute SVD
    [U,~,V] = svd(S);
    
    % compute homogeneous transformation matrix
    W2 = eye(size(U));
    W2(end,end) = det(V*U');
    R = V*W2*U';
    t = q_bar - R*p_bar;
    TF = [R, t; zeros(1,3), 1];

    % transform and compute residual RMSE (ie. fiducial registration error)
    p2 = hTF(p,TF,0);
    mse = mean(vecnorm(q-p2).^2);
    fre = sqrt(mse)

    % create figure and load original anatomy
    figure;
    set(gcf,'Position',[0488 1.242000e+02 6.634000e+02 6.378000e+02]);
    hold on; axis equal;
    for anat_idx = 1 : size(anatomy_list,1)
        this_stl = stlread([base_path 'Phantom_001_' anatomy_list{anat_idx,1} '.stl']);
        fv = reducepatch(this_stl.ConnectivityList,this_stl.Points,0.1);
        patch('Vertices',fv.vertices,'Faces',fv.faces,'EdgeColor',[0.2 0.2 0.2],'FaceColor',[0.7 0.7 0.7]);
        anat_orig(:,anat_idx) = mean(fv.vertices)';
    end
    
    % load and show new anatomy
    for anat_idx = 1 : size(anatomy_list,1)
        this_stl = stlread([base_path 'Phantom_' sprintf('%03d',moving_idx) '_' anatomy_list{anat_idx,1} '.stl']);
        fv = reducepatch(this_stl.ConnectivityList,this_stl.Points,0.1);
        fv.vertices = hTF(fv.vertices',TF,0)';
        patch('Vertices',fv.vertices,'Faces',fv.faces,'EdgeColor',[0.2 0.2 0.2],'FaceColor',anatomy_list{anat_idx,2});
        anat_disp(moving_idx,anat_idx) = norm(mean(fv.vertices)'-anat_orig(:,anat_idx));
    end
    
    % animate if desired
    first_plot_flag = 1;
    saveFrameIdx = 0;
    if( doAnimate || doMakeVideo )
        ang_vals = 0:-5:-360;
    else
        ang_vals = 0;  % 0, -90
    end
    for ang = ang_vals;
        view([ang,24]);
        if(first_plot_flag)
            axis vis3d;
            first_plot_flag = 0;
        end
        set(gca,'XColor','none','YColor','none','ZColor','none');
        drawnow;
        pause(0.1);
        
        % save frames for video if requested
        if(doMakeVideo)
            thisImgFile = sprintf('frame%03d.png',saveFrameIdx);
            saveFrameIdx = saveFrameIdx + 1;
            saveas(gcf,thisImgFile);
%             system(['convert -trim ' thisImgFile ' ' thisImgFile]);  % REQUIRES convert FROM IMAGEMAGICK!
        end
    end
    
    % generate movie with ffmpeg
    if(doMakeVideo)
        system(['ffmpeg -y -r ' num2str(videoFrameRate) ' -start_number 1 -i frame%03d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 ' videoFileName sprintf('%03d',moving_idx) '.mp4']);
        system('rm frame*.png');
    end
end