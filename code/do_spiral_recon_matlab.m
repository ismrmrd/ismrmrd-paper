% 
% For this code to work, you must add <ISMRMRD_SOURCE>/matlab to your path
% You also need to download the test data
%      https://github.com/ismrmrd/ismrmrd/releases/download/v1.2.3-data/ismrmrd_data.zip
% This code requires:
%   - vds.m for calculation of spiral trajectories (download from:
%   http://www-mrsrl.stanford.edu/~brian/vdspiral/)
%   - NUFFT toolbox for reconstruction (download from:
%   http://web.eecs.umich.edu/~fessler/irt/irt/nufft/)

%-----------------------------------------------------------------------------------------------

% This is an example of how to reconstruct images from data
% acquired with a fully sampled spiral k-space trajectories
%
% ISMRMRD trajectories: k-space trajectories stored in the ISMRMRD file are predicted using the
% gradient impulse response function for distortion correction
%
% Nominal trajectories: variable density spiral trajectories calculated using standard software

function do_spiral_recon_matlab() 

filename = 'spiral.h5';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Check that we have the necessary files  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist(filename, 'file')
    disp('data file exists.');
else
    error(['File ' filename ' does not exist.  Please generate it.'])
end

if (exist('vds.m','file')~=2)
    % try again using the shipped one
    addpath(fullfile('extern','vdspiral'))
    if (exist('vds.m','file')~=2)
        error('vds.m is not installed. Please visit http://www-mrsrl.stanford.edu/~brian/vdspiral/ to install');
    end
end
disp('vds.m is installed.'); 

if (exist('nufft_init.m','file')~=2 || exist('nufft_adj.m','file')~=2)
    % try again using the shipped one
    currfolder = cd(fullfile('extern','irt'));
    setup;
    cd(currfolder);
    if (exist('nufft_init.m','file')~=2 || exist('nufft_adj.m','file')~=2)
        error('NUFFT toolbox is not installed. Please visit http://web.eecs.umich.edu/~fessler/irt/irt/nufft/ to install');
    end
end
disp('NUFFT toolbox is installed.'); 
    
%%%%%%%%%%%%%%%
%  Read data  %
%%%%%%%%%%%%%%%

raw_data = h5read(filename,'/dataset/data');

interleaves = max(raw_data.head.idx.kspace_encode_step_1)+1;
samples = raw_data.head.number_of_samples(1);
channels =raw_data.head.active_channels(1);  

matrix_size = [192,192];

%--Calculate nominal spiral trajectories from freely available code--%%
dt = 2.6e-6; 
[k,g] = vds(14414.4,2.4,dt,double(interleaves),30,3.2);

%--Loop through spiral interleaves--% 
data = zeros(samples, interleaves, channels);
trajectory = zeros(samples,interleaves, 2); 
weights = zeros(samples,interleaves); 
trajectory_nominal = zeros(samples,interleaves,2); 
gradients_nominal = zeros(samples,interleaves,2); 
for i=1:interleaves;
        
    %%--Collect Data--%%
    d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
    d = reshape(d, samples, 1, channels);
    data(:,i,:) = d;  
        
    %%--Collect ISMRMRD trajectories and density compensation--%% 
    % trajectory_dimensions = 3; 
    %(1 = kx,   2 = ky,   3 = density compensation weights) 
    trajectory(:,i,1) = raw_data.traj{i}(1:3:end); 
    trajectory(:,i,2) = raw_data.traj{i}(2:3:end); 
    weights(:,i) = raw_data.traj{i}(3:3:end);        
        
    %%--Rotate nominal trajectories for each interleave--%
    rot = 2*pi*((double(i)-1)/double(interleaves));
    trajectory_nominal(:,i,1) = -(real(k)*cos(rot)+imag(k)*sin(rot));
    trajectory_nominal(:,i,2) =-(-real(k)*sin(rot)+imag(k)*cos(rot));
    gradients_nominal(:,i,1) = -(real(g)*cos(rot)+imag(g)*sin(rot));
    gradients_nominal(:,i,2) =-(-real(g)*sin(rot)+imag(g)*cos(rot));
end
data = reshape(data,interleaves*samples,channels);


%%%%%%%%%%%%%%%%%
%  Reconstruct  %
%%%%%%%%%%%%%%%%%

%%--ISMRMRD trajectories--% 
omega = trajectory*(pi/max(max(max(trajectory))));
omega = reshape(omega,interleaves*samples,2);

weights = reshape(weights, samples*interleaves,1);
    
st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
x = nufft_adj(data.*repmat(weights,[1,channels]), st);

img = squeeze(sqrt(sum(x.*conj(x),3)));
   

%%--Nominal trajectories--% 
omega = trajectory_nominal*(pi/max(max(max(trajectory_nominal))));
omega = reshape(omega,interleaves*samples,2);

    %Calculate density compensation for nominal trajectories%
gradients_nominal = reshape(gradients_nominal,interleaves*samples,2); 
grad = complex(gradients_nominal(:,1),gradients_nominal(:,2));
k = complex(omega(:,1),omega(:,2));
weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(k(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
      
st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
x = nufft_adj(data.*repmat(weights,[1,channels]), st);
    
img_nominal = squeeze(sqrt(sum(x.*conj(x),3)));
    

%%%%%%%%%%%%%%%%%%
% Display images %
%%%%%%%%%%%%%%%%%%

%Crop images
img_pl_n = img_nominal(55:154,25:144); 
img_pl = img(55:154,25:144); 

%Plot images
close all;
figure; 
axes('Position',[0.02 0.45 0.35 0.35]);
    imagesc(img_pl_n,[0 2.5]); axis image; axis off; colormap gray; title('Nominal Trajectories','FontSize',16) 
axes('Position',[0.32 0.45 0.35 0.35]);
    imagesc(img_pl,[0 2.5]); axis image; axis off; colormap gray; title('ISMRMRD Trajectories','FontSize',16) 
axes('Position',[0.62 0.45 0.35 0.35]);
    imagesc(img_pl_n-img_pl,[-1 1]); axis image; axis off; colormap gray; title('Difference image','FontSize',16) 
   
%%%%%%%%%%%%%%%%%
% Output images %
%%%%%%%%%%%%%%%%%
reconImages=zeros([size(img_pl_n),3]);
reconImages(:,:,1)=img_pl_n;
reconImages(:,:,2)=img_pl;
reconImages(:,:,3)=img_pl_n-img_pl;
% Append to the ismrmrd file as an NDArray
h5create(filename,'/dataset/matlab',size(reconImages),'Datatype','single')
h5write(filename,'/dataset/matlab',reconImages)

% close up shop
close all;

exit;

end
