%%% Cheng Huimin
%%% A0138497M
%%% EE4212 Assignment: Non-Parametric Sampling

close all; clear;

elapsed = tic;

%% STEP 0 Paramters
column = 200;
row = 120;
region_size = 11; % must be odd
half_region_size = floor(region_size/2);

%% STEP 1 Read and Create Image
texture = imread('texture2.jpg');
[tt_row, tt_column, num_channels] = size(texture);
I=zeros(row,column,3,'uint8');
i_offset=50;
j_offset=50;
I (i_offset+1:i_offset+tt_row, j_offset+1:j_offset+tt_column,:) = texture;

disp('texture size:'), disp(size(texture));
disp('generation size:'),disp(size(I));


%% STEP 2 Set Up Neighbouring Index Map
index_map = zeros(row,column,'int32');
index_map(i_offset+1:i_offset+tt_row,i_offset+1:i_offset+tt_column) = -ones(tt_row,tt_column); % initialise textured pixel to neigative one

% the higher the index, the more neibours the pixel has
for i=1:tt_row
    for j=1:tt_column
        index_map = update_neighbour(index_map,region_size,i,j);
    end
end

%% STEP 3 Construct Vectorised Regions

sliding_row = tt_row - region_size + 1;
sliding_col = tt_column - region_size + 1;

regions_mat = zeros(region_size*region_size,sliding_row*sliding_col,num_channels);
for ch = 1:num_channels
    regions_mat(:,:,ch) = im2col(texture(:,:,ch),[region_size region_size],'sliding');
end

% make channel dimension more major than rolling index
% pixel values => color => next rolling window
regions_mat = permute(regions_mat, [1 3 2]);
% stack RGB in one column, for each region
regions_mat2d = reshape(regions_mat,[],size(regions_mat,3),1);

gaussian = fspecial('gaussian', [region_size region_size], region_size / 5);
gaussian_vec = reshape(gaussian, [], 1);
gaussian_vec = repmat(gaussian_vec, num_channels, 1);

%imagesc(gaussian),colorbar;
figure;

count_goal = row*column - tt_row*tt_column;
count = 0;

while true
    %% STEP 4 Finding Pixel with Most Neighbour, To Be Generated
    t4 = tic;
    [M, idx] = max(index_map(:));
    [pi,pj] =  ind2sub(size(index_map),idx);
    if (M == -1)
        break; %% all pixel processed
    end
    td4 = toc(t4);
    
    %% STEP 5 find match for the pixel with maximum nearest neighbour
    t5 = tic;
    [selected_index,ssd_vec] = find_sampling_match(I,regions_mat2d,index_map,region_size,gaussian_vec,pi,pj);
    
    [pi_source, pj_source] = ind2sub([sliding_row, sliding_col],selected_index);
    pi_source =  pi_source + half_region_size;
    pj_source = pj_source + half_region_size;
    
    td5 = toc(t5);
    
    %% STEP 6 update index_map (counting new neighbours)
    t6 = tic;
    index_map(pi,pj) = -1; % mark as generated
    I(pi,pj,:) = texture(pi_source,pj_source,:);
    index_map = update_neighbour(index_map,region_size,pi,pj); % update index_map, with the new generated point
    
    td6 = toc(t6);
    
    count = count+1;
    if(mod(count,500)==0)
        disp( [num2str(count/count_goal*100,'%.1f'),'%']);
        subplot(1,3,1);
        imshow(I);
        hold on
        plot(i_offset+pi_source,j_offset + pj_source, 'rx','LineWidth' , 1.5, 'MarkerSize' , 10)
        subplot(1,3,2);
        imagesc(index_map), colorbar;
        subplot(1,3,3);
        ssd_2d = reshape(ssd_vec,[sliding_row,sliding_col]);
        imagesc(ssd_2d),colorbar;
        
        disp(td4*1000),disp(td5*1000),disp(td6*1000) % conclusion , step 5 takes 95%+ computation time
        pause(0.25)
    end
    
    
end


%%% Result

td_elapsed = toc(elapsed);
disp(['time elapsed: ' , num2str(td_elapsed,'%.2f') , 's'])

subplot(2,1,1);
imshow(I);
subplot(2,1,2);
imagesc(index_map);
colorbar;

imwrite(I,'result.png');


