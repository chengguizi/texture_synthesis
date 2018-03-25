%%% Cheng Huimin
%%% A0138497M
%%% EE4212 Assignment: Non-Parametric Sampling

%%% This function assume the point to be generated is at [pi,pj]. The match
%%% will be found using Gaussian-weighted SSD, based on region_size box.
%%% The optimal sampling pixel will be returned as [mi,mj]

%%% regions_mat2d consist of all vectorised sliding region (column major)
%%% in its columns. RGB channels are stacked in each column
function [selected_index, ssd_vec] = find_sampling_match(I,regions_mat2d,index_map,region_size,gaussian_vec,pi,pj)
    
    error_threshold = 0.1;
    [row, column, num_channels] = size(I);
    half_region_size = floor(region_size/2);
    top_left_i = pi - half_region_size;
    top_left_j = pj - half_region_size;
    i0 = max(top_left_i,1);
    j0 = max(top_left_j,1);
    
    bottom_right_i = pi + half_region_size;
    bottom_right_j = pj + half_region_size;
    i1 = min(bottom_right_i,row);
    j1 = min(bottom_right_j,column);
    
    %% STEP 1 construct the mask vector & neighbourhood pixel vector
    t1 = tic;
    mask = zeros(region_size,region_size,'logical');
    neighbours = zeros(region_size,region_size,num_channels);
    for i=i0:i1
        for j=j0:j1  
            if (index_map(i,j) == -1) % if the pixel is generated, mark on the corresponding mask position
                mask(i-top_left_i+1,j-top_left_j+1) = 1;
                %neighbours(i-top_left_i+1,j-top_left_j+1,:) = I(i,j,:);
            end
        end
    end
    neighbours (i0-top_left_i+1:i1-top_left_i+1,j0-top_left_j+1:j1-top_left_j+1,:)= I(i0:i1,j0:j1,:);
    td1 = toc(t1);
    %imshow(neighbours);
    %% STEP 1b reshape/repmat
    t1b = tic;
    mask_vec = reshape(mask,[],1); % vectorise mask into (* by 1) column vector
    mask_vec = repmat(mask_vec, num_channels, 1);
    
    neighbours_vec = reshape(neighbours,[],1);
    neighbours_mat = repmat(neighbours_vec,1,size(regions_mat2d,2));
    td1b = toc(t1b);
    %% STEP 2 obtain gaussian adjusted mask_vec
    t2 = tic;
    mask_norm = sum(mask_vec .* gaussian_vec);
    mask_vec = (mask_vec .* gaussian_vec) / mask_norm;
    td2 = toc(t2);
    %% STEP 3 obtain Sum of Sqaured Distance with the guassian mask   
    t3 = tic;
    ssd_vec = mask_vec' * ((regions_mat2d - neighbours_mat).^2);
    td3 = toc(t3);
    %% STEP 4 find randomised 'minimal' ssd distance, with a arbiturary threshold
    t4 = tic;
    ssd_cutoff = min(ssd_vec)*(1+error_threshold);
    ssd_good_index = find(ssd_vec <= ssd_cutoff);
    
    selected_index = ssd_good_index(randi(length(ssd_good_index)));
    %selected_ssd = ssd_vec(selected_index);
    td4 = toc(t4);
    
    %disp(td1*1000),disp(td1b*1000),disp(td2*1000),disp(td3*1000),disp(td4*1000),disp(' ')
end