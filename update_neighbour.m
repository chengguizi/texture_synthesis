%%% Cheng Huimin
%%% A0138497M
%%% EE4212 Assignment: Non-Parametric Sampling

%%% count neighbouring pixels in the regions
function index_map = update_neighbour(index_map,region_size,i,j)
    [row, column, dummy] = size(index_map);
    half_region_size = floor(region_size/2);

    if (index_map(i,j) ~= -1) % not a generated pixel, ignore
        return
    end
    for k_row=-half_region_size:half_region_size
        for k_col=-half_region_size:half_region_size
            if(i+k_row < 1 || i+k_row>row || j+k_col < 1 || j+k_col > column )
                continue
            end
            if (index_map(i+k_row,j+k_col)~=-1)
                index_map(i+k_row,j+k_col) = index_map(i+k_row,j+k_col)+1;
            end
        end
    end

    
end