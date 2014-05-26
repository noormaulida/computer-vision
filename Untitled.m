    img = imread('testing/BGN146.jpg');

    datarsz = imresize(img, [512 512]);
    sizes = size(datarsz);

    datadbl = im2double (datarsz);
    datagry = rgb2gray (datadbl);
    
    datahist = histeq(datagry);
    
    b = 255;
    trans_thresh = mrdivide(0:1:b,b);                                   

    for j=1:(b+1)
        value_trans1 = trans_thresh(1,j);
        result4_trans1(j,1) = bweuler(im2bw(datahist,value_trans1),4);
        %figure, imshow(im2bw(datahist,value_trans1));
        fprintf('j=%d value_trans=%d result_trans=%d\n', j, value_trans1, result4_trans1(j,1));
    end

    max_result4_trans1 = max(result4_trans1)
    min_result4_trans1 = min(result4_trans1)
    m = 0;
    n = 0;
    max_sum_result4_trans1 = 0;
    min_sum_result4_trans1 = 0;

    for j = 1:(b+1)
        if result4_trans1(j,1) == max_result4_trans1
            m = m+1;
            max_sum_result4_trans1 = max_sum_result4_trans1 + j;
        elseif result4_trans1(j,1) == min_result4_trans1
            n = n+1;
            min_sum_result4_trans1 = min_sum_result4_trans1 + j;    
        else
        end
    end
    max_sum_result4_trans1
    min_sum_result4_trans1
    m
    n
    threshold_I1 = ((max_sum_result4_trans1)+(min_sum_result4_trans1))/(m+n)
    dataseg =  im2bw(datahist,threshold_I1/(b+1));
    