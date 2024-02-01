%Code to generate T2 maps for a 140x100x9X15 MRI scan
%This is done using method 2 where  we first apply the ROI mask on the T2 scan 
%and then generate T2 maps only for the ROI's
%each pixel in the ROI is fit with a monoexponential decay for 15 time points. 
%Here, we DO NOT average T2 values for a ROI thereby producing a T2 map over 15 echo times for 
%every pixel in each ROI
%So we end up with ROI wise T2 maps 
%throughout this code, cx1,cy1,cz1,ck1 are the 4 counters for looping through the 4D image
% cj is used as a counter for animals

%Hungarian coding convention:
%c - count
%str - string
%d - double
%ui8 - uint8
%sing - single
%mat2d - 2 dimensional matrix
%mat3d - 3 dimensional matrix
%struct - structure
%char - character


%% 
%IMPORTING DATA
strNamesmaps = importfile('names_maps');
stranimal_ids = importfile('animal_ids');
strNameslabels = importfile('names_labels');
cellanimal_ids(:,1) = convertStringsToChars(stranimal_ids(:,1)); 
%initializing variables
%mat2d_T2_mapscan holds the final results, the dvec's hold the curve fitted
%t2 value for every pixel in the ROI dvec(a), dvec1(b) and dvec2(f)
%dmat2d_coeffvals hold the roi wise T2 maps that are plotted
mat2d_T2_mapscan = []; dvec = []; dvec1 = []; dvec2 = [];
dmat2d_coeffvals = zeros(60,80,9,5);

for cj = 120 : length(strNamesmaps) %This goes through the 143 animal scans

    
    charFilename = convertStringsToChars(strNamesmaps(cj,1)); 
    mat4d_t2 = niftiread(charFilename);

    %We take the segmentation mask and resize it from 280x280x59 to 140x100x59
    %the segmentation mask is loaded in the inverse direction along the z axis,so we flip
    %we downsample the segmentation mask using nearest neighbour interpolation
    charFilename1 = convertStringsToChars(strNameslabels(cj,1)); 
    mat3d_t2_label1 = imresize(flip(niftiread(charFilename1),3),1/2 , 'nearest');

    %We pick 9 precalculated slices from the segmentation mask
    ci1=1;
    for ci2 = 18:3:42       
        ui8mat3d_mask(:,:,ci1) = mat3d_t2_label1(:,:,ci2);
        ci1=ci1+1;
    end

%%
    %CROPPING AND VISUALIZATION
    %Function to crop scan and label based on center of mass algorithm
    [dmat4d_cropped,dmat3d_cropped_label] = cropping_images(mat4d_t2,ui8mat3d_mask); 
%         figure;
%     hold on;
%     for ci1 = 1:9
%         subplot(3,3,ci1)
%         imagesc(dmat3d_cropped_label(:,:,ci1))   
%                    caxis manual
%            caxis([1 10]);
%            %caxis([30 70]);
%            colorbar;
%     end
%     hold off;
    
%     %visualize if we got a good fit of the mask with the scan
%     %set ROI's to 0 and background to 1
    for cx1 = 1 : size(dmat3d_cropped_label,1)
            for cy1 = 1 : size(dmat3d_cropped_label,2)
                for cz1 = 1 : size(dmat3d_cropped_label,3)
                    if dmat3d_cropped_label(cx1,cy1,cz1) == 0
                        dmat3d_cropped_label2(cx1,cy1,cz1) = 1;
                    else
                        dmat3d_cropped_label2(cx1,cy1,cz1) = 0;
                    end

                end

            end
    end
    %Then multiply the binary mask with the T2 scan for any given echo time
    %Plot all slices to check fit
    mask8 = dmat4d_cropped(:,:,:,1) .* dmat3d_cropped_label2;
    figure;
    hold on;
    for ci1 = 1:9
        subplot(3,3,ci1)
        imagesc(mask8(:,:,ci1))   
           caxis manual
           %caxis([30 60]);
           caxis([30 300]);
           colorbar;
    end
    hold off;
       
%%
%ISOLATING ROI's FROM T2 SCAN
    %now we multiply each ROI with the T2 map 
    %first we set desired ROI to 1 and rest to 0 to create a binary mask
        dmat2d_coeffvalsb = zeros(size(dmat3d_cropped_label));   
        k = 0;  
        for dcount = 1:9 %number of ROI's    
            [dmat3d_cropped_label1] = Individual_ROI_mask(dmat3d_cropped_label, dcount);                                    
%%
        %CURVE FITTING
        %then we multiply the label for a specific ROI with the T2 map
        %we curve fit each pixel in the ROI and then average the fitted values for that specific ROI
        %Repeat for all ROI's
            dmat4d_ROI_T2 = dmat3d_cropped_label1.*dmat4d_cropped;
            [mat2d_T2_mapscan, dmat2d_coeffvals, k, dvec, dvec1, dvec2, cj] = curve_fitting21(dmat2d_coeffvals, k, dmat4d_ROI_T2, dcount, cj, dmat3d_cropped_label, mat2d_T2_mapscan, dvec, dvec1, dvec2);
            
       end 

       davg = 0; dnum = 0;
%        figure;
%        hold on;
%        for ci1 = 1:9
%            subplot(3,3,ci1)
%            imagesc(dmat2d_coeffvals(:,:,ci1,2,cj))  
%            caxis manual
%            caxis([40 60]);
%            colorbar;
%        end
%        hold off; 
        
end  

%% 
%EXPORT DATA TO EXCEL
%Output T2 values for all ROI's for all animals in an excel file
char_filename2 = 'T2values_method2_mono.xlsx';
measures = {' ','a', 'b', 'f', 'sse', 'rsquare', 'dfe', 'adjrsquare', 'rmse', 'sd_a', 'sd_b', 'sd_f', ...
    'max_a', 'max_b', 'max_f', 'peak_a', 'peak_b', 'peak_f'};
cellanimal_ids1 = repmat(cellanimal_ids, [9,1]);
writecell(measures, char_filename2, 'Sheet', 1 , 'Range' , 'A1');
writecell(cellanimal_ids1, char_filename2, 'Sheet' , 1, 'Range' , 'A2');
writematrix(mat2d_T2_mapscan, char_filename2, 'Sheet', 1 , 'Range' , 'B1');
writematrix(dvec, char_filename2, 'Sheet', 2 , 'Range' , 'B1');
writematrix(dvec1, char_filename2, 'Sheet', 3 , 'Range' , 'B1');
writematrix(dvec2, char_filename2, 'Sheet', 4 , 'Range' , 'B1');