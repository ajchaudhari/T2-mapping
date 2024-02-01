%Code to generate T2 maps for a 140x100x9X15 MRI scan
%This is done using method 3 where apply the ROI mask on the T2 scan and average the T2 values for 
%a given ROI and perform only one curve fit for the averaged values for each ROI
%throughout this code, cx1,cy1,cz1,ck1 are the 4 counters for looping through the 4D image
% cj is used as a counter for scans

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
mat2d_T2_mapscan = [];  %holds the results

for cj = 1 : length(strNamesmaps) %This goes through the 143 animal scans
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
    dmat3d_cropped_label2 = dmat3d_cropped_label;

    %checking fit of mask with scan
%     for cx1 = 1 : size(dmat3d_cropped_label,1)
%         for cy1 = 1 : size(dmat3d_cropped_label,2)
%             for cz1 = 1 : size(dmat3d_cropped_label,3)
%                 if dmat3d_cropped_label(cx1,cy1,cz1) == 0
%                      dmat3d_cropped_label2(cx1,cy1,cz1) = 1;
%                  else
%                      dmat3d_cropped_label2(cx1,cy1,cz1) = 0;
%                 end
%             end
%         end
%     end
%     visual = dmat3d_cropped_label2.*dmat4d_cropped(:,:,:,1);
%     figure;
%     hold on;
%     for ci1 = 1:9
%         subplot(3,3,ci1)
%         imagesc(visual(:,:,ci1))
%     end
%     hold off;
        
    
                    
   
%%   
    %ISOLATING ROI's FROM T2 SCAN
    [mat2d_T2_mapscan] = isolation_rois(dmat4d_cropped,dmat3d_cropped_label,cj,mat2d_T2_mapscan);
end 
    
 %%  
%CURVE FITTING 
%This function performs curve fitting using Matlab's inbuilt exponential curve fitting function
cj = 1;
[y, dvec_TE1, cj, dmat2d_T2_mapscan] = curve_fitting31(mat2d_T2_mapscan);

%% 
%EXPORT DATA TO EXCEL
%Output T2 values for all ROI's for all animals in an excel file
char_filename2 = 'T2values_method3_mono.xlsx';
cellanimal_ids1 = repmat(cellanimal_ids, [9,1]);
measures = {' ','a', 'b', 'f', 'sse', 'rsquare', 'dfe', 'adjrsquare', 'rmse'};
writecell(measures, char_filename2, 'Sheet', 1 , 'Range' , 'A2');
writecell(cellanimal_ids1, char_filename2, 'Sheet' , 1, 'Range' , 'A3');
writematrix(dmat2d_T2_mapscan, char_filename2, 'Sheet', 1 , 'Range' , 'B3');
