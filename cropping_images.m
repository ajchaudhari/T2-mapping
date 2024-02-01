function [dmat4d_cropped,dmat3d_cropped_label] = cropping_images(mat4d_t2,ui8mat3d_mask)

%This function performs a cropping using the center of mass algorithm. Since method 2 uses a pixel 
%wise curve fitting algorithm, cropping the scan reduces the computational load. 
%Cropping is performed for both the T2 scan as well as the label. 
%For the T2 scan it needs to be performed 15 times for all the TE’s whereas for the label it’s performed only once. 
%In both cases, the cropping is done over all the slices.

%First the center of mass is calculated along X and Y for the T2 scan for TE 1, 
%this same center of mass is used to crop the remaining 14 TE’s as well as the label. 
%This process is repeated for every slice.

%Incase, the center of mass exceeds the borders of the scan, 
%then the center of the scan is used as the center of mass and the middle portion of the scan is cropped.
%throughout this code, cx1,cy1,cz1,ck1 are the 4 counters for looping through the 4D image

cx1=1; cy1=1; cz1=1; ck1=1;
    %we crop the image with a center of gravity algorithm
    for cz1 = 1 : size(mat4d_t2,3)
        ui8mat2d_J2(:,:) = ui8mat3d_mask(:,:,cz1)';
    
        for ck1 = 1 : size(mat4d_t2, 4)
            singmat2d_scan(:,:) = mat4d_t2(:,:,cz1,1)'; 
            singmat2d_scan1(:,:) = mat4d_t2(:,:,cz1,ck1)'; 
            dvec_x = 1 : size(singmat2d_scan, 2); % Columns.
            dvec_y = 1 : size(singmat2d_scan, 1); % Rows.
            [dvec_X, dvec_Y] = meshgrid(dvec_x, dvec_y);
            sing_meanA = mean(singmat2d_scan(:));
            sing_centerOfMassX = round( mean(singmat2d_scan(:) .* dvec_X(:)) / sing_meanA);
            sing_centerOfMassY = round( mean(singmat2d_scan(:) .* dvec_Y(:)) / sing_meanA);

            %cropping the segmentation mask as well during the first TE

            if sing_centerOfMassY > 30 && sing_centerOfMassX >40
                if ck1 == 1
                dmat3d_cropped_label(:,:,cz1) = double(ui8mat2d_J2(sing_centerOfMassY-30:sing_centerOfMassY+29,sing_centerOfMassX-40:sing_centerOfMassX+39));
                end
                dmat4d_cropped(:,:,cz1,ck1) = singmat2d_scan1(sing_centerOfMassY-30:sing_centerOfMassY+29,sing_centerOfMassX-40:sing_centerOfMassX+39);



            else
                if ck1 == 1
                dmat3d_cropped_label(:,:,cz1) = double(ui8mat2d_J2(1:60,1:80));
                end
                dmat4d_cropped(:,:,cz1,ck1) = singmat2d_scan1(1:60,1:80); 

            end
            singmat2d_scan1 = [];
        end
    end
end

