function [mat2d_T2_mapscan,cj] = isolation_rois(dmat4d_cropped,dmat3d_cropped_label,cj,mat2d_T2_mapscan)

%For method 3
%In this function, we convert the ROI mask into a binary mask for a given ROI
%Then the T2 scan is multiplied with this binary mask to isolate the T2 values for that given ROI
%once all the T2 values for a given ROI are isolated, they are averaged and
%sent to the curve fitting function
%Input: T2 scan, labels
%output: averaged T2 values over a ROI for 15 TE's (9x15 matrix)

    for ck1 = 1 : size(dmat4d_cropped,4)
        davg = 0; dnum = 0;
        for dcount = 1:9 %number of ROI's
            for cx1 = 1 : size(dmat3d_cropped_label,1)
                for cy1 = 1 : size(dmat3d_cropped_label,2)
                    for cz1 = 1 : size(dmat3d_cropped_label,3)
                        if dmat3d_cropped_label(cx1,cy1,cz1) == dcount
                            dmat3d_cropped_label1(cx1,cy1,cz1) = 1;
                        else
                            dmat3d_cropped_label1(cx1,cy1,cz1) = 0;
                        end
                    end
                end
            end 

        %then we multiply the label for a specific ROI with the T2 map
        %we take an average for the T2 value in that specific ROI
        %Repeat for all ROI's
            dmat4d_ROI_T2 = dmat3d_cropped_label1.*dmat4d_cropped;
            for cx1=1:size(dmat4d_ROI_T2,1)
                for cy1 = 1:size(dmat4d_ROI_T2, 2)
                    for cz1 = 1:size(dmat4d_ROI_T2 , 3)
                        if dmat4d_ROI_T2(cx1,cy1,cz1) ~= 0
                            davg = davg + dmat4d_ROI_T2(cx1,cy1,cz1,ck1);
                            dnum = dnum + 1;
                        else
                            % do nothing
                        end
                    end            
                end
            end
            mat2d_T2_mapscan(dcount, ck1, cj) = davg/dnum;                
            dcount = dcount + 1;
        end 
    end
end

