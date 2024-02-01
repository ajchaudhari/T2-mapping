function [dmat3d_cropped_label1] = Individual_ROI_mask(dmat3d_cropped_label, dcount)

%For method 2
%This function sets the values for a particular ROI to 1 and the rest of it to 0
%Thereby creating a ROI specific binary mask 
%Input: dmat3d_cropped_label (cropped labels with all ROIs), dcount(desired ROI to be set to 1)
%Output: dmat3d_cropped_label1 (cropped label with only 1 ROI)

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
end

