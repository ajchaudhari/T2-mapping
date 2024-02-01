function [mat2d_T2_mapscan, dmat2d_coeffvals, k, dvec, dvec1, dvec2,cj] = curve_fitting21(dmat2d_coeffvals, k, dmat4d_ROI_T2, dcount, cj, dmat3d_cropped_label, mat2d_T2_mapscan, dvec, dvec1, dvec2);

%Input: T2 scan multiplied with ROI mask for a given ROI for 15 TE's
%Output: Monoexponential curve fitted T2 co-effecients for one ROI for one animal
%In method 2, we first curve fit every pixel in the ROI and then averaged the curve fit 
%T2 values for a ROI over all slices
%This is curve fitted accross 15 TE's 
%throughout this code, cx1,cy1,cz1,ck1 are the 4 counters for looping through the 4D image

    i = 1; davg = []; dnum = 0; dgof=[];l=0; x=[]; y1=[]; z=[];
    davg1 = 0; davg2 = 0; davg3 = 0; davg4 = 0; davg5 = 0; davg6 = 0; davg7 = 0; davg8 = 0;
    for cx1=1:size(dmat4d_ROI_T2,1)
        for cy1 = 1:size(dmat4d_ROI_T2, 2)
            for cz1 = 1:size(dmat4d_ROI_T2 , 3)
                if dmat4d_ROI_T2(cx1,cy1,cz1) ~=0
                    dvec_TE1 = transpose(10:10:150);
                    y(:,1) = dmat4d_ROI_T2(cx1,cy1,cz1,:);
                    charmodel = 'a*exp(-x/b)+e';
                    [dfitted_curve,dgof] = fit(dvec_TE1,y,charmodel,'Lower',[20,30,0], 'Upper' , [400, 90, 70] );
                    
                    % Save the coeffiecient values for b,c and d in a matrix
                    %corresponding to their position in the actual scan
                    % This can then be plotted as ROI wise T2 maps
                    dvec_coeffvals = coeffvalues(dfitted_curve);
                    dmat2d_coeffvals(cx1,cy1,cz1,1, cj) = dvec_coeffvals(1,1);
                    dmat2d_coeffvals(cx1,cy1,cz1,2, cj) = dvec_coeffvals(1,2);
                    dmat2d_coeffvals(cx1,cy1,cz1,3, cj) = dvec_coeffvals(1,3);
                    %saving the same co-efficient values as a list 
                    dvec(dnum+1, cj, dcount) = dvec_coeffvals(1,1);
                    dvec1(dnum+1, cj, dcount) = dvec_coeffvals(1,2);
                    dvec2(dnum+1, cj, dcount) = dvec_coeffvals(1,3);
                    
                    x(dnum+1,1) = dvec_coeffvals(1,1);
                    y1(dnum+1,1) = dvec_coeffvals(1,2);
                    z(dnum+1,1) = dvec_coeffvals(1,3);
                    %Collecting the fit values to average them over a ROI                    
                    davg1 = dvec_coeffvals(1,1) + davg1;
                    davg2 = dvec_coeffvals(1,2) + davg2;
                    davg3 = dvec_coeffvals(1,3) + davg3;
                    
                    %Measuring statistical parameters of the fit
                    %these parameters are also averaged over a ROI
                    davg4 = dgof.sse + davg4;
                    davg5 = dgof.rsquare + davg5;
                    davg6 = dgof.dfe + davg6;
                    davg7 = dgof.adjrsquare + davg7;
                    davg8 = dgof.rmse + davg8;
                    dnum = dnum + 1;
                end
            end            
        end
    end
    
    cj1 = cj + (130 * (dcount-1))+1;
    mat2d_T2_mapscan(cj1, k+1) = davg1/dnum;  
    mat2d_T2_mapscan(cj1, k+2) = davg2/dnum;
    mat2d_T2_mapscan(cj1, k+3) = davg3/dnum;
    mat2d_T2_mapscan(cj1, k+4) = davg4/dnum;
    mat2d_T2_mapscan(cj1, k+5) = davg5/dnum;
    mat2d_T2_mapscan(cj1, k+6) = davg6/dnum;
    mat2d_T2_mapscan(cj1, k+7) = davg7/dnum;
    mat2d_T2_mapscan(cj1, k+8) = davg8/dnum;
    
    %Measuring the standard deviation and maximum of all fit values in a ROI
    mat2d_T2_mapscan(cj1, k+9) = std(x(:,1));
    mat2d_T2_mapscan(cj1, k+10) = std(y1(:,1));
    mat2d_T2_mapscan(cj1, k+11) = std(z(:,1));
    mat2d_T2_mapscan(cj1, k+12) = max(x(:,1));
    mat2d_T2_mapscan(cj1, k+13) = max(y1(:,1));
    mat2d_T2_mapscan(cj1, k+14) = max(z(:,1));

    %measuring the % of fit values within 5% of the maximum value in a ROI
    x1 = mat2d_T2_mapscan(cj1, k+12) - (0.05 * mat2d_T2_mapscan(cj1, k+12));
    x2 = mat2d_T2_mapscan(cj1, k+13) - (0.05 * mat2d_T2_mapscan(cj1, k+13));
    x3 = mat2d_T2_mapscan(cj1, k+14) - (0.05 * mat2d_T2_mapscan(cj1, k+14));
    count3 = 0; count1 = 0; count2 = 0;
    l = dvec(:,cj);
    for j = 1:length(l)
        if dvec(j,cj) > x1
            count1 = count1 +1;
        elseif dvec1(j,cj) > x2
            count2 = count2 +1;
        elseif dvec2(j,cj) > x3
            count3 = count3 +1;  
        end
    end

    mat2d_T2_mapscan(cj1, k+15) = 100 * (count1/dnum);
    mat2d_T2_mapscan(cj1, k+16) = 100 * (count2/dnum);
    mat2d_T2_mapscan(cj1, k+17) = 100 * (count3/dnum);
    
    %measuring kurtosis and skew of all the fit values in a ROI
    mat2d_T2_mapscan(cj1, k+18) = kurtosis(x(:,1));
    mat2d_T2_mapscan(cj1, k+19) = kurtosis(y1(:,1));
    mat2d_T2_mapscan(cj1, k+20) = kurtosis(z(:,1));
    mat2d_T2_mapscan(cj1, k+21) = skewness(x(:,1));
    mat2d_T2_mapscan(cj1, k+22) = skewness(y1(:,1));
    mat2d_T2_mapscan(cj1, k+23) = skewness(z(:,1));
    mat2d_T2_mapscan(cj1, k+24) = var(x(:,1));
    mat2d_T2_mapscan(cj1, k+25) = var(y1(:,1));
    mat2d_T2_mapscan(cj1, k+26) = var(z(:,1));
            
end