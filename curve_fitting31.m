function [y, dvec_TE1, cj, dmat2d_T2_mapscan] = curve_fitting31(mat2d_T2_mapscan);

%Input: T2 scan multiplied with ROI mask for a given ROI for 15 TE's
%Output: Monoexponential curve fitted T2 co-effecients for one ROI for one animal
%In method 3, we have averaged the T2 values for a ROI over all slices for
%a single TE, then we perform one curve fit for each ROI over 15 TE's
%throughout this code, cx1,cy1,cz1,ck1 are the 4 counters for looping through the 4D image

k=0;
for cj = 1:size(mat2d_T2_mapscan,3)
    for cx1 = 1:size(mat2d_T2_mapscan,1)
                dvec_TE1 = transpose(10:10:150);
                y(:,1) = mat2d_T2_mapscan(cx1,:,cj);
                charmodel = 'a*exp(-x/b)+e';
                [dfitted_curve,dgof] = fit(dvec_TE1,y,charmodel,'Lower',[20,30,0], 'Upper' , [400, 90, 70]);

                % Save the coeffiecient values for b,c and d in a vector
                dvec_coeffvals = coeffvalues(dfitted_curve);
                
                %Storing the results in such a way that ROI values
                %vertically stack over each other for all animals
                cj1 = cj + (130 * (cx1 -1));
                dmat2d_T2_mapscan(cj1, k+1) = dvec_coeffvals(1,1);  
                dmat2d_T2_mapscan(cj1, k+2) = dvec_coeffvals(1,2);
                dmat2d_T2_mapscan(cj1, k+3) = dvec_coeffvals(1,3);
                
                %statistical measures from each fit
                dmat2d_T2_mapscan(cj1, k+4) = dgof.sse;
                dmat2d_T2_mapscan(cj1, k+5) = dgof.rsquare;
                dmat2d_T2_mapscan(cj1, k+6) = dgof.dfe;
                dmat2d_T2_mapscan(cj1, k+7) = dgof.adjrsquare;
                dmat2d_T2_mapscan(cj1, k+8) = dgof.rmse;
    end
end
    
end