% We compute the exact integrals for all the genz functions for dimensions
% 1,2,3,5,10,20 and then store them as CSV format
% we use the notation from https://www.sfu.ca/~ssurjano/disc.html
% thus alpha corresponds to the a's, and beta corresponds to the b's
% we fix a = [5,5,5,5,5...]
% and    u = [0.5,....]
% we integrate from 0 to 1

num_genz = 6;

dimensions = [1, 2, 3, 5, 10, 20];
integrals = zeros(length(dimensions), num_genz);

for indx=1:num_genz
    
    for k=1:length(dimensions); dim=dimension(k);
        
        % set the integration parameters     
        alpha = ones(dim, 1)*5;
        beta = ones(dim, 1)*0.5;
        lower_lim = zeros(dim, 1);
        upper_lim = ones(dim, 1);
        
        % compute the integral
        integral_val = genz_integral(indx, dim, 0, 1, alpha, beta);
        
        % store the integral in a 2d-matrix
        integrals(indx, k) = integral_val;
        integral_val
        
    end
end

integrals;
csvwrite("integrals.csv", integrals);