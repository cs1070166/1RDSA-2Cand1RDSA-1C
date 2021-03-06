function y = loss_myexample(p, theta, type)
% type 1 = quadratic, 2 = fourth order function from Spall's 2000 paper
% Create upper triangular matrix with entries 1/p

tmpMat = zeros(1,p);
tmpMat(1,1) = 1;
Amat=1/p*toeplitz(tmpMat,ones(1,p));

% based on type, return the appropriate loss - quadratic if 1, 4th order
% otherwise
if type == 1
    b = ones(p, 1);
    y=theta'*Amat*theta + theta'* b;

elseif type == 2 % Spall 2000
     y=theta'*(Amat'*Amat)*theta + 0.1* sum((Amat*theta).^3) + 0.01 * sum((Amat*theta).^4);

elseif type == 3 % Rosenbrock   
    psum = 0;
    for j = 1:p-1;
        psum = psum+0.01*(theta(j,1)^2-theta(j+1,1))^2+(theta(j,1)-1)^2;
    end
    y = psum;  
else % Powell
    psum=0;
    for j = 2:p-2
        psum = psum + (theta(j-1,1)+10*(theta(j,1)))^2 + 5*(theta(j+1,1)-theta(j+2,1))^2 + 0.01*(theta(j,1)-2*(theta(j+1)))^4 + 0.01*(theta(j-1,1)-theta(j+2,1))^4;
    end;
    y = psum;
end
