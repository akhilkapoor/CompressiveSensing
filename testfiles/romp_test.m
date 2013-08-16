addpath('./../algorithms');

vecSize     = 1000; % Size of the vector in the 'original' input space
nbMeasures  = 50; % size of the measurement vector

A = randn( nbMeasures, vecSize );
x = zeros(vecSize, 1);
s = 25;
x(1:s) = rand(s,1);
y = A*x;

n = 1;
disp(['Test ', num2str(n)]); n = n+1;
[xStar1, NbIter] = romp(y,A);
res = norm(x - xStar1);
disp(['Using romp, the result is: ', num2str(res), ' in ', num2str(NbIter), ' iterations.']);

disp(['Test ', num2str(n)]); n = n+1;
[xStar2, NbIter] = romp(y,A, {'maxnbiter',10});
res = norm(x - xStar2);
disp(['Using romp, the result is: ', num2str(res), ' in ' , num2str(NbIter), ' iterations.']);

disp(['Test ', num2str(n)]); n = n+1;
[xStar3, NbIter] = romp(y,A, {'maxnbiter',25});
res = norm(x - xStar3);
disp(['Using romp, the result is: ', num2str(res), ' in ' , num2str(NbIter), ' iterations.']);

disp(['Test ', num2str(n)]); n = n+1;
[xStar4, NbIter] = romp(y,A, {'maxnbiter',50});
res = norm(x - xStar4);
disp(['Using romp, the result is: ', num2str(res), ' in ' , num2str(NbIter), ' iterations.']);

disp(['Test ', num2str(n)]); n = n+1;
[xStar5, NbIter] = romp(y,A, {'maxnbiter',100});
res = norm(x - xStar5);
disp(['Using romp, the result is: ', num2str(res), ' in ' , num2str(NbIter), ' iterations.']);