vecSize     = 1000; % Size of the vector in the 'original' input space
nbMeasures  = 200; % size of the measurement vector

A = randn( nbMeasures, vecSize )/sqrt(nbMeasures);
x = zeros(vecSize, 1);
s = 25;
x(1:s) = rand(s,1);
y = A*x;

n = 1;
disp(['Test ', num2str(n)]); n = n+1;
[xStar1,~,~,~,~,~] = htp_back(y,A,s);
res = norm(x - xStar1);
disp(['Using original htp, the result is: ', num2str(res)]);

disp(['Test ', num2str(n)]); n = n+1;
[xStar2,~,~,~,~,~] = htp(y,A,s,{});
res = norm(x - xStar2);
disp(['Using modified htp, the result is: ', num2str(res)]);

disp(['Test ', num2str(n)]); n = n+1;
[xStar3,~,~,~,~,~] = htp(y,A,s,{'MaxNbIter', 500});
res = norm(x - xStar3);
disp(['Using modified htp, the result is: ', num2str(res)]);

disp(['Test ', num2str(n)]); n = n+1;
[xStar4,~,~,~,~,~] = htp(y,A,s,{'maxbIter', 500});
res = norm(x - xStar4);
disp(['Using modified htp, the result is: ', num2str(res)]);

disp(['Test ', num2str(n)]); n = n+1;
[xStar5,~,~,~,~,~] = htp(y,A,s,{'maxNbIter'});
res = norm(x - xStar5);
disp(['Using modified htp, the result is: ', num2str(res)]);
