function [vOut, numIts] = romp_Vershynin(n, Phi, x)
% [vOut] = romp_Vershynin(n, Phi, x)
%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
% d = ambient dimension of the signal v
% N = number of measurements
% n = sparsity level of n
% Phi = N by d measurement matrix
% x = measurement vector (Phi * v)
% vOut = reconstructed signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FUNCTION DESCRIPTION %%%%%%%%%%%%%%%%%
% romp takes parameters as described
% above. Given the sparsity level n and
% the N by d measurement matrix Phi, and
% the measurement vector x = Phi * v, romp
% reconstructs the original signal v.
% This reconstruction is the output.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear r I J J0 u b ix numIts Jvals 
warning off all

N = size(Phi, 1);
d = size(Phi, 2);

% Set residual
r = x;

%Set index set to "empty"
I = zeros(1,1);

%Counter (to be used optionally)
numIts = 0;

%Run ROMP
while length(I)-1 < 2*n && norm(r) > 10^(-6)

   numIts = numIts + 1;

   %Find J, the biggest n coordinates of u
   u = Phi' * r;
   absu = abs(u);
   [b, ix] = sort(absu, 'descend');
   J = ix(1:n);
   Jvals = b(1:n);

   %Find J0, the set of comparable coordinates with maximal energy
   w=1;
   best = -1;
   J0 = zeros(1);
   while w <= n
       first = Jvals(w);
       firstw = w;
       energy = 0;
       while ( w <= n ) && ( Jvals(w) >= 1/2 * first )
           energy = energy + Jvals(w)^2;
           w = w+1;
       end
       if energy > best
           best = energy;
           J0 = J(firstw:w-1);
       end
   end

   %Add J0 to the index set I
   I(length(I)+1: length(I)+length(J0)) = J0;

   %Update the residual
   PhiSubI = Phi(:, I(2));
   for c=3:length(I)
       if ~ismember(I(2:c-1),I(c))
         PhiSubI(:,c-1) = Phi(:, I(c));
       end
   end
   y = lscov(PhiSubI, x);
   r = x - PhiSubI * y;

end % end Run IRA  

vSmall = PhiSubI \ x;
vOut = zeros(d, 1);
for c=2:length(I)
    vOut(I(c)) = vSmall(c-1);
end
               