%KUNG Ho-Kalman-Kung Algorithm for approximate indification of linear 
%   discrete time-invariant SISO systems based on impulse response samples.
%
%   G = KUNG(h,r,ts) determines an r-th order approximate linear discrete
%   time-invariant SISO realization of the impulse response samples h with
%   sampling time ts.
function G = kung(h,r,ts)
D = h(1);
n = length(h);

l = floor((n-2)/2);

H = hankel(h(2:l+2),h(l+2:2*l+1));
Hb = hankel(h(3:l+3),h(l+3:2*l+2));

[U,S,V] = svd(H);
S = diag(S);
S(r+1:end) = 0;
Cr = diag(sqrt(S(1:r)))*V(:,1:r)';
Or = U(:,1:r)*(diag(sqrt(S(1:r))));

A = Or\Hb/Cr;
B = Cr(:,1);
C = Or(1,:);
G = ss(A,B,C,D,ts);

end