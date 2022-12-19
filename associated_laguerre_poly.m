% credit:  
% David Terr (2022). 
% AssociatedLaguerrePoly.m 
% (https://www.mathworks.com/matlabcentral/fileexchange/4914-associatedlaguerrepoly-m), 
% MATLAB Central File Exchange.

function pol = associated_laguerre_poly(n,a,x)
  if a==0
    pol=polyval(laguerre_poly(n), x);
  else
    Lnk = zeros(n+1,1);
    for m=0:n
        Lnk(n+1-m) = (-1)^m * binomial(a+n,n-m) / factorial(m);
    end
    pol=polyval(Lnk, x);
  end
end