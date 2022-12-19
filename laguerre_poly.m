% credit:  
% David Terr (2022). 
% AssociatedLaguerrePoly.m 
% (https://www.mathworks.com/matlabcentral/fileexchange/4914-associatedlaguerrepoly-m), 
% MATLAB Central File Exchange.

function Lk = laguerre_poly(n)
    if n==0 
      Lk = 1;
    elseif n==1
      Lk = [-1 1];
    else
      
      Lkm2 = zeros(n+1,1);
      Lkm2(n+1) = 1;
      Lkm1 = zeros(n+1,1);
      Lkm1(n) = -1;
      Lkm1(n+1) = 1;
      for k=2:n
        Lk = zeros(n+1,1);
        for e=n-k+1:n
          Lk(e) = (2*k-1)*Lkm1(e) - Lkm1(e+1) + (1-k)*Lkm2(e);
        end
        Lk(n+1) = (2*k-1)*Lkm1(n+1) + (1-k)*Lkm2(n+1);
        Lk = Lk/k;
        if k<n
           Lkm2 = Lkm1;
           Lkm1 = Lk;
        end 
      end
    end
end
