function [Xk] = DFT(xn, N)
Xk = zeros(1, N);
for n = 1:N
   for k = 1:N
      Xk(n) = Xk(n) + xn(k) * exp(-1j * 2 * pi * (n-1) * (k-1) / N); 
   end
end
end