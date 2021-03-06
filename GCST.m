function [x,Jx,GJx,nit] = GCST(J,GJ,x0,pas,epsil,nitmax)
error = 10000;
nit = 0;
while (error > epsil) && (nit < nitmax)
   error = abs(norm(GJ(x0)));
   x = x0 - pas .* GJ(x0);
   nit = nit +1;
   x0=x;
   Jx = J(x);
   GJx = GJ(x);
end

