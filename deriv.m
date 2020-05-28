function deriv=deriv(X,mu,gamma,ifr,N)
ncat=length(N);
S=X(1:ncat); I=X(ncat+1:2*ncat); 
R=X(2*ncat+1:3*ncat); D=X(3*ncat+1:4*ncat);

del=mu*(I./N).*S;

dSdt = - del;
dIdt = + del - gamma *          I;
dRdt =       + gamma *(1-ifr).* I;
dDdt =       + gamma *  ifr  .* I;

deriv=[dSdt;dIdt;dRdt;dDdt];
       
