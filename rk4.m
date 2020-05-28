function Xnext=rk4(Xnow,dt,mu,gamma,ifr,N)

k1=deriv(Xnow,mu,gamma,ifr,N)*dt;
k2=deriv(Xnow+k1/2,mu,gamma,ifr,N)*dt;
k3=deriv(Xnow+k2/2,mu,gamma,ifr,N)*dt;
k4=deriv(Xnow+k3,mu,gamma,ifr,N)*dt;
Xnext=Xnow+(k1+2*k2+2*k3+k4)/6;


