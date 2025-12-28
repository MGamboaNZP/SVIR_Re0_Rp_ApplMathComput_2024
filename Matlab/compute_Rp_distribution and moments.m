
% =================================================================================
% Moments and distribution of the bidimensional random variable (Rp^S,Rp^V)
% Based on:
% Gamboa, López-García y Lopez-Herrero (2024) 
% «On the exact and population bi-dimensional reproduction numbers in a 
% stochastic SVIR model with imperfect vaccine», 
% Applied Mathematics and Computation, 468.
% doi:10.1016/J.AMC.2023.128526.
%
% ============================================================

[mvsikp,fprobabilidad,fmasa_marginal_Rp_S,fmasa_marginal_Rp_V]=fmasa_Conjunta_Rp(8,1.2,1,0.05,0.01,4,3,1);

function [mvsikp,fprobabilidad,fmasa_marginal_Rp_S,fmasa_marginal_Rp_V]=fmasa_Conjunta_Rp(N,beta,gamma,h,xi,v0,s0,i0)
% INPUT:
%   v0, s0, i0 : initial conditions
%   N          : population size
%   gamma     : recovery rate
%   xi        : external infection rate
%   beta      : transmission rate
%   h         : vaccine failure probability
%
% OUTPUT:
%   mvsikp: auxiliary probabilities
%   fprobabilidad (Rp_S,Rp_V): joint probability distribution
%   fmasa_marginal_Rp_S: marginal distributin of Rp^S
%   fmasa_marginal_Rp_V: marginal distributin of Rp^V
% =================================================================================


fprobabilidad=zeros(s0+1,v0+1);

%%%%% k=0 y p=0 para (v0,s0,i0)
k=0;
p=0;
v=0;
mvsikp(:,:,v+1,p+1,k+1)=zeros(s0+1,N);
mvsikp(1,:,v+1,p+1,k+1)=1;
for s=1:s0
    for i=1:(N-s-v)
        mvsikp(s+1,i,v+1,p+1,k+1)=i*gamma/q(beta,xi,gamma,h,N,s,i,v);
    end
end 
for v=1:v0
    for s=0:(s0)
       for i=1:(N-s-v)
        mvsikp(s+1,i,v+1,p+1,k+1)=i*gamma/q(beta,xi,gamma,h,N,s,i,v);
        
        end
    end 
end 

fprobabilidad(k+1,p+1)= mvsikp(s0+1,i0,v0+1,p+1,k+1);

for p=1:v0
  for v=0:(p-1)
    mvsikp(:,:,v+1,p+1,k+1)=zeros(s0+1,N);
  end
  
  for v=p:v0
    mvsikp(:,:,v+1,p+1,k+1)=zeros(s0+1,N);
    
    for s=0:(s0)
       for i=1:(N-s-v)
        mvsikp(s+1,i,v+1,p+1,k+1)=(eta(beta,xi,h,N,v,i)*mvsikp(s+1,i+1,v,p,k+1))/q(beta,xi,gamma,h,N,s,i,v);
       end
    end 
  end 
fprobabilidad(k+1,p+1)= mvsikp(s0+1,i0,v0+1,p+1,k+1);
end

for k=1:s0
p=0;
v=0;
mvsikp(:,:,v+1,p+1,k+1)=zeros(s0+1,N);
for s=k:(s0)
       for i=1:(N-s-v)
        mvsikp(s+1,i,v+1,p+1,k+1)=(lambda(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k))/q(beta,xi,gamma,h,N,s,i,v);
       end
end 

for v=1:v0
    mvsikp(:,:,v+1,p+1,k+1)=zeros(s0+1,N);
    
    for s=k:(s0)
       for i=1:(N-s-v)
        mvsikp(s+1,i,v+1,p+1,k+1)=(lambda(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k))/q(beta,xi,gamma,h,N,s,i,v);
       end
    end 
end 
fprobabilidad(k+1,p+1)= mvsikp(s0+1,i0,v0+1,p+1,k+1);

for p=1:v0 
for v=0:(p-1)
mvsikp(:,:,v+1,p+1,k+1)=zeros(s0+1,N);
end
for v=p:v0
    mvsikp(:,:,v+1,p+1,k+1)=zeros(s0+1,N);
    
    for s=k:(s0)
       for i=1:(N-s-v)
        mvsikp(s+1,i,v+1,p+1,k+1)=(eta(beta,xi,h,N,v,i)*mvsikp(s+1,i+1,v,p,k+1)+lambda(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k))/q(beta,xi,gamma,h,N,s,i,v);
       end
    end 
end 
fprobabilidad(k+1,p+1)= mvsikp(s0+1,i0,v0+1,p+1,k+1);
end

end

fmasa_marginal_Rp_V=sum(fprobabilidad);
fmasa_marginal_Rp_S=sum(fprobabilidad,2);

end
function [m1,m2,varianza]=Rp_S(v0,s0,i0,N,beta,gamma,xi,h)

% INPUT:
%   v0, s0, i0 : initial conditions
%   N          : population size
%   gamma     : recovery rate
%   xi        : external infection rate
%   beta      : transmission rate
%   h         : vaccine failure probability
%
% =================================================================================

v=0;
s=0;
Mv_k2=zeros(N+1,N+1);
Mv_k1=zeros(N+1,N+1); 
Mv_k0=ones(N+1,N+1);

s=1;
Mv_k1(s+1,1)=lambda(beta,xi,N,s,0)*Mv_k0(s,2)/q(beta,xi,gamma,h,N,s,0,v);
Mv_k2(s+1,1)=lambda(beta,xi,N,s,0)*2*Mv_k1(s,2)/q(beta,xi,gamma,h,N,s,0,v);
for i=1:(N-v-s)
    Mv_k1(s+1,i+1)=(lambda(beta,xi,N,s,i)*Mv_k0(s,i+2))/q(beta,xi,gamma,h,N,s,i,v);
    Mv_k2(s+1,i+1)=(lambda(beta,xi,N,s,i)*2*Mv_k1(s,i+2))/q(beta,xi,gamma,h,N,s,i,v);
end 
for s=2:s0
    Mv_k1(s+1,1)=(lambda(beta,xi,N,s,0)*(Mv_k0(s,2)+Mv_k1(s,2)))/q(beta,xi,gamma,h,N,s,0,v);
    Mv_k2(s+1,1)=(lambda(beta,xi,N,s,0)*(2*Mv_k1(s,2)+Mv_k2(s,2)))/q(beta,xi,gamma,h,N,s,0,v);
    for i=1:(N-v-s)
          Mv_k1(s+1,i+1)=(lambda(beta,xi,N,s,i)*(Mv_k0(s,i+2)+Mv_k1(s,i+2)))/q(beta,xi,gamma,h,N,s,i,v);
          Mv_k2(s+1,i+1)=(lambda(beta,xi,N,s,i)*(2*Mv_k1(s,i+2)+Mv_k2(s,i+2)))/q(beta,xi,gamma,h,N,s,i,v);
    end
end 
Mvmenos1_k2=Mv_k2;
Mvmenos1_k1=Mv_k1;
Mvmenos1_k0=Mv_k0;
for v=1:v0
    s=0;
    Mv_k2=zeros(N+1,N+1);
    Mv_k1=zeros(N+1,N+1);
    Mv_k0=ones(N+1,N+1);
    for s=1:s0
        Mv_k1(s+1,1)=(lambda(beta,xi,N,s,0)*(Mv_k0(s,2)+Mv_k1(s,2))+eta(beta,xi,h,N,v,0)*Mvmenos1_k1(s+1,2))/q(beta,xi,gamma,h,N,s,0,v);
        Mv_k2(s+1,1)=(lambda(beta,xi,N,s,0)*(2*Mv_k1(s,2)+Mv_k2(s,2))+eta(beta,xi,h,N,v,0)*Mvmenos1_k2(s+1,2))/q(beta,xi,gamma,h,N,s,0,v);
        for i=1:(N-v-s)
             Mv_k1(s+1,i+1)=(lambda(beta,xi,N,s,i)*(Mv_k0(s,i+2)+Mv_k1(s,i+2))+eta(beta,xi,h,N,v,i)*Mvmenos1_k1(s+1,i+2))/q(beta,xi,gamma,h,N,s,i,v);
             Mv_k2(s+1,i+1)=(lambda(beta,xi,N,s,i)*(2*Mv_k1(s,i+2)+Mv_k2(s,i+2))+eta(beta,xi,h,N,v,i)*Mvmenos1_k2(s+1,i+2))/q(beta,xi,gamma,h,N,s,i,v);
        end
    end
   Mvmenos1_k2=Mv_k2; 
   Mvmenos1_k1=Mv_k1;
   Mvmenos1_k0=Mv_k0;
    
end
m1=Mv_k1(s0+1,i0+1);
m2=Mv_k2(s0+1,i0+1);
varianza=m2+m1-m1*m1;


end
function [m1,m2,varianza]=Rp_V(v0,s0,i0,N,beta,gamma,xi,h)

% INPUT:
%   v0, s0, i0 : initial conditions
%   N          : population size
%   gamma     : recovery rate
%   xi        : external infection rate
%   beta      : transmission rate
%   h         : vaccine failure probability
%
% =================================================================================


v=0;
s=0;
Mv_k2=zeros(N+1,N+1);
Mv_k1=zeros(N+1,N+1); 
Mv_k0=ones(N+1,N+1);

v=1;
Mv_k1(v+1,1)=eta(beta,xi,h,N,v,0)*Mv_k0(v,2)/q(beta,xi,gamma,h,N,s,0,v);
Mv_k2(v+1,1)=eta(beta,xi,h,N,v,0)*2*Mv_k1(v,2)/q(beta,xi,gamma,h,N,s,0,v);
for i=1:(N-v-s)
    Mv_k1(v+1,i+1)=(eta(beta,xi,h,N,v,i)*Mv_k0(v,i+2))/q(beta,xi,gamma,h,N,s,i,v);
    Mv_k2(v+1,i+1)=(eta(beta,xi,h,N,v,i)*2*Mv_k1(v,i+2))/q(beta,xi,gamma,h,N,s,i,v);
end 
for v=2:v0
    Mv_k1(v+1,1)=(eta(beta,xi,h,N,v,0)*(Mv_k0(v,2)+Mv_k1(v,2)))/q(beta,xi,gamma,h,N,s,0,v);
    Mv_k2(v+1,1)=(eta(beta,xi,h,N,v,0)*(2*Mv_k1(v,2)+Mv_k2(v,2)))/q(beta,xi,gamma,h,N,s,0,v);
    for i=1:(N-v-s)
          Mv_k1(v+1,i+1)=(eta(beta,xi,h,N,v,i)*(Mv_k0(v,i+2)+Mv_k1(v,i+2)))/q(beta,xi,gamma,h,N,s,i,v);
          Mv_k2(v+1,i+1)=(eta(beta,xi,h,N,v,i)*(2*Mv_k1(v,i+2)+Mv_k2(v,i+2)))/q(beta,xi,gamma,h,N,s,i,v);
    end
end 
Mvmenos1_k2=Mv_k2;
Mvmenos1_k1=Mv_k1;
Mvmenos1_k0=Mv_k0;
for s=1:s0
    v=0;
    Mv_k2=zeros(N+1,N+1);
    Mv_k1=zeros(N+1,N+1);
    Mv_k0=ones(N+1,N+1);
    for v=1:v0
        Mv_k1(v+1,1)=(eta(beta,xi,h,N,v,0)*(Mv_k0(v,2)+Mv_k1(v,2))+lambda(beta,xi,N,s,0)*Mvmenos1_k1(v+1,2))/q(beta,xi,gamma,h,N,s,0,v);
        Mv_k2(v+1,1)=(eta(beta,xi,h,N,v,0)*(2*Mv_k1(v,2)+Mv_k2(v,2))+lambda(beta,xi,N,s,0)*Mvmenos1_k2(v+1,2))/q(beta,xi,gamma,h,N,s,0,v);
        for i=1:(N-v-s)
             Mv_k1(v+1,i+1)=(eta(beta,xi,h,N,v,i)*(Mv_k0(v,i+2)+Mv_k1(v,i+2))+lambda(beta,xi,N,s,i)*Mvmenos1_k1(v+1,i+2))/q(beta,xi,gamma,h,N,s,i,v);
             Mv_k2(v+1,i+1)=(eta(beta,xi,h,N,v,i)*(2*Mv_k1(v,i+2)+Mv_k2(v,i+2))+lambda(beta,xi,N,s,i)*Mvmenos1_k2(v+1,i+2))/q(beta,xi,gamma,h,N,s,i,v);
        end
    end
   Mvmenos1_k2=Mv_k2; 
   Mvmenos1_k1=Mv_k1;
   Mvmenos1_k0=Mv_k0;
    
end
m1=Mv_k1(v0+1,i0+1);
m2=Mv_k2(v0+1,i0+1);
varianza=m2+m1-m1*m1;


end

%% ============================================================
% Auxiliary rate functions
%% ============================================================

function [x]=lambda(beta,xi,N,s,i)
x=beta*s*i/N+xi*s;
end
function [x]=eta(beta,xi,h,N,v,i)
x=beta*v*i*h/N+xi*h*v;
end
function [x]=q(beta,xi,gamma,h,N,s,i,v)
x=beta*s*i/N+xi*s+beta*v*i*h/N+xi*h*v+gamma*i;
end