
% =================================================================================
% Moments and distribution of the bidimensional random variable (Re0^S,Re0^V)
% Based on:
% Gamboa, López-García y Lopez-Herrero (2024) 
% «On the exact and population bi-dimensional reproduction numbers in a 
% stochastic SVIR model with imperfect vaccine», 
% Applied Mathematics and Computation, 468.
% doi:10.1016/J.AMC.2023.128526.
%

% ============================================================
function [mvsikp,fprobabilidad,fmasa_marginal_Re0_S,fmasa_marginal_Re0_V] = fmasa_Conjunta_Re0(N,beta,gamma,h,xi,v0,s0,i0)
fprobabilidad = zeros(s0+1,v0+1);

% k=0 y p=0
k=0; p=0; v=0;
mvsikp(:,:,v+1,p+1,k+1) = zeros(s0+1,N);
mvsikp(1,:,v+1,p+1,k+1) = 1;

for s=1:s0
    i=1;
    mvsikp(s+1,i,v+1,p+1,k+1) = (gamma + lambda_tilda(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k+1))/q(beta,xi,gamma,h,N,s,i,v);
    for i=2:(N-s-v)
        mvsikp(s+1,i,v+1,p+1,k+1) = (gamma + lambda_tilda(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k+1) + (i-1)*gamma*mvsikp(s+1,i-1,v+1,p+1,k+1))/q(beta,xi,gamma,h,N,s,i,v);
    end
end

for v=1:v0
    s=0; i=1;
    mvsikp(s+1,i,v+1,p+1,k+1) = (gamma + eta_tilda(beta,xi,h,N,v,i)*mvsikp(s+1,i+1,v,p+1,k+1))/q(beta,xi,gamma,h,N,s,i,v);
    for i=2:(N-v)
        mvsikp(1,i,v+1,p+1,k+1) = ((i-1)*gamma*mvsikp(1,i-1,v+1,p+1,k+1) + gamma + eta_tilda(beta,xi,h,N,v,i)*mvsikp(s+1,i+1,v,p+1,k+1))/q(beta,xi,gamma,h,N,s,i,v);
    end
    for s=1:s0
        i=1;
        mvsikp(s+1,1,v+1,p+1,k+1) = (eta_tilda(beta,xi,h,N,v,i)*mvsikp(s+1,i+1,v,p+1,k+1) + lambda_tilda(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k+1) + gamma)/q(beta,xi,gamma,h,N,s,1,v);
        for i=2:(N-s-v)
            mvsikp(s+1,i,v+1,p+1,k+1) = (eta_tilda(beta,xi,h,N,v,i)*mvsikp(s+1,i+1,v,p+1,k+1) + lambda_tilda(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k+1) + (i-1)*gamma*mvsikp(s+1,i-1,v+1,p+1,k+1) + gamma)/q(beta,xi,gamma,h,N,s,i,v);
        end
    end
end

fprobabilidad(k+1,p+1) = mvsikp(s0+1,i0,v0+1,p+1,k+1);

% p>0, k=0
for p=1:v0
    for v=0:(p-1)
        mvsikp(:,:,v+1,p+1,k+1) = zeros(s0+1,N);
    end
    for v=p:v0
        mvsikp(:,:,v+1,p+1,k+1) = zeros(s0+1,N);
        s=0;
        mvsikp(s+1,1,v+1,p+1,k+1) = (eta_ast(beta,xi,h,N,v,1)*mvsikp(s+1,2,v,p,k+1) + eta_tilda(beta,xi,h,N,v,1)*mvsikp(s+1,2,v,p+1,k+1))/q(beta,xi,gamma,h,N,s,1,v);
        for i=2:(N-s-v)
            mvsikp(s+1,i,v+1,p+1,k+1) = ((i-1)*gamma*mvsikp(s+1,i-1,v+1,p+1,k+1) + eta_ast(beta,xi,h,N,v,i)*mvsikp(s+1,i+1,v,p,k+1) + eta_tilda(beta,xi,h,N,v,i)*mvsikp(s+1,i+1,v,p+1,k+1))/q(beta,xi,gamma,h,N,s,i,v);
        end
        for s=1:s0
            mvsikp(s+1,1,v+1,p+1,k+1) = (lambda_tilda(beta,xi,N,s,1)*mvsikp(s,2,v+1,p+1,k+1) + eta_ast(beta,xi,h,N,v,1)*mvsikp(s+1,2,v,p,k+1) + eta_tilda(beta,xi,h,N,v,1)*mvsikp(s+1,2,v,p+1,k+1))/q(beta,xi,gamma,h,N,s,1,v);
            for i=2:(N-s-v)
                mvsikp(s+1,i,v+1,p+1,k+1) = (lambda_tilda(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k+1) + (i-1)*gamma*mvsikp(s+1,i-1,v+1,p+1,k+1) + eta_ast(beta,xi,h,N,v,i)*mvsikp(s+1,i+1,v,p,k+1) + eta_tilda(beta,xi,h,N,v,i)*mvsikp(s+1,i+1,v,p+1,k+1))/q(beta,xi,gamma,h,N,s,i,v);
            end
        end
    end
    fprobabilidad(k+1,p+1) = mvsikp(s0+1,i0,v0+1,p+1,k+1);
end

% k>0
for k=1:s0
    p=0; v=0;
    mvsikp(:,:,v+1,p+1,k+1) = zeros(s0+1,N);
    for s=k:s0
        i=1;
        mvsikp(s+1,i,v+1,p+1,k+1) = (lambda_ast(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k) + lambda_tilda(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k+1))/q(beta,xi,gamma,h,N,s,i,v);
        for i=2:(N-s-v)
            mvsikp(s+1,i,v+1,p+1,k+1) = ((i-1)*gamma*mvsikp(s+1,i-1,v+1,p+1,k+1) + lambda_ast(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k) + lambda_tilda(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k+1))/q(beta,xi,gamma,h,N,s,i,v);
        end
    end
    for v=1:v0
        mvsikp(:,:,v+1,p+1,k+1) = zeros(s0+1,N);
        for s=k:s0
            i=1;
            mvsikp(s+1,i,v+1,p+1,k+1) = (eta_tilda(beta,xi,h,N,v,i)*mvsikp(s+1,i+1,v,p+1,k+1) + lambda_ast(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k) + lambda_tilda(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k+1))/q(beta,xi,gamma,h,N,s,i,v);
            for i=2:(N-s-v)
                mvsikp(s+1,i,v+1,p+1,k+1) = (eta_tilda(beta,xi,h,N,v,i)*mvsikp(s+1,i+1,v,p+1,k+1) + (i-1)*gamma*mvsikp(s+1,i-1,v+1,p+1,k+1) + lambda_ast(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k) + lambda_tilda(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k+1))/q(beta,xi,gamma,h,N,s,i,v);
            end
        end
    end
    fprobabilidad(k+1,p+1) = mvsikp(s0+1,i0,v0+1,p+1,k+1);
    for p=1:v0
        for v=0:(p-1)
            mvsikp(:,:,v+1,p+1,k+1) = zeros(s0+1,N);
        end
        for v=p:v0
            mvsikp(:,:,v+1,p+1,k+1) = zeros(s0+1,N);
            for s=k:s0
                i=1;
                mvsikp(s+1,i,v+1,p+1,k+1) = (eta_ast(beta,xi,h,N,v,i)*mvsikp(s+1,i+1,v,p,k+1) + eta_tilda(beta,xi,h,N,v,i)*mvsikp(s+1,i+1,v,p+1,k+1) + lambda_ast(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k) + lambda_tilda(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k+1))/q(beta,xi,gamma,h,N,s,i,v);
                for i=2:(N-s-v)
                    mvsikp(s+1,i,v+1,p+1,k+1) = (eta_ast(beta,xi,h,N,v,i)*mvsikp(s+1,i+1,v,p,k+1) + eta_tilda(beta,xi,h,N,v,i)*mvsikp(s+1,i+1,v,p+1,k+1) + (i-1)*gamma*mvsikp(s+1,i-1,v+1,p+1,k+1) + lambda_ast(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k) + lambda_tilda(beta,xi,N,s,i)*mvsikp(s,i+1,v+1,p+1,k+1))/q(beta,xi,gamma,h,N,s,i,v);
                end
            end
        end
        fprobabilidad(k+1,p+1) = mvsikp(s0+1,i0,v0+1,p+1,k+1);
    end
end

fmasa_marginal_Re0_V = sum(fprobabilidad);
fmasa_marginal_Re0_S = sum(fprobabilidad,2);

end
function [m1,m2,varianza]=Re0_S(v0,s0,i0,N,beta,gamma,xi,h)

v=0;
s=0;
Mv_k2=zeros(N+1,N+1);
Mv_k1=zeros(N+1,N+1); %Deberia de ser (s0+1)x(i0+1)
Mv_k0=ones(N+1,N+1);

s=1;

for i=1:(N-v-s)
   Mv_k1(s+1,i+1)=(lambda_ast(beta,xi,N,s,i)*Mv_k0(s,i+2)+lambda(beta,xi,N,s,i)*Mv_k1(s,i+2)+gamma*(i-1)*Mv_k1(s+1,i))/q(beta,xi,gamma,h,N,s,i,v);
   Mv_k2(s+1,i+1)=(2*lambda_ast(beta,xi,N,s,i)*Mv_k1(s,i+2)+lambda(beta,xi,N,s,i)*Mv_k2(s,i+2)+gamma*(i-1)*Mv_k2(s+1,i))/q(beta,xi,gamma,h,N,s,i,v);
end 
for s=2:s0
    for i=1:(N-v-s)   
     Mv_k1(s+1,i+1)=(lambda_ast(beta,xi,N,s,i)*Mv_k0(s,i+2)+lambda(beta,xi,N,s,i)*Mv_k1(s,i+2)+gamma*(i-1)*Mv_k1(s+1,i))/q(beta,xi,gamma,h,N,s,i,v);
     Mv_k2(s+1,i+1)=(2*lambda_ast(beta,xi,N,s,i)*Mv_k1(s,i+2)+lambda(beta,xi,N,s,i)*Mv_k2(s,i+2)+gamma*(i-1)*Mv_k2(s+1,i))/q(beta,xi,gamma,h,N,s,i,v);
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
        
        for i=1:(N-v-s)
             Mv_k1(s+1,i+1)=(eta(beta,xi,h,N,v,i)*Mvmenos1_k1(s+1,i+2)+lambda_ast(beta,xi,N,s,1)*Mv_k0(s,i+2)+lambda(beta,xi,N,s,i)*Mv_k1(s,i+2)+gamma*(i-1)*Mv_k1(s+1,i))/q(beta,xi,gamma,h,N,s,i,v);
             Mv_k2(s+1,i+1)=(eta(beta,xi,h,N,v,i)*Mvmenos1_k1(s+1,i+2)+2*lambda_ast(beta,xi,N,s,1)*Mv_k1(s,i+2)+lambda(beta,xi,N,s,i)*Mv_k2(s,i+2)+gamma*(i-1)*Mv_k2(s+1,i))/q(beta,xi,gamma,h,N,s,i,v);
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
function [m1,m2,varianza]=Re0_V(v0,s0,i0,N,beta,gamma,xi,h)

s=0;
v=0;
Ms_k2=zeros(N+1,N+1);
Ms_k1=zeros(N+1,N+1); %Deberia de ser (s0+1)x(i0+1)
Ms_k0=ones(N+1,N+1);

v=1;

for i=1:(N-v-s)                          
   Ms_k1(v+1,i+1)=(eta_ast(beta,xi,h,N,v,i)*Ms_k0(v,i+2)+eta(beta,xi,h,N,v,i)*Ms_k1(v,i+2)+gamma*(i-1)*Ms_k1(v+1,i))/q(beta,xi,gamma,h,N,s,i,v);
   Ms_k2(v+1,i+1)=(2*eta_ast(beta,xi,h,N,v,i)*Ms_k1(v,i+2)+eta(beta,xi,h,N,v,i)*Ms_k2(v,i+2)+gamma*(i-1)*Ms_k2(v+1,i))/q(beta,xi,gamma,h,N,s,i,v);
end 
for v=2:v0
    for i=1:(N-v-s)   
     Ms_k1(v+1,i+1)=(eta_ast(beta,xi,h,N,v,i)*Ms_k0(v,i+2)+eta(beta,xi,h,N,v,i)*Ms_k1(v,i+2)+gamma*(i-1)*Ms_k1(v+1,i))/q(beta,xi,gamma,h,N,s,i,v);
     Ms_k2(v+1,i+1)=(2*eta_ast(beta,xi,h,N,v,i)*Ms_k1(v,i+2)+eta(beta,xi,h,N,v,i)*Ms_k2(v,i+2)+gamma*(i-1)*Ms_k2(v+1,i))/q(beta,xi,gamma,h,N,s,i,v);
    end
end 
Msmenos1_k2=Ms_k2;
Msmenos1_k1=Ms_k1;
Msmenos1_k0=Ms_k0;

for s=1:s0
    v=0;
    Ms_k2=zeros(N+1,N+1);
    Ms_k1=zeros(N+1,N+1);
    Ms_k0=ones(N+1,N+1);
    for v=1:v0
        
        for i=1:(N-v-s) 
             Ms_k1(v+1,i+1)=(lambda(beta,xi,N,s,i)*Msmenos1_k1(v+1,i+2)+eta_ast(beta,xi,h,N,v,i)*Ms_k0(v,i+2)+eta(beta,xi,h,N,v,i)*Ms_k1(v,i+2)+gamma*(i-1)*Ms_k1(v+1,i))/q(beta,xi,gamma,h,N,s,i,v);
             Ms_k2(v+1,i+1)=(lambda(beta,xi,N,s,i)*Msmenos1_k2(v+1,i+2)+2*eta_ast(beta,xi,h,N,v,i)*Ms_k1(v,i+2)+eta(beta,xi,h,N,v,i)*Ms_k2(v,i+2)+gamma*(i-1)*Ms_k2(v+1,i))/q(beta,xi,gamma,h,N,s,i,v);
        end
    end
   Msmenos1_k2=Ms_k2; 
   Msmenos1_k1=Ms_k1;
   Msmenos1_k0=Ms_k0;
    
end
m1=Ms_k1(v0+1,i0+1);
m2=Ms_k2(v0+1,i0+1);
varianza=m2+m1-m1*m1;


end
%% ============================================================
% Auxiliary rate functions
%% ============================================================

function [x]=lambda(beta,xi,N,s,i)
x=beta*s*i/N+xi*s;
end
function [x]=lambda_ast(beta,xi,N,s,i)
x=beta*s/N;
end
function [x]=lambda_tilda(beta,xi,N,s,i)
x=(beta*(i-1)/N+xi)*s;
end
function [x]=eta(beta,xi,h,N,v,i)
x=beta*v*i*h/N+xi*h*v;
end
function [x]=eta_ast(beta,xi,h,N,v,i)
x=beta*v*h/N;
end
function [x]=eta_tilda(beta,xi,h,N,v,i)
x=(beta*(i-1)/N+xi)*v*h;
end
function [x]=q(beta,xi,gamma,h,N,s,i,v)
x=beta*s*i/N+xi*s+beta*v*i*h/N+xi*h*v+gamma*i;
end