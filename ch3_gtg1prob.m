%parameters
beta=0.2;
gamma=1;
rho=0.8;
lambda=rho;
mu=1;
cl=2*pi/gamma; %length of a cycle
run=5000;
b=20;
defy=10^(-4); % accuracy for tabling functions
defx=10^(-5); % accuracy for tabling functions
delta=0.95;
pos=(0:cl/40:cl); %y's in a cycle
P=zeros(length(pos),run);
number=floor(10*b/(1/lambda-1));
X=-1/lambda*log(rand(number,run));
Y=-log(rand(number,run));

%for each position y in a cycle
for r=1:length(pos)

%table fn. \Lambda_y^{-1} (normalized) for one cycle
xarray1=(0:defx:cl);
yarray1=xarray1+beta/gamma*(cos(gamma*(xarray1-pos(r)))-cos(gamma*pos(r)));
yvec1=(0:defy:cl);
xvec1=zeros(1,length(yvec1));
i=1;
j=1;
while j<length(xarray1)+1 && i<length(yvec1)+1
y=yvec1(i);
x=xarray1(j);
if y>yarray1(j)
j=j+1;
else
xvec1(i)=x;
i=i+1;
end
end

for s=1:run
S1=0;
S2=0;
S3=0;
j=0;
while S1-S3<b
j=j+1;
if j>number
Z=rand(1,2);
Z(1)=-1/lambda*log(Z(1));
Z(2)=-log(Z(2));
S1=S1+Z(1);
S2=S2+Z(2);
else
S1=S1+X(j,s);
S2=S2+Y(j,s);
end
S3=xvec1(1+floor(rem(S2,cl)/defy))+S2-rem(S2,cl);
end
P(r,s)=exp(-(mu-lambda)*(S1-S2));
end
end

P_mean=mean(P');
P_std=std(P')/sqrt(run);
P_hw=norminv(1-(1-delta)/2)*P_std;
P_re=P_std./P_mean;
P_approx=exp(-(mu-lambda)*b)*ones(1, length(pos));
A=P_mean./P_approx;
A_approx=rho*exp(-(mu-lambda)*beta/gamma*cos(gamma*pos));
A_LB=rho*exp(-(mu-lambda)*beta/gamma*(cos(gamma*pos)+1));
A_UB=rho*exp(-(mu-lambda)*beta/gamma*(cos(gamma*pos)-1));
          
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
