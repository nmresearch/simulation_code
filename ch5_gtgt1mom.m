%parameters
beta=0.2;
gamma=0.1;
rho=0.8; %rho for rate matching control
lambda=rho;
mu=1;
cl=2*pi/gamma; %length of a cycle
run=40000;
defy=10^(-4); % accuracy for tabling functions
defx=10^(-5); % accuracy for tabling functions
Delta=0.95;
delta=0.001;
b=[0:35000]*delta; %35 rho=0.8; 43 rho0.84; 86 rho0.92; 173 rho0.96; 345 rho0.98; 691 rho0.99

pos=[0:cl/40:cl];

P_positive=zeros(length(pos),run);
EW=zeros(length(pos),run);
EW2=zeros(length(pos),run);
number=floor(10*b(end)/(1/lambda-1));


%table fn. \Lambda_y^{-1} (normalized) for one cycle
xarray1=(0:defx:cl);
yarray=zeros(length(pos),length(xarray1));
yvec1=(0:defy:cl);
xvec=zeros(length(pos),length(yvec1));

for k=1:length(pos)
yarray(k,:)=xarray1+beta/gamma*(cos(gamma*(xarray1-pos(k)))-cos(gamma*pos(k)));
i=1;
j=1;
while j<length(xarray1)+1 && i<length(yvec1)+1
y=yvec1(i);
x=xarray1(j);
if y>yarray(k,j)
j=j+1;
else
xvec(k,i)=x;
i=i+1;
end
end
end
clearvars yarray

for s=1:run
P=zeros(length(b),length(pos));
X=-1/lambda*log(rand(number,1));
Y=-log(rand(number,1));
S1=cumsum(X);
S2=cumsum(Y);
ind=1+floor(rem(S2,cl)/defy);
S3=zeros(number,length(pos));

for k=1:length(pos)
tic
S3(:,k)=xvec(k,ind)'+S2-rem(S2,cl);
S3(:,k)=S3(:,k)-beta/gamma*(cos(gamma*pos(k))-cos(gamma*(S3(:,k)-pos(k))));
Mb=b-beta/gamma*(cos(gamma*(pos(k)+b))-cos(gamma*pos(k)));
i=1;
j=1;
while j<length(b)+1;
if S1(i)-S3(i,k)<Mb(j)
i=i+1;
else
P(j,k)=exp(-(mu-lambda)*(S1(i)-S2(i)));
j=j+1;
end
end
EW(k,s)=sum(delta*P(:,k))+P(length(b),k)/(mu-lambda);
EW2(k,s)=sum(2*P(:,k)'.*b*delta)+2*(b(end)/(mu-lambda)+1/(mu-lambda)^2)*P(length(b),k);
end
P_positive(:,s)=P(1,:)';
end
             
Ppos_final=mean(P_positive')
Ppos_std=std(P_positive')/sqrt(run)
Ppos_hw=norminv(1-(1-Delta)/2)*Ppos_std;
Ppos_CI=[Ppos_final-Ppos_hw,Ppos_final+Ppos_hw]
EW_final=mean(EW')
EW_std=std(EW')/sqrt(run)
EW_hw=norminv(1-(1-Delta)/2)*EW_std;
EW_CI=[EW_final-EW_hw,EW_final+EW_hw]
EW2_final=mean(EW2')
EW2_std=std(EW2')/sqrt(run)
EW2_hw=norminv(1-(1-Delta)/2)*EW2_std;
EW2_CI=[EW2_final-EW2_hw,EW2_final+EW2_hw]
SDW=(EW2_final-EW_final.^2).^0.5
