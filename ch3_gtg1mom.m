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
pos=(0:cl/40:cl);

P_positive=zeros(length(pos),run);
P_estimate=zeros(2,length(b));
EW=zeros(length(pos),run);
EW2=zeros(length(pos),run);
number=floor(10*b(end)/(1/lambda-1));

%table fn. \Lambda_y^{-1} (normalized) for one cycle
xarray1=(0:defx:cl);
yarray1=xarray1+beta/gamma*(cos(gamma*(xarray1-pos(1)))-cos(gamma*pos(1)));
yarray2=xarray1+beta/gamma*(cos(gamma*(xarray1-pos(2)))-cos(gamma*pos(2)));
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

xvec2=zeros(1,length(yvec1));
i=1;
j=1;
while j<length(xarray1)+1 && i<length(yvec1)+1
y=yvec1(i);
x=xarray1(j);
if y>yarray2(j)
j=j+1;
else
xvec2(i)=x;
i=i+1;
end
end

for s=1:run
P=zeros(length(b),2);
X=-1/lambda*log(rand(number,1));
Y=-log(rand(number,1));
S1=cumsum(X);
S2=cumsum(Y);
ind=1+floor(rem(S2,cl)/defy);
S3_1=xvec1(ind)'+S2-rem(S2,cl);
S3_2=xvec2(ind)'+S2-rem(S2,cl);

i=1;
j=1;
while j<length(b)+1;
if S1(i)-S3_1(i)<b(j)
i=i+1;
else
P(j,1)=exp(-(mu-lambda)*(S1(i)-S2(i)));
j=j+1;
end
end

i=1;
j=1;
while j<length(b)+1;
if S1(i)-S3_2(i)<b(j)
i=i+1;
else
P(j,2)=exp(-(mu-lambda)*(S1(i)-S2(i)));
j=j+1;
end
end

P_estimate=P_estimate+P';
P_positive(:,s)=P(1,:)';
EW(1,s)=sum(delta*P(:,1))+P(length(b),1)/(mu-lambda);
EW(2,s)=sum(delta*P(:,2))+P(length(b),2)/(mu-lambda);
EW2(1,s)=sum(2*P(:,1)'.*b*delta)+2*(b(end)/(mu-lambda)+1/(mu-lambda)^2)*P(length(b),1);
EW2(2,s)=sum(2*P(:,2)'.*b*delta)+2*(b(end)/(mu-lambda)+1/(mu-lambda)^2)*P(length(b),2);
end
P_est=P_estimate/run;
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
                                                                                                           
