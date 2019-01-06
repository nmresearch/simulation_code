%parameters
time=20000;
nogenerate=40000; %number of iid variables we first generate to simulate arrivals
beta=0.2;
gamma=0.001;
rho=0.8; %rho for rate matching control
nu1=0.2; %constant for the first square root service rate control
nu2=1; %constant for PSA-based square root service rate control
cl=2*pi/gamma; %length of a cycle
dt=cl/1000;
l=floor(time/dt); %length of time frame
run=10000;
delta=0.95; % confidence interval
defy=10^(-6); % accuracy for tabling functions
defx=10^(-7); % accuracy for tabling functions

%performance measures
W1t_a{1}=zeros(run,l);
W2t_a{1}=zeros(run,l);
W3t_a{1}=zeros(run,l);
Q1t_a{1}=zeros(run,l);
Q2t_a{1}=zeros(run,l);
Q3t_a{1}=zeros(run,l);

%table fn. \Lambda^{-1} for one cycle for generating arrivals and for generating service times under rate matching control
xarray1=(0:defx:cl);
yarray1=xarray1-beta/gamma*(cos(gamma*xarray1)-1);
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

%table fn. M and M^{-1} for one cycle for generating service times under first square root control (M(t) is the integral of \mu(s))
xarray2=(0:defx:cl);
yarray2=nu1*sqrt(1+beta*sin(gamma*xarray2));
yarray2=cumsum(defx*yarray2)+yarray1;
yvec2=(0:defy:yarray2(end));
xvec2=zeros(1,length(yvec2));
i=1;
j=1;
while j<length(xarray2)+1 && i<length(yvec2)+1
y=yvec2(i);
x=xarray2(j);
if y>yarray2(j)
j=j+1;
else
xvec2(i)=x;
i=i+1;
end
end

%table fn. M and M^{-1} for one cycle for generating service times under second square root control
xarray3=(0:defx:cl);
lambda3=1+beta*sin(gamma*xarray3);
yarray3=lambda3/2.*(sqrt(1+nu2./lambda3)+1);
yarray3=cumsum(defx*yarray3);
yvec3=(0:defy:yarray3(end));
xvec3=zeros(1,length(yvec3));
i=1;
j=1;
while j<length(xarray3)+1 && i<length(yvec3)+1
y=yvec3(i);
x=xarray3(j);
if y>yarray3(j)
j=j+1;
else
xvec3(i)=x;
i=i+1;
end
end

%run independent replications
for s=1:run
%simulate the customer arrivals A_k
U1=rand(1,nogenerate);
X=-log(U1); %exponential distribution
T=cumsum(X);
q=floor(T/cl);
r=T-cl*q;
A{1}=xvec1(1+floor(r/defy))+cl*q;
A{1}=A{1}(A{1}<=time); %we only want those arrivals that occur within the time interval we consider
num=length(A{1}); %number of arrivals we consider

%generate service requirements S_k
U2=rand(num,1);
S=-log(U2); %exponential distribution

%simulate begin service time B_k, departure time D_k, service time V_k, waiting time W_k
D1{1}=zeros(1,num+1);
D2{1}=zeros(1,num+1);
D3{1}=zeros(1,num+1);
B1{1}=zeros(1,num);
B2{1}=zeros(1,num);
B3{1}=zeros(1,num);
V1{1}=zeros(1,num);
V2{1}=zeros(1,num);
V3{1}=zeros(1,num);
W1{1}=zeros(1,num);
W2{1}=zeros(1,num);
W3{1}=zeros(1,num);
for i=1:num
B1{1}(i)=max(D1{1}(i),A{1}(i));
B2{1}(i)=max(D2{1}(i),A{1}(i));
B3{1}(i)=max(D3{1}(i),A{1}(i));
sum1=rho*S(i)+B1{1}(i)-beta/gamma*(cos(gamma*B1{1}(i))-1);
V1{1}(i)=xvec1(1+floor(rem(sum1,cl)/defy))+sum1-rem(sum1,cl)-B1{1}(i); %rate matching control
sum2=S(i)+yarray2(1+floor(rem(B2{1}(i),cl)/defx))+floor(B2{1}(i)/cl)*yarray2(end);
V2{1}(i)=xvec2(1+floor(rem(sum2,yvec2(end))/defy))+floor(sum2/yvec2(end))*cl-B2{1}(i); %first square root service rate control
sum3=S(i)+yarray3(1+floor(rem(B3{1}(i),cl)/defx))+floor(B3{1}(i)/cl)*yarray3(end);
V3{1}(i)=xvec3(1+floor(rem(sum3,yvec3(end))/defy))+floor(sum3/yvec3(end))*cl-B3{1}(i); %second square root service rate control
D1{1}(i+1)=B1{1}(i)+V1{1}(i);
D2{1}(i+1)=B2{1}(i)+V2{1}(i);
D3{1}(i+1)=B3{1}(i)+V3{1}(i);
W1{1}(i)=B1{1}(i)-A{1}(i);
W2{1}(i)=B2{1}(i)-A{1}(i);
W3{1}(i)=B3{1}(i)-A{1}(i);
end

%convert to At, Dt, Qt, Wt
At{1}=num*ones(1,l);
i=1;
k=0;
while k<num
if i*dt<A{1}(k+1)
At{1}(i)=k;
i=i+1;
else
k=k+1;
end
end
D1t{1}=num*ones(1,l);
i=1;
k=1;
while i<l+1 && k<num+1
if i*dt<D1{1}(k+1)
D1t{1}(i)=k-1;
i=i+1;
else
k=k+1;
end
end
D2t{1}=num*ones(1,l);
i=1;
k=1;
while i<l+1 && k<num +1
if i*dt<D2{1}(k+1)
D2t{1}(i)=k-1;
i=i+1;
else
k=k+1;
end
end
D3t{1}=num*ones(1,l);
i=1;
k=1;
while i<l+1 && k<num +1
if i*dt<D3{1}(k+1)
D3t{1}(i)=k-1;
i=i+1;
else
k=k+1;
end
end
Q1t{1}=At{1}-D1t{1};
Q2t{1}=At{1}-D2t{1};
Q3t{1}=At{1}-D3t{1};
i=floor(A{1}(1)/dt); %time before the first arrival
W1t{1}=zeros(1,l);
W2t{1}=zeros(1,l);
W3t{1}=zeros(1,l);
i=i+1;
W1t{1}(i:l)=max(W1{1}(At{1}(i:l))+V1{1}(At{1}(i:l))-((i:l)*dt-A{1}(At{1}(i:l))),0);
W2t{1}(i:l)=max(W2{1}(At{1}(i:l))+V2{1}(At{1}(i:l))-((i:l)*dt-A{1}(At{1}(i:l))),0);
W3t{1}(i:l)=max(W3{1}(At{1}(i:l))+V3{1}(At{1}(i:l))-((i:l)*dt-A{1}(At{1}(i:l))),0);

%performance measures
W1t_a{1}(s,:)=W1t{1};
Q1t_a{1}(s,:)=Q1t{1};
W2t_a{1}(s,:)=W2t{1};
Q2t_a{1}(s,:)=Q2t{1};
W3t_a{1}(s,:)=W3t{1};
Q3t_a{1}(s,:)=Q3t{1};
end %independent replication loop ends

%calculate mean and construct confidence intervals from independent experiments
W1t_m{1}=mean(W1t_a{1});
Q1t_m{1}=mean(Q1t_a{1});
W2t_m{1}=mean(W2t_a{1});
Q2t_m{1}=mean(Q2t_a{1});
W3t_m{1}=mean(W3t_a{1});
Q3t_m{1}=mean(Q3t_a{1});
d1_m{1}=sum(W1t_a{1}>0)/run;
d2_m{1}=sum(W2t_a{1}>0)/run;
d3_m{1}=sum(W3t_a{1}>0)/run;
W1t_hw{1}= tinv(1-(1-delta)/2,run-1)*std(W1t_a{1})/sqrt(run); %half-width of the confidence interval
Q1t_hw{1}= tinv(1-(1-delta)/2,run-1)*std(Q1t_a{1})/sqrt(run);
W2t_hw{1}= tinv(1-(1-delta)/2,run-1)*std(W2t_a{1})/sqrt(run);
Q2t_hw{1}= tinv(1-(1-delta)/2,run-1)*std(Q2t_a{1})/sqrt(run);
W3t_hw{1}= tinv(1-(1-delta)/2,run-1)*std(W3t_a{1})/sqrt(run);
Q3t_hw{1}= tinv(1-(1-delta)/2,run-1)*std(Q3t_a{1})/sqrt(run);
W1t_CIu{1}=W1t_m{1}+W1t_hw{1}; %upper bound of confidence interval
W1t_CIl{1}=W1t_m{1}-W1t_hw{1}; %lower bound of confidence interval
Q1t_CIu{1}=Q1t_m{1}+Q1t_hw{1};
Q1t_CIl{1}=Q1t_m{1}-Q1t_hw{1};
W2t_CIu{1}=W2t_m{1}+W2t_hw{1};
W2t_CIl{1}=W2t_m{1}-W2t_hw{1};
Q2t_CIu{1}=Q2t_m{1}+Q2t_hw{1};
Q2t_CIl{1}=Q2t_m{1}-Q2t_hw{1};
W3t_CIu{1}=W3t_m{1}+W3t_hw{1};
W3t_CIl{1}=W3t_m{1}-W3t_hw{1};
Q3t_CIu{1}=Q3t_m{1}+Q3t_hw{1};
Q3t_CIl{1}=Q3t_m{1}-Q3t_hw{1};
