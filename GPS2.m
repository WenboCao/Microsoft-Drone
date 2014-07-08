nlen=20;
dt=1.5;
a=[0 1;0 0];
h=[1 0];
F=[1 dt;0 1];
Q=[10 0;0 10];

R=1;

N=20;
x = cell(1,N);
xapriori = cell(1,N);
xaposteriori = cell(1,N);
residual = cell(1,N);
papriori = cell(1,N);
paposteriori = cell(1,N);
k = cell(1,N);
z={557987.3760211447 557990.291374629 557992.883970632 557998.6481182646 558003.789113737 558010.1207816512 558016.7975010423 558025.2657416272 558038.0070575431 558048.923546594 558058.9831481396 558065.6932202347 51.5218424 558072.32532515 558074.216969771 558080.548678947 558084.955197257 558090.118478903 558093.801718182 558096.728198621};

%Calculate the process and measurement noise.
% w1=randn(1,nlen);
v1=randn(1,nlen);
v=v1*sqrt(R);

%w1=randn([2 1]); %just random a value, need to change after this program run
%w=w1*sqrt(Q);   %How to figure out w when Q is known?
w1=randn;
w=w1*[sqrt(10);sqrt(10)];
%Initial condition on the state, x.
x_0=[557990;1.5];

%Initial guesses for state and a posteriori covariance.
xaposteriori_0=[558000;1.5];
paposteriori_0=[10 0;0 10];

%Calculate the state and the output
% x{1}=a*x_0+w;
% z{1}=h*x{1}+v(1);
%Predictor equations
xapriori{1}=F*xaposteriori_0;
residual{1}=z{1}-h*xapriori{1};
papriori{1}=F*paposteriori_0*F'+Q;
%Corrector equations
k{1}=papriori{1}*h'/(h*papriori{1}*h'+R);
paposteriori{1}=(1-k{1}*h)*papriori{1};
xaposteriori{1}=xapriori{1}+k{1}*residual{1};
%Calculate the rest of the values.
for j=2:nlen,
    %Calculate the state and the output
    %x(j)=a*x(j-1)+w(j);
%     x{j}=a*x{j-1}+w;
%     z{j}=h*x{j}+v(j);
    %Predictor equations
    xapriori{j}=F*xaposteriori{j-1};
    residual{j}=z{j}-h*xapriori{j};
    papriori{j}=F*paposteriori{j-1}*F'+Q;
    %Corrector equations
    k{j}=papriori{j}*h'/(h*papriori{j}*h'+R);
    paposteriori{j}=(1-h*k{j})*papriori{j};
    xaposteriori{j}=xapriori{j}+k{j}*residual{j};
end


%display the state 
j=1:nlen;

%subplot(221);
L=cell2mat(z);
h1=stem(j,L,'g');
hold on 
M=cell2mat(xapriori);
M1=M(1,:);
h2=stem(j+0.5,M1,'b');

N=cell2mat(xaposteriori);
N1=N(1,:);
h3=stem(j+1,N1,'r');

hold off
legend([h1(1),h2(1),h3(1)],'exact','priori','posteriori');
%legend([h1(1)],'exact');
axis([0 50 557980 559000])
% %Plot covariance
% subplot(222);
% L=cell2mat(papriori);
% h1=stem(j,L,'g');
% hold on
% M=cell2mat(paposteriori);
% h2=stem(j+0.5,M,'r');
% hold off
% legend([h1(1),h2(1)],'papriori','paposteriori');