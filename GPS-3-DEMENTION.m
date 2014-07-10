nlen=20;
n_parameter=100;
dt=1.5;
a=[0 1 0 0 0 0;0 0 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 0;0 0 0 0 0 1;0 0 0 0 0 0];
h=[1 0 0 0 0 0;0 0 1 0 0 0;0 0 0 0 1 0];
F=[1 dt 0 0 0 0;0 1 0 0 0 0;0 0 1 dt 0 0;0 0 0 1 0 0;0 0 0 0 1 dt;0 0 0 0 0 1];
T=[n_parameter n_parameter n_parameter n_parameter n_parameter n_parameter];
Q=diag(T);
R=[1 0 0;0 1 0;0 0 1];

N=20;
xapriori = cell(1,N);
xaposteriori = cell(1,N);
residual = cell(1,N);
papriori = cell(1,N);
paposteriori = cell(1,N);
k = cell(1,N);
z={[557987.3760211447;9985674.841317466;39] [557990.291374629;9985666.845917454;39] [557992.883970632;9985656.550789112;39] [557998.6481182646;9985650.686192768;39] [558003.789113737;9985645.40498652;39] [558010.1207816512;9985639.337731386;39] [558016.7975010423;9985637.418589126;39] [558025.2657416272;9985631.824132374;39] [558038.0070575431;9985623.782472568;39] [558048.923546594;9985623.186587939;39] [558058.9831481396;9985620.806814278;39] [558065.6932202347;9985617.36781923;39] [51.5218424;2222;39] [558072.32532515;134545143;39] [558074.216969771;44444;39] [558080.548678947;5555555;39] [558084.955197257;23545465;39] [558090.118478903;2345453245;39] [558093.801718182;45234234;39] [558096.728198621;3423423;39]};

%Calculate the process and measurement noise.
% w1=randn(1,nlen);
v1=randn;
v=v1*sqrt(R);

%w1=randn([2 1]); %just random a value, need to change after this program run
%w=w1*sqrt(Q);   %How to figure out w when Q is known?
w1=randn;
w=w1*[sqrt(n_parameter);sqrt(n_parameter);sqrt(n_parameter);sqrt(n_parameter);sqrt(n_parameter);sqrt(n_parameter)];

%Initial guesses for state and a posteriori covariance.
xaposteriori_0=[558000;1.5;998569;1.5;39;0.1];
paposteriori_0=[10 0 0 0 0 0;0 10 0 0 0 0;0 0 10 0 0 0;0 0 0 10 0 0;0 0 0 0 10 0;0 0 0 0 0 10];

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
    %Predictor equations
    xapriori{j}=F*xaposteriori{j-1};
    residual{j}=z{j}-h*xapriori{j};
    papriori{j}=F*paposteriori{j-1}*F'+Q;
    %Corrector equations
    k{j}=papriori{j}*h'/(h*papriori{j}*h'+R);
    paposteriori{j}=(1-k{j}*h)*papriori{j};
    xaposteriori{j}=xapriori{j}+k{j}*residual{j};
end


%display the state 
j=1:nlen;

z_1=cell2mat(xaposteriori);
z_x=z_1(1,:);
z_y=z_1(2,:);
z_z=z_1(3,:);
figure
plot3(z_x,z_y,z_z);

% hold on 
% M=cell2mat(xapriori);
% M1=M(1,:);
% h2=stem(j+0.5,M1,'b');
% 
% N=cell2mat(xaposteriori);
% N1=N(1,:);
% h3=stem(j+1,N1,'r');
% 
% hold off
% legend([h1(1),h2(1),h3(1)],'exact','priori','posteriori');
% 
% axis([0 50 557980 558500])
