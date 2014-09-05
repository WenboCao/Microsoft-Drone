nlen=length(VarName1);
p_parameter=1;
dt=1;
a=[0 1 0 0 ;0 0 0 0 ;0 0 0 1;0 0 0 0];
h=[1 0 0 0;0 0 1 0];
F=[1 dt 0 0;0 1 0 0;0 0 1 dt;0 0 0 1];

q1=4;   
q2=4;      

Q1=[q1 0 0 0;0 q1 0 0;0 0 q2 0;0 0 0 q2];
Q=[1/3 1/2 0 0;1/2 1 0 0;0 0 1/3 1/2;0 0 1/2 1]*Q1;
%Kalman Filter preallocation
xapriori = cell(1,nlen);
xaposteriori = cell(1,nlen);
residual = cell(1,nlen);
papriori = cell(1,nlen);
paposteriori = cell(1,nlen);
NEES=zeros(1,nlen);
k = cell(1,nlen);
z=cell(1,nlen);
%Kalman Filter preallocation
xapriori2 = cell(1,nlen);
xaposteriori2 = cell(1,nlen);
residual2 = cell(1,nlen);
papriori2 = cell(1,nlen);
paposteriori2 = cell(1,nlen);
k2 = cell(1,nlen);
z2=cell(1,nlen);
Dkl_0=zeros(1,nlen);
Dkl=zeros(1,nlen);


for j=1:nlen,
z{j}=[x(j);y(j)];
end
R=[0.1446 0;0 0.2470];

paposteriori{1}=[p_parameter 0 0 0;0 p_parameter 0 0;0 0 p_parameter 0;0 0 0 p_parameter];
paposteriori2{1}=[p_parameter 0 0 0;0 p_parameter 0 0;0 0 p_parameter 0;0 0 0 p_parameter];

iZ = z{1};

xaposteriori{1}=[iZ(1);0;iZ(2);0];
xaposteriori2{1}=[iZ(1);0;iZ(2);0];

for j=2:nlen,
    %KF1
    %Predictor equations
    xapriori{j}=F*xaposteriori{j-1};
    papriori{j}=F*paposteriori{j-1}*F'+Q;
    
    % Gating
    residual{j}=z{j}-h*xapriori{j};
    S=h*papriori{j}*h'+R; 

    NEES(j)=residual{j}'*inv(S)*residual{j};
    % An "innovation gate" is used to throw out obviously bad measurements
    if (NEES(j) < 2) 
        k{j}=papriori{j}*h'*inv(S);
        paposteriori{j}=papriori{j}-k{j}*S*k{j}';
        xaposteriori{j}=xapriori{j}+k{j}*residual{j};
    else
        paposteriori{j}=papriori{j};
        xaposteriori{j}=xapriori{j};        
    end
    %KF2
    xapriori2{j}=F*xaposteriori2{j-1};
    papriori2{j}=F*paposteriori2{j-1}*F'+Q;
    
    residual2{j}=z{j}-h*xapriori2{j};
    S=h*papriori2{j}*h'+R; 
    
    k2{j}=papriori2{j}*h'*inv(S);
    paposteriori2{j}=papriori2{j}-k2{j}*S*k2{j}';
    xaposteriori2{j}=xapriori2{j}+k2{j}*residual2{j};
     
    %KL divergency         
    M_a=xaposteriori{j};
    M_b=xapriori2{j};
    C_a=paposteriori{j};
    C_b=papriori2{j};
    
    Dkl_0(j)=trace(C_a\C_b)+((M_a-M_b)'/C_a)*(M_a-M_b)-4-log(det(C_b)/det(C_a));
    Dkl(j)=Dkl_0(j)*0.5;
    Dkl_0(j)=Dkl(j); 
    
    
end

xxx=x';
yyy=y';
%plot 
figure(1)
z_1_es=cell2mat(xaposteriori);
z_x_es=z_1_es(1,:);
z_y_es=z_1_es(3,:);
plot(z_x_es,z_y_es,'r:*');
hold off
title('Posteriori State')

figure(2)
z_1=cell2mat(xapriori2);
z_x=z_1(1,:);
z_y=z_1(3,:);
plot(z_x,z_y,'b*');
title('Estimation')

figure(3)
hold off
plot(z_1_es(1,:), z_1_es(3,:), ':r','LineWidth',2)
hold on
plot(z_1(1,:), z_1(3,:), ':b','LineWidth',2)
legend({'Observation','Estimation'})

figure(4)
hold off
plot(NEES,'r:*','LineWidth',2)
hold on
%plot([1 nlen],2*[1 1],'b:','LineWidth',2)
legend({'Normalized Innovation','Ideal average value'});
%kullback divergency
figure(5)
Mean=mean(Dkl);
hold on
plot(Dkl,'b:*');
hold on
plot([1 nlen],[Mean Mean],'r-');
legend({'KL divergence','mean'})
%observation
figure(6)
zz_1_es=cell2mat(z);
zz_x_es=zz_1_es(1,:);
zz_y_es=zz_1_es(2,:);
plot(zz_x_es,zz_y_es,'r:*');
hold off
title('Observation without NIS')
%plot NEES and KL divergence together
figure(7)
hold off
plot(NEES,'r:*','LineWidth',2)
hold on
plot(Dkl,'b:*');
legend({'Normalized Innovation','KL divergence'});