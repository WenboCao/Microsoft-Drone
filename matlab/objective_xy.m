nlen=length(VarName1);
p_parameter=1;
dt=1;
a=[0 1 0 0 ;0 0 0 0 ;0 0 0 1;0 0 0 0];
h=[1 0 0 0;0 0 1 0];
F=[1 dt 0 0;0 1 0 0;0 0 1 dt;0 0 0 1];

% q1=9.9;   
% q2=9.9; 
q1=4;   
q2=4; 

Q1=[q1 0 0 0;0 q1 0 0;0 0 q2 0;0 0 0 q2];
Q=[1/3 1/2 0 0;1/2 1 0 0;0 0 1/3 1/2;0 0 1/2 1]*Q1;
%KF1 preallocation
xapriori = cell(1,nlen);
xaposteriori = cell(1,nlen);
residual = cell(1,nlen);
papriori = cell(1,nlen);
paposteriori = cell(1,nlen);
NEES=zeros(1,nlen);
k = cell(1,nlen);
z=cell(1,nlen);
%KF2 preallocation
xapriori2 = cell(1,nlen);
xaposteriori2 = cell(1,nlen);
residual2 = cell(1,nlen);
papriori2 = cell(1,nlen);
paposteriori2 = cell(1,nlen);
k2 = cell(1,nlen);
z2=cell(1,nlen);
Dkl_0=zeros(1,nlen);
Dkl=zeros(1,nlen);
yy=y';

for j=1:nlen,
z{j}=[x(j);y(j)];
end
R=[0.1446 0;0 0.2470];

% SJJ: NOT PLAUSIBLE;
paposteriori{1}=[p_parameter 0 0 0;0 p_parameter 0 0;0 0 p_parameter 0;0 0 0 p_parameter];
paposteriori2{1}=[p_parameter 0 0 0;0 p_parameter 0 0;0 0 p_parameter 0;0 0 0 p_parameter];

% SJJ: NEEDED BECAUSE CELL ARRAYS UNNECESSARILY COMPLICATE THINGS
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
    %keyboard

    NEES(j)=residual{j}'*inv(S)*residual{j};
    % An "innovation gate" is used to throw out obviously bad measurements
    if (NEES(j) < 1) 
        k{j}=papriori{j}*h'*inv(S);
        paposteriori{j}=papriori{j}-k{j}*S*k{j}';
        xaposteriori{j}=xapriori{j}+k{j}*residual{j};
    else
        paposteriori{j}=papriori{j};
        xaposteriori{j}=xapriori{j};
        disp('*');
    end
    %KF2
    %Predictor equations
    xapriori2{j}=F*xaposteriori2{j-1};
    papriori2{j}=F*paposteriori2{j-1}*F'+Q;
    
    % Gating
    residual2{j}=z{j}-h*xapriori2{j};
    S=h*papriori{j}*h'+R; 
    
    k2{j}=papriori2{j}*h'*inv(S);
    paposteriori2{j}=papriori2{j}-k2{j}*S*k2{j}';
    xaposteriori2{j}=xapriori2{j}+k2{j}*residual2{j};
     
    %KL divergency
%     xaposteriori_a=(cell2mat(xaposteriori))';
%     xapriori2(cellfun(@isempty, xapriori2) ) = xaposteriori(1);
%     xapriori2_a=cell2mat(xapriori2)';

    
%     M_0=mean(xaposteriori_a);
%     M_1=mean(xapriori2_a);
%     C_0=cov(xaposteriori_a);
%     C_1=cov(xapriori2_a);
    
    M_a=xaposteriori{j};
    M_b=xaposteriori2{j};
    C_a=paposteriori{j};
    C_b=paposteriori2{j};
    
    %Dkl_0(j)=trace(inv(C_b)*C_a)+(M_b-M_a)'*inv(C_b)*(M_b-M_a)-4-exp(det(C_a)/det(C_b));
    Dkl_0(j)=trace(C_a\C_b)+((M_a-M_b)'/C_a)*(M_a-M_b)-4-log(det(C_b)/det(C_a));
    %Dkl_0(j)=trace(C_1\C_0)+((M_1-M_0)/C_b)*(M_1-M_0)'-4-log(det(C_0)/det(C_1));
    Dkl(j)=Dkl_0(j)*0.5;
    
   
end

% figure(1)
% 
% z_1_es=cell2mat(xapriori);
% z_x_es=z_1_es(1,:);
% z_y_es=z_1_es(3,:);
% plot(z_x_es,z_y_es,'-');
% hold off
% title('Estimate')
% 
% figure(2)
% z_1=cell2mat(z);
% %z_1=cell2mat(xapriori2);
% z_x=z_1(1,:);
% z_y=z_1(2,:);
% plot(z_x,z_y,':');
% title('Observations')
% 
% % Offset
% % offx=median(z_1_es(1,:));
% % offy=median(z_1_es(3,:));
% 
% figure(3)
% hold off
% % plot(z_1_es(1,:)-offx, z_1_es(3,:)-offy, 'r','LineWidth',2)
% plot(z_1_es(1,:), z_1_es(3,:), 'r*','LineWidth',2)
% hold on
% % plot(z_1(1,:)-offx, z_1(2,:)-offy, 'b--','LineWidth',2)
% plot(z_1(1,:), z_1(2,:), 'b*','LineWidth',2)
% legend({'Estimate','Observation'})
% 
% %Should be around 2
% figure(4)
% hold off
% plot(NEES,'r','LineWidth',2)
% hold on
% plot([1 nlen],2*[1 1],'b--','LineWidth',2)
% legend({'Normalized Innovation','Ideal average value'});
% 
% figure(5)
% plot(Dkl,':*');
% hold off
% title('kull-back divergency')
% 
