nlen=245;
p_parameter=1;
dt=1;
a=[0 1 0 0 0 0;0 0 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 0;0 0 0 0 0 1;0 0 0 0 0 0];
h=[1 0 0 0 0 0;0 0 1 0 0 0;0 0 0 0 1 0];
F=[1 dt 0 0 0 0;0 1 0 0 0 0;0 0 1 dt 0 0;0 0 0 1 0 0;0 0 0 0 1 dt;0 0 0 0 0 1];

q1=5;
q2=1;
q3=100;
Q1=[q1 0 0 0 0 0;0 q1 0 0 0 0;0 0 q2 0 0 0;0 0 0 q2 0 0;0 0 0 0 q3 0;0 0 0 0 0 q3];
Q=[1/3 1/2 0 0 0 0;1/2 1 0 0 0 0;0 0 1/3 1/2 0 0;0 0 1/2 1 0 0;0 0 0 0 1/3 1/2;0 0 0 0 1/2 1]*Q1;
% SJJ: THIS EXPRESSION IS WRONG - NEEDS TO BE COMPUTED EMPIRICALLY
xapriori = cell(1,nlen);
xaposteriori = cell(1,nlen);
residual = cell(1,nlen);
papriori = cell(1,nlen);
paposteriori = cell(1,nlen);
NEES=zeros(1,nlen);
k = cell(1,nlen);
z=cell(1,nlen);

for j=1:nlen,
z{j}=[x(j);y(j);alt(j)];
end
R=[0.0627 0 0;0 0.2557 0;0 0 1.3855];


%Calculate the process and measurement noise.

% SJJ: NOT PLAUSIBLE;
paposteriori{1}=[p_parameter 0 0 0 0 0;0 p_parameter 0 0 0 0;0 0 p_parameter 0 0 0;0 0 0 p_parameter 0 0;0 0 0 0 p_parameter 0;0 0 0 0 0 p_parameter];

% SJJ: NEEDED BECAUSE CELL ARRAYS UNNECESSARILY COMPLICATE THINGS
iZ = z{1};

xaposteriori{1}=[iZ(1);0;iZ(2);0;iZ(3);0];

for j=2:nlen,
    %Predictor equations
    xapriori{j}=F*xaposteriori{j-1};
    papriori{j}=F*paposteriori{j-1}*F'+Q;
    
    % Gating
    residual{j}=z{j-1}-h*xapriori{j};
    S=h*papriori{j}*h'+R;

    NEES(j)=residual{j}'*inv(S)*residual{j};
    % An "innovation gate" is used to throw out obviously bad measurements
    if (NEES(j) < 10)
        k{j}=papriori{j}*h'*inv(S);
        paposteriori{j}=papriori{j}-k{j}*S*k{j}';
        xaposteriori{j}=xapriori{j}+k{j}*residual{j};
    else
        paposteriori{j}=papriori{j};
        xaposteriori{j}=xapriori{j};
    end
end


%display the state 
j=1:nlen;

figure(1)

z_1_es=cell2mat(xaposteriori);
z_x_es=z_1_es(1,:);
z_y_es=z_1_es(3,:);
z_z_es=z_1_es(5,:);
plot3(z_x_es,z_y_es,z_z_es,'r');
hold off
title('Estimate')

figure(2)
z_1=cell2mat(z);
z_x=z_1(1,:);
z_y=z_1(2,:);
z_z=z_1(3,:);
plot3(z_x,z_y,z_z,'b');
title('Observations')

% Offset
offx=median(z_1_es(1,:));
offy=median(z_1_es(3,:));

figure(3)
hold off
plot(z_1_es(1,:)-offx, z_1_es(3,:)-offy, 'r','LineWidth',2)
hold on
plot(z_1(1,:)-offx, z_1(2,:)-offy, 'b--','LineWidth',2)
legend({'Estimate','Observation'})

%Should be around 2
figure(4)
hold off
% plot(log(NEES),'r','LineWidth',2)
plot(NEES,'r','LineWidth',2)
hold on
plot([1 nlen],log(2)*[1 1],'b--','LineWidth',2)
legend({'Normalized Innovation','Ideal average value'});

