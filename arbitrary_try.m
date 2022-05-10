Re = 400.0;%雷诺数
itermax = 50;
nitermax = 1;
Lambda = 0.000001;%Tikhonov正则化参数
disthreshold = 0.5;
record_res = inf;
alpha_steps = 9;%余弦退火方法
alpha_list = 0:1/alpha_steps:1-1/alpha_steps;
alpha_list = (cos(pi*alpha_list)+1)/2;
N = 25;
col = 1./41:1./41:1-1./41;
lc = length(col);
contour = [[col',ones(size(col))'];[col',zeros(size(col))'];[ones(size(col))',col'];[zeros(size(col))',col']];
% col = 1./20:1./20:1-1./20;
% col_2 = 1/30:1/15:1-1/30;
contour = [contour;[0,0];[1,0]];
% interior = [generate_couple(col,col)';];%
% f = @(r)(((1-abs(r)).^4.*(4*abs(r)+1)).*(abs(r)<=1)+(abs(r)>1).*(0));
% LF = @(r)((100*abs(r).^3 - 240*abs(r).^2 + 180.*abs(r) - 40).*(abs(r)<=1)+(abs(r)>1).*(0));
% Gf = @(r)((4*(abs(r) - 1).^4 + 4*(4*abs(r) + 1).*(abs(r) - 1).^3).*(abs(r)<=1)+(abs(r)>1).*(0));
load('six_2_sample.mat');
% ci = 0.15;
% R = 2.0;
% f = @(r)(sqrt(r.^2+ci.^2).*(r<R) + 0.*(r>R));
% LF = @(r)(((2*ci.^2+r.^2)./(sqrt(r.^2+ci.^2).^3)).*(r<R)+0.*(r>R));
% Gf = @(r)((1./(sqrt(r.^2+ci.^2))).*(r<R) + 0.*(r>R));
ci = 0.120;
R = 2.0;
f = @(r)(sqrt(r.^2+ci.^2).*(r<R) + 0.*(r>R));
LF = @(r)(((2*ci.^2+r.^2)./(sqrt(r.^2+ci.^2).^3)).*(r<R)+0.*(r>R));
Gf = @(r)((1./(sqrt(r.^2+ci.^2))).*(r<R) + 0.*(r>R));
% f = @(r)(r.^6.*log(r+eps));
% LF = @(r)(30.*r.^4.*log(r+eps) + (6.*r.^5.*log(r+eps) + r.^5)./r + 11.*r.^4);
% Gf = @(r)(6*r.^4.*log(r+eps) + r.^4);
% f = @(r)(1./sqrt(r.^2+ci.^2));
% LF = @(r)((-2*ci.^2+r.^2)./(sqrt(r.^2+ci.^2).^5));
% Gf = @(r)(-1./sqrt(r.^2+ci.^2).^3);
% col = 1./20:1./20:1-1./20;
% interior = generate_couple(col,col)';
%% 局部加细:
IC = [interior;contour];% 全部插值点
[num,~] = size(IC);
IND = generate_couple(1:num,1:num);
H = reshape(f(sqrt(sum((IC(IND(1,:),:)'-IC(IND(2,:),:)').^2))),num,num);%插值点阵

col = 0:0.05:1;
Tr = generate_couple(col,col);
TND = generate_couple(1:length(col)^2,1:num);
Tr = Tr';
TH = reshape(f(sqrt(sum((Tr(TND(1,:),:)'-IC(TND(2,:),:)').^2))),length(col)^2,num);
% trial = [generate_couple(0:1./18:1,0:1./18:1)'];%测试点
% trial = IC(1:length(interior),:);
% trial = trial(randperm(length(interior),length(interior)-10),:);
% trial = [trial;contour];
trial = IC;
size(trial);
[tri,~]=size(trial);
IND = generate_couple(1:tri,1:num);
% DIS = reshape(sqrt(sum((trial(IND(1,:),:)'-IC(IND(2,:),:)').^2)),tri,num);
% DIS(find(DIS<disthreshold)) = -1.0;
% DIS(find(DIS>=0)) = 0;
% DIS(find(DIS<0)) = 1;
% H = H.*DIS;
cond(H)
[vv,ee] = eig(H);
INV = diag(ee);
INV = vv*diag(1./INV)*vv';
clear vv ee
norm(INV*H-eye(size(H)))
% VH = reshape(f(sqrt(sum((trial(IND(1,:),:)'-IC(IND(2,:),:)').^2))),tri,num);
FH =  reshape(LF(sqrt(sum((trial(IND(1,:),:)'-IC(IND(2,:),:)').^2))),tri,num);
 Gx = reshape((trial(IND(1,:),1)'-IC(IND(2,:),1)').*...
Gf(sqrt(sum((trial(IND(1,:),:)'-IC(IND(2,:),:)').^2))),tri,num);
Gy = reshape((trial(IND(1,:),2)'-IC(IND(2,:),2)').*...
Gf(sqrt(sum((trial(IND(1,:),:)'-IC(IND(2,:),:)').^2))),tri,num);
% Gx = Gx.*DIS;
% Gy = Gy.*DIS;
% FH = FH.*DIS;
Ubdy = zeros(length(contour),1);
Ubdy(1:lc,:) = 1.0;
Vbdy = zeros(length(contour),1);
IND = generate_couple(1:length(contour),1:num);
% BD = reshape(f(sqrt(sum((contour(IND(1,:),:)'-IC(IND(2,:),:)').^2))),length(contour),num);
clear TND IND
%% 变量部分:
u = rand(num,1)-0.5;
v = rand(num,1)-0.5;
p = rand(num,1)-0.5;
ori_uvp = [u;v;p];
% load('Re400_ini.mat');
% u = ori_uvp(1:num);
% v = ori_uvp(num+1:2*num);
% p = ori_uvp(2*num+1:3*num);
H(isnan(H)) = 0;
Gx(isnan(Gx)) = 0;
Gy(isnan(Gy)) = 0;
FH(isnan(FH)) = 0;
%% 采用边界内部不同元部分:
% IND = generate_couple(1:length(IC),1:length(contour));
% H(:,end-length(contour)+1:end) = reshape(BF(sqrt(sum((IC(IND(1,:),:)'-contour(IND(2,:),:)').^2))),num,length(contour));
% FH(:,end-length(contour)+1:end) = reshape(BL(sqrt(sum((IC(IND(1,:),:)'-contour(IND(2,:),:)').^2))),num,length(contour));
% Gx(:,end-length(contour)+1:end) = reshape((IC(IND(1,:),1)'-contour(IND(2,:),1)').*...
% BG(sqrt(sum((IC(IND(1,:),:)'-contour(IND(2,:),:)').^2))),num,length(contour));
% Gy(:,end-length(contour)+1:end) = reshape((IC(IND(1,:),2)'-contour(IND(2,:),2)').*...
% BG(sqrt(sum((IC(IND(1,:),:)'-contour(IND(2,:),:)').^2))),num,length(contour));
% IND = generate_couple(1:length(contour),1:length(contour));
% BD(:,end-length(contour)+1:end) = reshape(BF(sqrt(sum((contour(IND(1,:),:)'-contour(IND(2,:),:)').^2))),length(contour),length(contour));
% clear IND
% % 采用01化元:
% INV = inv(H);
% RH = H;
% % 检验病态方程组求解误差:
% T = H\VH';
% fprintf('VH求解误差:%.10f,', norm(H*T-VH'));
% T = H\FH';
% fprintf('FH求解误差:%.10f,', norm(H*T-FH'));
% T = H\Gx';
% fprintf('Gx求解误差:%.10f,', norm(H*T-Gx'));
% T = H\Gy';
% fprintf('Gx求解误差:%.10f \n', norm(H*T-Gy'));
% VH = Gauss_eliminate(H,VH')';
% FH = Gauss_eliminate(H,FH')';
% Gx = Gauss_eliminate(H,Gx')';
% Gy = Gauss_eliminate(H,Gy')';
VH = eye(size(H));%
 FH = FH*INV;
  Gx = Gx*INV;
  Gy =  Gy*INV;
%   FH = FH.*DIS;
%   Gx = Gx.*DIS;
%   Gy = Gy.*DIS;
% BD = BD*INV;
% BD(find(abs(BD)<1.0e-1)) = 0;
% H = eye(size(H));
% % 引入Adam优化器
% alpha = 1.0;
% beta1 = 0.9;
% beta2 = 0.999;
% epsilon = 1.0e-8;
% mm = zeros(3*num,1);% 一阶动量
% mv = zeros(3*num,1);% 二阶动量
%% 迭代部分:
for iter = 1:itermax
u(length(interior)+1:end) = Ubdy;%分离出的边界点
v(length(interior)+1:end) = Vbdy;%分离出的边界点
% 计算残差:
resu = zeros(tri,1);
resv = zeros(tri,1);
% 雅可比阵
Ju = zeros(tri,3*num);
Jv = zeros(tri,3*num);
Ju(1:tri,1:num) =   VH.*sum(((ones(tri,1)*u').*Gx)')'+...
Gx.*sum(((ones(tri,1)*u').*VH)')' + Gy.*sum(((ones(tri,1)*v').*VH)')';
Ju(1:tri,num+1:2*num) =  VH.*sum(((ones(tri,1)*u').*Gy)')';
Jv(1:tri,1:num) = VH.*sum(((ones(tri,1)*v').*Gx)')';
 Jv(1:tri,num+1:2*num) =  VH.*sum(((ones(tri,1)*v').*Gy)')'+...
+ Gy.*sum(((ones(tri,1)*v').*VH)')' + Gx.*sum(((ones(tri,1)*u').*VH)')';
Ju(1:tri,1:num) = Ju(1:tri,1:num) - FH/Re;
Jv(1:tri,num+1:2*num) = Jv(1:tri,num+1:2*num) - FH/Re;
Ju(1:tri,2*num+1:3*num) = Ju(1:tri,2*num+1:3*num) + Gx;
Jv(1:tri,2*num+1:3*num) = Jv(1:tri,2*num+1:3*num) + Gy;
resu = (Gx*u).*(VH*u) +(Gy*u).*(VH*v) - FH*u/Re + Gx*p;
resv = (Gx*v).*(VH*u) +(Gy*v).*(VH*v) - FH*v/Re + Gy*p;
resdiv = Gx*u+Gy*v;
Jd = [Gx,Gy,zeros(size(Gx))];
Ju(:,length(interior)+1:num) = 0;%常量部分为0
Ju(:,num+length(interior)+1:2*num) = 0;%常量部分为0
Jv(:,num+length(interior)+1:2*num) = 0;%常量部分为0
Jv(:,length(interior)+1:num) = 0;%常量部分为0
Jd(:,length(interior)+1:num) = 0;%常量部分为0
Jd(:,num+length(interior)+1:2*num)=0;%常量部分为0
% resbdyu = BD*u-Ubdy;
% resbdyv = BD*v-Vbdy;
zero_point = find(IC(:,1) == 0 & IC(:,2) == 0);
pressure_cond = p(zero_point);%以求和平均值代替积分
pressure_cond_list = zeros(1,num);
pressure_cond_list(1,zero_point) = 1;
% pressure_cond_list = sum(VH);
Jall = [Ju;Jv;Jd;[zeros(1,2*num), pressure_cond_list]];
resall = [resu;resv;resdiv;pressure_cond];
% Jall = [Ju;Jv;Jd;[[BD,zeros(size(BD)),zeros(size(BD))];[zeros(size(BD)),BD,zeros(size(BD))]];[zeros(1,2*num), pressure_cond_list]];
% resall = [resu;resv;resdiv;resbdyu;resbdyv;pressure_cond];
rp = zeros(3*num,1);
ins = resall - Jall*rp;
HH = Jall'*Jall + Lambda*eye(3*num);
for nter = 1:nitermax
% fprintf('epoch = %d tikhonov = %d新的残差norm = %.5f\n',iter,nter,log10(norm(ins)));
rp = rp + HH\(Jall'*ins);
ins = resall - Jall*rp;
end
% rp = 2*Jall'*resall;%梯度
% mm = beta1*mm+(1-beta1)*rp;
% mv = beta2*mv+(1-beta2)*rp.^2;
% rp = (mm./(1-beta1.^iter))./sqrt((mv./(1-beta2.^iter))+epsilon);
% rp = 5*rp/norm(rp);
fprintf('epoch = %d,log10(residual norm) = %.5f,norm(gradient) = %f\n',iter,log10(norm(resall)),norm(2*Jall'*resall));
% rp = rp/norm(rp);
if iter == 1
record_res = norm(resall);
end
%% 由于高雷诺数下产生非数值的震荡,采取谨慎方法
preserved_uvp = [u,v,p];%保存当前极小值
for try_k = 1:alpha_steps
    u = preserved_uvp(:,1);
    v = preserved_uvp(:,2);
    p = preserved_uvp(:,3);
    u = u - alpha_list(try_k)*rp(1:num);
    v = v - alpha_list(try_k)*rp(num+1:2*num);
    p = p - alpha_list(try_k)*rp(2*num+1:3*num);
    u(length(interior)+1:end) = Ubdy;%分离出的边界点
    v(length(interior)+1:end) = Vbdy;%分离出的边界点
resu = (Gx*u).*(VH*u) +(Gy*u).*(VH*v) - FH*u/Re + Gx*p;
resv = (Gx*v).*(VH*u) +(Gy*v).*(VH*v) - FH*v/Re + Gy*p;
resdiv = Gx*u+Gy*v;
pressure_cond = p(zero_point);%以求和平均值代替积分
resall = [resu;resv;resdiv;pressure_cond];
if (norm(resall)<record_res)
record_res = norm(resall);
fprintf('u摄动 = %f,v摄动 = %f,p摄动 = %f \n',norm(preserved_uvp(:,1)-u),norm(preserved_uvp(:,2)-v),norm(preserved_uvp(:,3)-p));
% fprintf('epoch = %d,annealing_epoch = %d,new record:log10(residual norm) = %.5f,norm(gradient) = %f\n',iter,try_k,log10(norm(resall)),norm(2*Jall'*resall));
break;
end
end
if try_k == alpha_steps + 1
fprintf('epoch = %d,annealing failed \n',iter);
end
end
u(length(interior)+1:end) = Ubdy;%分离出的边界点
v(length(interior)+1:end) = Vbdy;%分离出的边界点
[row_Jall,col_Jall] = size(Jall)
% fprintf('为了保证解的唯一性,雅可比矩阵的秩应该达到:%d,实际为:%d\n',(col_Jall - length(contour)*2),rank(Jall));
% % % % % 绘制图像:
% % % % subplot(2,2,1),
% % % %  mesh(0:0.01:1,0:0.01:1,reshape(TH*INV*u,101,101));
% % % %  xlabel('x'),ylabel('y'),
% % % %  title('流速u');
% % % %  subplot(2,2,2),
% % % %  mesh(0:0.01:1,0:0.01:1,reshape(TH*INV*v,101,101));
% % % %  xlabel('x'),ylabel('y'),
% % % %  title('流速v');
%  re = H\TH';
%  fprintf('H|TH转置的求解误差:%.10f\n',norm(H*re-TH'));
vel_2y = [0.0,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5,0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766,1.0];
vel_2u = [0.0,-0.08186,-0.09266,-0.10338,-0.14612,-0.24299,-0.32726,-0.17119,-0.11477,0.02135,0.16256,0.29093,0.55892,0.61756,0.68439,0.75837,1.0];
figure(1),Centreline = reshape(TH*INV*u,21,21)';
plot(Centreline(:,11),0:0.05:1,'b')
hold on;
 scatter(vel_2u,vel_2y,'rx')
Centreline = reshape(TH*INV*v,21,21)';
xlabel('x=0.5上流速分量u的值');ylabel('纵坐标');
legend('Kansa法结果','ghia计算结果');
figure(2),plot(0:0.05:1,Centreline(11,:),'b')
hold on;
vel_1x = [0 0.0625 0.0703 0.0781 0.0938 0.1563 0.2266 0.2344 0.5 0.8047 0.8594 0.9063 0.9453 0.9531 0.9609 0.9688 1.0];
vel_lv = [0.0,0.18360,0.19713,0.20920,0.22965,0.28124,0.30203,0.30174,0.05186,-0.38598,-0.44993,-0.23827,-0.22847,-0.19254,-0.15663,-0.12146,0];
scatter(vel_1x,vel_lv,'rx')
xlabel('横坐标');ylabel('y=0.5上流速分量v的值');
legend('Kansa法结果','ghia计算结果')
figure(3),scatter(interior(:,1),interior(:,2),'filled','r')
hold on;
scatter(contour(:,1),contour(:,2),'filled','b');
xlabel('x');ylabel('y');title('插值点')
figure(4),
 subplot(1,2,1),
 mesh(0:0.05:1,0:0.05:1,reshape(TH*INV*p,21,21)');
 xlabel('x'),ylabel('y'),
 title('压强p');
U = reshape(TH*INV*u,21,21); 
V = reshape(TH*INV*v,21,21); 
 x = 0:0.05:1;
y = 0:0.05:1;
[x,y] = meshgrid(x,y);
subplot(1,2,2),
 [startx,starty] = meshgrid([0.0:0.05:0.2-0.05,0.2:0.2:0.8-0.2,0.8:0.05:1],[0.0:0.05:0.2-0.05,0.2:0.2:0.8-0.2,0.8:0.05:1]);
streamline(x,y,U',V',startx,starty)
