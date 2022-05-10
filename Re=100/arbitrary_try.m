Re = 100.0;%雷诺数
itermax = 20;
nitermax = 1;
Lambda = 0.000001;%Tikhonov正则化参数
alpha = 1.0;%步长
col = 1./32:1./32:1-1./32;
lc = length(col);
contour = [[col',ones(size(col))'];[col',zeros(size(col))'];[ones(size(col))',col'];[zeros(size(col))',col']];
% col = 1./30:1./30:1-1./30;
contour = [contour;[0,0];[1,0]];
% interior = generate_couple(col,col)';
% f = @(r)(((1-r).^4.*(4*r+1)).*(r<=1)+(r>1).*(0));
% LF = @(r)((100*r.^3 - 240*r.^2 + 180.*r - 40).*(r<=1)+(r>1).*(0));
% Gf = @(r)((4*(r - 1).^4 + 4*(4*r + 1).*(r - 1).^3).*(r<=1)+(r>1).*(0));
load('six_2_sample.mat');
% ci = 0.125;
% R = 2.0;
% f = @(r)(sqrt(r.^2+ci.^2).*(r<R) + 0.*(r>R));
% LF = @(r)(((2*ci.^2+r.^2)./(sqrt(r.^2+ci.^2).^3)).*(r<R)+0.*(r>R));
% Gf = @(r)((1./(sqrt(r.^2+ci.^2))).*(r<R) + 0.*(r>R));
% ci = 0.30;
% R = 2.0;
% f = @(r)(sqrt(r.^2+ci.^2).*(r<R) + 0.*(r>R));
% LF = @(r)(((2*ci.^2+r.^2)./(sqrt(r.^2+ci.^2).^3)).*(r<R)+0.*(r>R));
% Gf = @(r)((1./(sqrt(r.^2+ci.^2))).*(r<R) + 0.*(r>R));
f = @(r)( r.^7);
LF = @(r)(49*r.^5);
Gf = @(r)(7*r.^5);
% BF = @(r)(1./sqrt(r.^2+ci.^2));
% BL = @(r)((-2*ci.^2+r.^2)./(sqrt(r.^2+ci.^2).^5));
% BG = @(r)(-1./sqrt(r.^2+ci.^2).^3);
% col = 1./20:1./20:1-1./20;
% interior = generate_couple(col,col)';
%% 局部加细:
IC = [interior;contour];% 全部插值点
[num,~] = size(IC);
IND = generate_couple(1:num,1:num);
H = reshape(f(sqrt(sum((IC(IND(1,:),:)'-IC(IND(2,:),:)').^2))),num,num);%插值点阵
cond(H)
[vv,ee] = eig(H);
INV = diag(ee);
INV = vv*diag(1./INV)*vv';
clear vv ee
norm(INV*H-eye(size(H)))
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
VH = reshape(f(sqrt(sum((trial(IND(1,:),:)'-IC(IND(2,:),:)').^2))),tri,num);
FH =  reshape(LF(sqrt(sum((trial(IND(1,:),:)'-IC(IND(2,:),:)').^2))),tri,num);
 Gx = reshape((trial(IND(1,:),1)'-IC(IND(2,:),1)').*...
Gf(sqrt(sum((trial(IND(1,:),:)'-IC(IND(2,:),:)').^2))),tri,num);
Gy = reshape((trial(IND(1,:),2)'-IC(IND(2,:),2)').*...
Gf(sqrt(sum((trial(IND(1,:),:)'-IC(IND(2,:),:)').^2))),tri,num);
Ubdy = zeros(length(contour),1);
Ubdy(1:lc,:) = 1.0;
Vbdy = zeros(length(contour),1);
IND = generate_couple(1:length(contour),1:num);
BD = reshape(f(sqrt(sum((contour(IND(1,:),:)'-IC(IND(2,:),:)').^2))),length(contour),num);
clear TND IND
%% 变量部分:
u = rand(num,1);
v = rand(num,1);
p = rand(num,1);
% load('uvp_refined.mat');
H(isnan(H)) = 0;
Gx(isnan(Gx)) = 0;
Gy(isnan(Gy)) = 0;
FH(isnan(FH)) = 0;
VH = eye(size(H));%
 FH = FH*INV;
  Gx = Gx*INV;
  Gy =  Gy*INV;
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
zero_point = find(IC(:,1) == 0 & IC(:,2) == 0);
pressure_cond = p(zero_point);%以求和平均值代替积分
pressure_cond_list = zeros(1,num);
pressure_cond_list(1,zero_point) = 1;
Jall = [Ju;Jv;Jd;[zeros(1,2*num), pressure_cond_list]];
resall = [resu;resv;resdiv;pressure_cond];
rp = zeros(3*num,1);
ins = resall - Jall*rp;
HH = Jall'*Jall + Lambda*eye(3*num);
for nter = 1:nitermax
rp = rp + HH\(Jall'*ins);
ins = resall - Jall*rp;
end
fprintf('epoch = %d,log10(residual norm) = %.5f,norm(gradient) = %f\n',iter,log10(norm(resall)),norm(2*Jall'*resall));
u = u - alpha*rp(1:num);
v = v - alpha*rp(num+1:2*num);
p = p - alpha*rp(2*num+1:3*num);
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
vel_2u = [0.0,-0.03717,-0.04192,-0.04775,-0.06434,-0.10150,-0.15662,-0.21090,-0.20581,-0.13641,0.00332,0.23151,0.68717,0.73722,0.78871,0.84123,1.0];
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
vel_1v = [0 0.09233 0.10091 0.10890 0.12317 0.16077 0.17507 0.17527 0.05454 -0.24533 -0.22445 -0.16914 -0.10313 -0.08864 -0.07391 -0.05906 0];
scatter(vel_1x,vel_1v,'rx')
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
