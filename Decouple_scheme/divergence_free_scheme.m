Re = 1000;%雷诺数
itermax = 30;
nitermax = 1;
disthreshold = 0.5;%迎风机制,设定上游半径
record_res = inf;
alpha_steps = 9;%余弦退火方法
alpha_list = 0:1/alpha_steps:1-1/alpha_steps;
alpha_list = (cos(pi*alpha_list)+1)/2;
col = 1./42:1./42:1-1./42;
lc = length(col);
contour = [[col',ones(size(col))'];[col',zeros(size(col))'];[ones(size(col))',col'];[zeros(size(col))',col']];

contour = [contour;[0,0];[1,0]];
load('six_2_sample.mat');
% col = 1./25:1./25:1-1./25;
% interior = [generate_couple(col,col)';];%
ind = find(interior(:,1)<1-0.015 & interior(:,2)<1-0.015 & interior(:,1)>0.015 & interior(:,2)>0.015);%去除距离边界太近的点
interior = interior(ind,:);

ci = 0.10;
R = 2.0;
f = @(r)(sqrt(r.^2+ci.^2).*(r<R) + 0.*(r>R));
LF = @(r)(((2*ci.^2+r.^2)./(sqrt(r.^2+ci.^2).^3)).*(r<R)+0.*(r>R));
Gf = @(r)((1./(sqrt(r.^2+ci.^2))).*(r<R) + 0.*(r>R));

%% 局部加细:
IC = [interior;contour];% 全部插值点
[num,~] = size(IC);
IND = generate_couple(1:num,1:num);
H = reshape(f(sqrt(sum((IC(IND(1,:),:)'-IC(IND(2,:),:)').^2))),num,num);%插值点阵
trial = IC;
col = 0:0.05:1;
Tr = generate_couple(col,col);
TND = generate_couple(1:length(col)^2,1:num);
Tr = Tr';
TH = reshape(f(sqrt(sum((Tr(TND(1,:),:)'-IC(TND(2,:),:)').^2))),length(col)^2,num);
size(trial);
[tri,~]=size(trial);
IND = generate_couple(1:tri,1:num);
DIS = reshape(sqrt(sum((trial(IND(1,:),:)'-IC(IND(2,:),:)').^2)),tri,num);
DIS(find(DIS<disthreshold)) = -1.0;
DIS(find(DIS~=-1.0)) = 0;
DIS = -DIS;
% H = H.*DIS;
fprintf('H条件数=%f\n',cond(H))
[vv,ee] = eig(H);
INV = diag(ee);
INV = vv*diag(1./INV)*vv';
clear vv ee
fprintf('H求逆误差=%f\n',norm(INV*H-eye(size(H))))
% VH = reshape(f(sqrt(sum((trial(IND(1,:),:)'-IC(IND(2,:),:)').^2))),tri,num);
FH =  reshape(LF(sqrt(sum((trial(IND(1,:),:)'-IC(IND(2,:),:)').^2))),tri,num);
 Gx = reshape((trial(IND(1,:),1)'-IC(IND(2,:),1)').*...
Gf(sqrt(sum((trial(IND(1,:),:)'-IC(IND(2,:),:)').^2))),tri,num);
Gy = reshape((trial(IND(1,:),2)'-IC(IND(2,:),2)').*...
Gf(sqrt(sum((trial(IND(1,:),:)'-IC(IND(2,:),:)').^2))),tri,num);

Ubdy = zeros(length(contour),1);
Ubdy(1:lc,:) = 1.0;
Vbdy = zeros(length(contour),1);
IND = generate_couple(1:length(contour),1:num);

clear TND IND
%% 变量部分:
u = zeros(num,1);
v = zeros(num,1);
p = zeros(num,1);
% load('solution_1.mat');
u(num-length(Ubdy)+1:end) = Ubdy;
v(num-length(Vbdy)+1:end) = Vbdy;
H(isnan(H)) = 0;
Gx(isnan(Gx)) = 0;
Gy(isnan(Gy)) = 0;
FH(isnan(FH)) = 0;
VH = eye(size(H));%
 FH = FH*INV;
  Gx = Gx*INV;
  Gy =  Gy*INV;
    Gx = Gx.*DIS;
Gy = Gy.*DIS;
FH = FH.*DIS;
 %% 零散度零边值试探函数空间:
TA = [Gx(:,1:(num-length(Ubdy))),Gy(:,1:(num-length(Vbdy)))];
% TA = [Gx',Gy'];
X = null(TA);
[rx,cx] = size(X);
X = [X(1:rx/2,:);zeros(length(Ubdy),cx);X(rx/2+1:rx,:);zeros(length(Vbdy),cx)];
% X(1:num,:) = INV*X(1:num,:);
% X(num+1:2*num,:) = INV*X(num+1:2*num,:);
% for iter = 1:itermax
%  A_1 = X(1:num,:)'.*(ones(cx,1)*(Gx*u)')-X(1:num,:)'*FH/Re + X(num+1:2*num,:)'.*(ones(cx,1)*(Gx*v)');
%  A_2 = X(num+1:2*num,:)'.*(ones(cx,1)*(Gy*v)')-X(num+1:2*num,:)'*FH/Re + X(1:num,:)'.*(ones(cx,1)*(Gy*u)');
%  W = [A_1(:,1:num-length(Ubdy)),A_2(:,1:num-length(Ubdy));Gx(:,1:num-length(Ubdy)),Gy(:,1:num-length(Ubdy))];
%  rhf = -[A_1(:,num-length(Ubdy)+1:end),A_2(:,num-length(Ubdy)+1:end);Gx(:,num-length(Ubdy)+1:end),Gy(:,num-length(Ubdy)+1:end)]*([Ubdy;Vbdy]);
%  rp = W\rhf;
%  rpu = [rp(1:length(rp)/2);Ubdy];
%  rpv = [rp(length(rp)/2+1:end);Vbdy];
%  fprintf('round = %d,norm(u-u^) = %f,norm(v-v^)=%f \n',iter,norm(u-rpu),norm(v-rpv));
%  u = rpu;
%  v = rpv;
% %  cond(W)
% end
%% 高雷诺数
for iter = 1:itermax
Ju = (-FH'/Re +(Gy'.*(ones(num,1)*v'))+(Gx'.*(ones(num,1)*u'))+(eye(num).*(ones(num,1)*u'*Gx')))*X(1:num,:)+...
(eye(num).*(ones(num,1)*v'*Gx'))*X(num+1:2*num,:);
Jv = (-FH'/Re +(Gx'.*(ones(num,1)*u'))+(Gy'.*(ones(num,1)*v'))+(eye(num).*(ones(num,1)*v'*Gy')))*X(num+1:2*num,:)+...
(eye(num).*(ones(num,1)*u'*Gy'))*X(1:num,:);   
Jdiv = [Gx,Gy];
Ju = Ju';
Jv = Jv';
resuv = X(1:num,:)'*((Gx*u).*u+(Gy*u).*v-FH*u/Re)+X(num+1:2*num,:)'*((Gx*v).*u+(Gy*v).*v-FH*v/Re);
resdiv = Gx*u+Gy*v;
Jall = [Ju(:,1:num-length(Ubdy)),Jv(:,1:num-length(Vbdy));Gx(:,1:num-length(Ubdy)),Gy(:,1:num-length(Vbdy))];
resall = [resuv;resdiv];
rp = Jall\resall;
rpu = [rp(1:length(rp)/2);Ubdy];
rpv = [rp(length(rp)/2+1:end);Vbdy];
preserved_uvp = [u,v];%保存当前极小值
for try_k = 1:alpha_steps
    u = preserved_uvp(:,1);
    v = preserved_uvp(:,2);
    u(1:length(rp)/2) =  u(1:length(rp)/2) - alpha_list(try_k)*rp(1:length(rp)/2);
    v(1:length(rp)/2) =  v(1:length(rp)/2) - alpha_list(try_k)*rp(length(rp)/2+1:end);
    u(length(rp)/2+1:end) = Ubdy;
    v(length(rp)/2+1:end) = Vbdy;
resuv = X(1:num,:)'*((Gx*u).*u+(Gy*u).*v-FH*u/Re)+X(num+1:2*num,:)'*((Gx*v).*u+(Gy*v).*v-FH*v/Re);
resdiv = Gx*u+Gy*v;
resall = [resuv;resdiv];
if (norm(resall)<record_res)
record_res = norm(resall);
fprintf('u摄动 = %f,v摄动 = %f\n',norm(preserved_uvp(:,1)-u),norm(preserved_uvp(:,2)-v));
% fprintf('epoch = %d,annealing_epoch = %d,new record:log10(residual norm) = %.5f,norm(gradient) = %f\n',iter,try_k,log10(norm(resall)),norm(2*Jall'*resall));
break;
end
end
if try_k == alpha_steps 
fprintf('epoch = %d,annealing failed \n',iter);
end

 fprintf('epoch = %d,log10(residual norm) = %.5f,norm(gradient) = %f\n',iter,log10(norm(resall)),norm(2*Jall'*resall));
end
% size(Ju),size(Jv),size(Jdiv)
%% graph:
U = reshape(TH*INV*u,21,21); 
V = reshape(TH*INV*v,21,21); 
 x = 0:0.05:1;
y = 0:0.05:1;
[x,y] = meshgrid(x,y);
subplot(1,2,2),
 [startx,starty] = meshgrid([0.0:0.05:0.2-0.05,0.2:0.2:0.8-0.2,0.8:0.05:1],[0.0:0.05:0.2-0.05,0.2:0.2:0.8-0.2,0.8:0.05:1]);
streamline(x,y,U',V',startx,starty)
