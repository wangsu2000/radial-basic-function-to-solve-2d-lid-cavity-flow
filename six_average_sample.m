%% 六边形均匀插值
function CList = six_average_sample(outer_radius,r)
ZList = zeros(1,1000);
w = exp(pi/3*sqrt(-1));
numberptr = 1;
searchptr = 0;
while(1)
    searchptr = searchptr + 1;
    for j = 0:5
        z = ZList(searchptr)+r*w.^j;
        ind = find(abs(ZList(1:numberptr)-z)<r/2 | abs(z)>=outer_radius*0.975);%判断该点是已经存在于ZList中的点或者是半径outer_radius之外的点
        if(isempty(ind))
        numberptr = numberptr + 1;
        ZList(1,numberptr) = z;
        %fprintf('search = %d ,number = %d\n',searchptr,numberptr);
        end
    end
    if searchptr == numberptr
        break;
    end
end
searchptr = 0;
ZList = [0,ZList(find(abs(ZList)>r/2))];
CList = [real(ZList)',imag(ZList)'];