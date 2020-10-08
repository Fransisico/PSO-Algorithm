%%生成50个粒子
%%求解一个函数的最小值
clc;clear;

Nq = 50;    %%生成的粒子个数
Nd = 1;     %%维数
iter = 3;    %%迭代次数为3
w = 0.8;     %%惯性系数
c1= 0.5;
c2 = 0.5;
limit = [0,50];  %%在0-50中求最优解
vlimit = [-1,1];  %%速度边界为-1-1

p = zeros(Nq,Nd);  %%生成初始位置，Nq个粒子，Nd维度
v = zeros(Nq,Nd);  %%生成初始速度

y_best = inf*ones(Nq,Nd);  %%初始化个体最佳适应度为无穷大
x_best = zeros(Nq,Nd);    %%初始化个体最优解
qy_best = inf;           %%初始化全局最优适应度
qx_best = zeros(Nq,Nd);            %%初始化全局最优解

%%生成50个粒子，初始化每个粒子对应的速度，计算每个粒子的适应度，求得全局最优解和局部最优解
for i = 1:Nq          
    p(i,:) = limit(1,1) + (0.5*rands(1,Nd) + 0.5)*(limit(1,2) - limit(1,1));
    v(i,:) = vlimit(1,1) + (0.5*rands(1,Nd) + 0.5)*(vlimit(1,2) - vlimit(1,1));
    y1 = cots_1(p(i,:));
    
    if  y1 < y_best(i)
        y_best(i) = y1;
        x_best(i,:) = p(i,:);
    end
    
    if  qy_best > y1
        qy_best = y1;
        qx_best(i,:) = p(i,:);
    end
    
%     plot(p(i,:),y1,'o');
%     hold on
     
    
end

 
%%迭代次数
for j = 1:iter
    for i = 1:Nq
        v(i,:) = w*v(i,:) + rand*c1*(x_best(i,:) - p(i,:)) + rand*c2*(qx_best(i,:) - p(i,:));
        
        %位置限幅
        p(i,:) = p(i,:) + v(i,:);
        if p(i,:) > limit(1,2)
            p(i,:) = limit(1,2);
        end
        if p(i,:) < limit(1,1)
            p(i,:) = limit(1,1);
        end
        
        % 速度限幅
        for k = 1:Nd
            if  v(i,k) > vlimit(1,2)
                v(i,k) = vlimit(1,2);
            end
            if  v(i,k) < vlimit(1,1)
                v(i,k) = vlimit(1,1);
            end
        end
        
        
        y1 = cots_1(p(i,:));
        if  y1 < y_best(i)
            y_best(i) = y1;
            x_best = p(i,:);
        end
        
        if qy_best > y1
            qy_best = y1;
            qx_best = p(i,:);
            
        end
        plot(p(i,:),y1,'o');
        hold on
    end
end



qy_best 
qx_best



