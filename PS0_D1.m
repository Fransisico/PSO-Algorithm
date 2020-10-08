%%����50������
%%���һ����������Сֵ
clc;clear;

Nq = 50;    %%���ɵ����Ӹ���
Nd = 1;     %%ά��
iter = 3;    %%��������Ϊ3
w = 0.8;     %%����ϵ��
c1= 0.5;
c2 = 0.5;
limit = [0,50];  %%��0-50�������Ž�
vlimit = [-1,1];  %%�ٶȱ߽�Ϊ-1-1

p = zeros(Nq,Nd);  %%���ɳ�ʼλ�ã�Nq�����ӣ�Ndά��
v = zeros(Nq,Nd);  %%���ɳ�ʼ�ٶ�

y_best = inf*ones(Nq,Nd);  %%��ʼ�����������Ӧ��Ϊ�����
x_best = zeros(Nq,Nd);    %%��ʼ���������Ž�
qy_best = inf;           %%��ʼ��ȫ��������Ӧ��
qx_best = zeros(Nq,Nd);            %%��ʼ��ȫ�����Ž�

%%����50�����ӣ���ʼ��ÿ�����Ӷ�Ӧ���ٶȣ�����ÿ�����ӵ���Ӧ�ȣ����ȫ�����Ž�;ֲ����Ž�
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

 
%%��������
for j = 1:iter
    for i = 1:Nq
        v(i,:) = w*v(i,:) + rand*c1*(x_best(i,:) - p(i,:)) + rand*c2*(qx_best(i,:) - p(i,:));
        
        %λ���޷�
        p(i,:) = p(i,:) + v(i,:);
        if p(i,:) > limit(1,2)
            p(i,:) = limit(1,2);
        end
        if p(i,:) < limit(1,1)
            p(i,:) = limit(1,1);
        end
        
        % �ٶ��޷�
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



