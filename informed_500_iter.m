%% informed RRT*
Imp=im2bw(imread('ral_jincou.bmp')); 
source=[25 25]; % source position in Y, X format
goal=[384 690]; % goal position in Y, X format
op_dis=759.2*1.01;
% x_I=24; y_I=24;           % 锟斤拷锟矫筹拷始锟斤拷
% x_G=838; y_G=838;       % 锟斤拷锟斤拷目锟斤拷锟?
% x_I=1; y_I=1;           % 锟斤拷锟矫筹拷始锟斤拷
% x_G=700; y_G=700;       % 锟斤拷锟斤拷目锟斤拷锟?    
%% 锟斤拷始锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷业锟斤拷锟斤拷
figure(1);
% ImpRgb=imread('newmap.png');
% % ImpRgb=imread('aaa111a.jpg');
% Imp=rgb2gray(ImpRgb);
imshow(Imp)
hold on
plot(source(1,1), source(1,2), 'ro', 'MarkerSize',10, 'MarkerFaceColor','r');
plot(goal(1,1), goal(1,2), 'go', 'MarkerSize',10, 'MarkerFaceColor','g');% 锟斤拷锟斤拷锟斤拷锟斤拷目锟斤拷锟?
data_500=[];
for i=1:500
    data=informed_rrt_1(source,goal,op_dis,Imp);
    data(:,3)=i;
    data_500=[data_500;data];
end

function s=informed_rrt_1(source,goal,op_dis,Imp)
    Thr=10;             
    Delta= 10;
    x_I=source(1,1); y_I=source(1,2);     
    x_G=goal(1,1); y_G=goal(1,2);     
    T.v(1).x = x_I;        
    T.v(1).y = y_I; 
    T.v(1).xPrev = x_I;   
    T.v(1).yPrev = y_I;
    T.v(1).dist=0;       
    T.v(1).indPrev = 0;     
    xL=size(Imp,2);%锟斤拷图x锟结长锟斤拷
    yL=size(Imp,1);%锟斤拷图y锟结长锟斤拷
    count=1;
    %goal = [x_G,y_G];
    start_goal_dist = 1000000;
    path.pos(1).x = (source(1,1)+goal(1,1))/2;%350;
    path.pos(1).y = (source(1,2)+goal(1,2))/2;%350;
    knt=5000;
    total_dis=zeros(knt,2);
    data_=[];
    Time=tic;
    for iter = 1:5000
        x_rand=[];
        %%=========閲囨牱鎯硏_rand========%%
        if start_goal_dist < 1000000
            while 1
                x_rand(1) = xL*rand; 
                x_rand(2) = yL*rand; 
                if new_node(x_rand(1),x_rand(2),start_goal_dist,source,goal)
                    break;
                end
            end
        else
            if rand < 0.5
                x_rand(1) = xL*rand; 
                x_rand(2) = yL*rand;
            else
                x_rand=goal;
            end
        end
        %%=======瀵绘壘x_near===========%%
        x_near=[];
        min_dist = 1000000;
        near_iter = 1;
        %near_iter_tmp = 1;
        [~,N]=size(T.v);
        for j = 1:N
           x_near(1) = T.v(j).x;
           x_near(2) = T.v(j).y;
           dist = norm(x_rand - x_near);
           if min_dist > dist
               min_dist = dist;
               near_iter = j;
           end
        end
        x_near(1) = T.v(near_iter).x;
        x_near(2) = T.v(near_iter).y;
        %%========鑾峰彇x_new============%%
        x_new=[];
        near_to_rand = [x_rand(1)-x_near(1),x_rand(2)-x_near(2)];
        normlized = near_to_rand / norm(near_to_rand) * Delta;
        x_new = x_near + normlized;
        %plot([x_near(1),x_new(1)],[x_near(2),x_new(2)],'color',[0.7,0.7,0.7], 'Linewidth', 0.5);
        %%=======闅滅妫?娴?===============%%
        if ~collisionChecking_deep(x_near,x_new,Imp) 
           continue;
        end
        %%=======  nearC && chooseParent  =========%%
        near_iter_tmp = near_iter;
        nearptr = [];
        nearcount = 0;
        neardist =  norm(x_new - x_near) + T.v(near_iter_tmp).dist;  
        for j = 2:N
           if j == near_iter_tmp
               continue;
           end
           x_neartmp(1) = T.v(j).x;
           x_neartmp(2) = T.v(j).y;
           dist = norm(x_new - x_neartmp) + T.v(j).dist;
           norm_dist = norm(x_new - x_neartmp);
           if norm_dist < 50
               %nearC
               if collisionChecking_deep(x_neartmp,x_new,Imp)
                    nearcount = nearcount + 1;
                    nearptr(nearcount,1) = j;
                    if neardist > dist 
                        neardist = dist;
                        near_iter = j;
                    end
               end
           end
        end

        x_near(1) = T.v(near_iter).x;
        x_near(2) = T.v(near_iter).y;
        count=count+1;
        %%========灏哫_NEW澧炲姞鍒版爲涓?========%%
        T.v(count).x = x_new(1);
        T.v(count).y = x_new(2); 
        T.v(count).xPrev = x_near(1);     
        T.v(count).yPrev = x_near(2);
        T.v(count).dist= norm(x_new - x_near) + T.v(near_iter).dist;          
        T.v(count).indPrev = near_iter;   
        %%========  rewirte  =========%%
        [M,~] = size(nearptr);
        for k = 1:M
            x_1(1) = T.v(nearptr(k,1)).x;
            x_1(2) = T.v(nearptr(k,1)).y;
            x1_prev(1) = T.v(nearptr(k,1)).xPrev;
            x1_prev(2) = T.v(nearptr(k,1)).yPrev;
            if T.v(nearptr(k,1)).dist >  (T.v(count).dist + norm(x_1-x_new))
                T.v(nearptr(k,1)).dist = T.v(count).dist + norm(x_1-x_new);
                T.v(nearptr(k,1)).xPrev = x_new(1);    
                T.v(nearptr(k,1)).yPrev = x_new(2);
                T.v(nearptr(k,1)).indPrev = count;
    %             plot([x_1(1),x1_prev(1)],[x_1(2),x1_prev(2)],'-w');
    %             hold on;
    %             plot([x_1(1),x_new(1)],[x_1(2),x_new(2)],'-g');
    %             hold on;
            end
        end

    %     plot([x_near(1),x_new(1)],[x_near(2),x_new(2)],'-r');
    %     hold on;
    %     plot(x_new(1),x_new(2),'*r');
    %     hold on;
        if norm(x_new - goal) < 30 || collisionChecking_deep(x_new,goal,Imp)%2*Thr
            if (T.v(count).dist + norm(x_new - goal)) < start_goal_dist
                start_goal_dist = (T.v(count).dist + norm(x_new - goal));
                total_dis(iter,1)=start_goal_dist;
                tEnd = toc(Time); 
                xxxx=[start_goal_dist,tEnd];
                data_=[data_;xxxx];
                %disp(['运行时间: ',num2str(toc)]);
                if length(path.pos) > 2
    %                 for j = 2 : length(path.pos)
    %                     plot([path.pos(j).x; path.pos(j-1).x;], [path.pos(j).y; path.pos(j-1).y], 'w', 'Linewidth', 3);
    %                 end
                end
%                 path.pos = [];
%                 if iter < 5000
%                     path.pos(1).x = x_G; path.pos(1).y = y_G;
%                     path.pos(2).x = T.v(end).x; path.pos(2).y = T.v(end).y;
%                     pathIndex = T.v(end).indPrev; % 锟秸碉拷锟斤拷锟铰凤拷锟?
%                     j=0;
%                     while 1
%                         path.pos(j+3).x = T.v(pathIndex).x;
%                         path.pos(j+3).y = T.v(pathIndex).y;
%                         pathIndex = T.v(pathIndex).indPrev;
%                         if pathIndex == 1
%                             break
%                         end
%                         j=j+1;
%                     end  % 锟斤拷锟秸碉拷锟斤拷莸锟斤拷锟斤拷
%                     path.pos(end+1).x = x_I; path.pos(end).y = y_I; % 锟斤拷锟斤拷锟斤拷路锟斤拷
%     %                 for j = 2:length(path.pos)
%     %                     plot([path.pos(j).x; path.pos(j-1).x;], [path.pos(j).y; path.pos(j-1).y], 'b', 'Linewidth', 3);
%     %                 end
%                 else
%                     disp('Error, no path found!');
%                 end
            end
            %break;
            continue;
        end
        if start_goal_dist<op_dis
            s=data_;
            break;
        end
        %pause(0.01); 
    end
end

function feasible = collisionChecking_deep(startPose, goalPose, map)
    % 起点和终点转为整型坐标
    x1 = round(startPose(1)); y1 = round(startPose(2));
    x2 = round(goalPose(1)); y2 = round(goalPose(2));
    
    % 获取图像尺寸
    [rows, cols] = size(map);
    
    % 边界检查（如果超出图像范围直接判为不可行）
    if x1<1 || x1>cols || y1<1 || y1>rows || x2<1 || x2>cols || y2<1 || y2>rows
        feasible = false;
        return;
    end
    
    % 利用 Bresenham 思想或线性插值得到两点之间所有的坐标点
    % 这里的 max 就是直线的像素步数，完美替代了您那个 0.5 步长的 for 循环
    dx = abs(x2 - x1);
    dy = abs(y2 - y1);
    steps = max(dx, dy); 
    
    if steps == 0
        % 如果两点重合，检查该点是否在障碍物上
        feasible = (map(y1, x1) ~= 0); 
        return;
    end
    
    % 核心提速操作：生成全部 X 和 Y 的向量（一次性生成，不再走循环）
    x_points = round(linspace(x1, x2, steps+1));
    y_points = round(linspace(y1, y2, steps+1));
    
    % 将 X和Y 转化为 MATLAB 图像矩阵的一维线性索引
    % 这一步直接替代了上面的 4次 ceil/floor 判断和 4次 feasiblePoint 调用
    idx = y_points + (x_points - 1) * rows;
    
    % 检查这些像素点中，是否含有障碍物（图像中黑色为0，白色为1或255）
    if any(map(idx) == 0)
        feasible = false;
    else
        feasible = true;
    end
end