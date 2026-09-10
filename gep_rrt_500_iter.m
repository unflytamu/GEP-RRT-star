
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Imp=im2bw(imread('ral_jincou.bmp')); 
source=[25 25]; % source position in Y, X format
goal=[384 690]; % goal position in Y, X format
op_dis=759.2*1.01;
Delta= 10;    
mindis=15;
imshow(Imp)
hold on
% plot(source(1,1), source(1,2), 'ro', 'MarkerSize',10, 'MarkerFaceColor','r');
% plot(goal(1,1), goal(1,2), 'go', 'MarkerSize',10, 'MarkerFaceColor','g');% 锟斤拷锟斤拷锟斤拷锟斤拷目锟斤拷锟?
GA_rrt_500=[];
%tic
attempt=1;
while attempt < 2
    s_time=tic;
    %data=rrt_1(source,goal,op_dis,Imp,Delta);%与RRT*结合
    data=informed_rrt_1(source,goal,op_dis,Imp);%与Informed RRT*结合
    %data=Q_rrt_informed_1(source,goal,op_dis,Imp);%与Q_RRT*结合
    init_path=getInitPath(data,mindis);
    GA_path=GA_optmisation(init_path,op_dis,Imp,mindis,Delta);
    tEnd = toc(s_time); 
%     plot(data(:,1), data(:,2), 'r.');
%     plot(data(:,1), data(:,2), 'r');
%     plot(init_path(:,1), init_path(:,2), 'g*');
%     plot(init_path(:,1), init_path(:,2), 'g');
    if isempty(GA_path)
        continue
    end
    if GA_path(1,3)<761%op_dis
    %if GA_path(end,1)<761
        dd=[GA_path(1,3),tEnd];
        GA_rrt_500=[GA_rrt_500;dd];
        attempt=attempt+1;
    end
end
% path=data;
% init_path=getInitPath(path,mindis);
% plot(init_path(:,1), init_path(:,2), 'r.');
% plot(init_path(:,1), init_path(:,2), 'r');
% GA_path=GA_optmisation(init_path,op_dis,Imp,mindis,Delta);
% plot(init_path(:,1), init_path(:,2), 'r.');
% plot(init_path(:,1), init_path(:,2), 'r');

function s = GA_optmisation(init_path, op_dis, Imp, mindis, Delta)
    NP = 30;       % 种群数量
    max_gen = 30;  % 最大进化代数
    pc = 0.9;      % 交叉概率
    pm = 0.9;      % 变异概率
    
    new_pop1 = cell(NP, 1); 
    best_path_ever = init_path; 
    best_dis_ever = cal_dis1(init_path); 
    rec_data=zeros(max_gen,2);
    Time=tic;
    
    %% 1. 种群初始化
    for i = 1:NP
       if i == 1
           path = init_path;
       else
           path = getInitGAPath(init_path, Delta, Imp, 0.8); % 用0.5保证初始多样性
       end
       dis = cal_dis1(path);
       if dis < best_dis_ever
           best_dis_ever = dis;
           best_path_ever = path;
       end
       new_pop1(i, 1) = {path};
    end   
    
    mean_path_value = zeros(1, max_gen);
    min_path_value = zeros(1, max_gen);
    
    %% 2. GA 主循环 (加入局部最优突围)
    stagnation_count = 0; % 记录“停滞”的代数
    best_dis_last = best_dis_ever; % 记录上一代的最优值

    for gen = 1 : max_gen
        % ================== (1) 计算路径长度和适应度 ==================
        path_value = cal_path_value1(new_pop1); 
        [min_dis, m] = min(path_value);         
        
        mean_path_value(gen) = mean(path_value);
        min_path_value(gen) = min_dis;
        
        % 记录和判断是否陷入停滞
        if abs(min_dis - best_dis_last) < 0.5
            stagnation_count = stagnation_count + 1;
        else
            stagnation_count = 0;
        end
        best_dis_last = min_dis;
        
        % 更新全局最优
        if min_dis < best_dis_ever
            best_dis_ever = min_dis;
            best_path_ever = new_pop1{m, 1};
%             if best_dis_ever < op_dis
%                 tEnd = toc(Time); 
%                 best_path_ever(:,3)=best_dis_ever;
%                %best_path_ever(:,4)=best_dis_ever;
%                 s = best_path_ever;
%                 fprintf('第 %d 代: 当前最好路径长度 = %.2f (时间 %0.2f )\n', gen, best_dis_ever, tEnd); 
%                 %disp(['达到目标路径长度，提前结束，运行时间: ', num2str(toc)]);
%                 return;
%             end
        end 
        % 打印当前状态
        %tEnd = toc(Time); 
        %fprintf('第 %d 代: 当前最好路径长度 = %.2f (停滞 %d 代)\n', gen, best_dis_ever, stagnation_count); 
        %fprintf('第 %d 代: 当前最好路径长度 = %.2f (时间 %0.2f )\n', gen, best_dis_ever, tEnd); 
        % 计算适应度
        fit_value = path_value .^ -1;     
        % ================== (2) 遗传操作 ==================
        new_pop2 = selection1(new_pop1, fit_value);
        new_pop2 = crossover1(new_pop2, pc); 
        % ================== (3) 【突围机制】：连续 5 代不优化，触发“大跳跃” ==================
        if stagnation_count >=2 % 停滞 3 代就触发突围，给点压力
            tmp_path = best_path_ever;
            [len, ~] = size(tmp_path);   
            % 确保路径够长才触发突围
            if len > 20
                % 【核心修改】：设定突围跨度。避免跨越半数路径，改成尝试 15%~25% 的总长度。
                min_span = max(3, round(len * 0.10)); % 最小跨度（至少 5 个点）
                max_span = max(8, round(len * 0.20)); % 最大跨度（不超过总长的 25%）
                % 随机生成跨度长度
                span = randi([min_span, max_span]);  
                % 根据生成的跨度，随机安全地选择起始点 idx1
                max_start_idx = len - span - 1;
                if max_start_idx > 2
                    idx1 = randi([2, max_start_idx]);
                    idx2 = idx1 + span; % 终点自动等于起始点+跨度                  
                    p1 = tmp_path(idx1, :);
                    p2 = tmp_path(idx2, :);                 
                    % 检查直线是否可行
                    if collisionChecking_deep(p1, p2, Imp)
                        x_pts = round(linspace(p1(1), p2(1), idx2 - idx1 + 1));
                        y_pts = round(linspace(p1(2), p2(2), idx2 - idx1 + 1));
                        tmp_path(idx1:idx2, 1) = x_pts;
                        tmp_path(idx1:idx2, 2) = y_pts;
                        
                        new_len = cal_dis1(tmp_path);
                        if new_len < best_dis_ever
                            best_dis_ever = new_len;
                            best_path_ever = tmp_path;
%                             if best_dis_ever<op_dis
%                                % tEnd = toc(Time); 
%                                 best_path_ever(:,3)=best_dis_ever;
%                                % best_path_ever(:,4)=best_dis_ever;
%                                 s = best_path_ever;
%                                 tEnd = toc(Time); 
%                                 fprintf('  ---> 路径长度跳变到 %.3f,时间：%.3f\n',  best_dis_ever,tEnd);
%                                 %disp(['范围突围达到目标路径长度提前结束，运行时间: ', num2str(toc)]);
%                                 return;
%                             end
                            % 把这个突变异种强制塞入种群
                            new_pop2{randi(NP), 1} = tmp_path;
                            %fprintf('  ---> 触发范围突围(跨度 %d 点)！路径长度跳变到 %.2f\n', span, best_dis_ever);
                        end
                    end
                end
            end
        end
        % ================== (4) 变异操作 ==================
        % 这里的变异尝试次数已经优化为 5 次
        new_pop2 = mutation1(new_pop2, pm, Imp, 3, m);%
        
        % ================== (5) 精英保留策略 ==================
        new_pop2{m, 1} = best_path_ever;
        
        % ================== (6) 更新种群 ==================
        new_pop1 = new_pop2;
%         tEnd = toc(Time); 
%         cur_dis=cal_dis1(best_path_ever(:,1:2));
%         rec_data(gen,:)=[cur_dis,tEnd];
    end
    best_path_ever(:,3)=cal_dis1(best_path_ever(:,1:2));
    %best_path_ever(:,4)=;
    %disp(['GA 迭代全部结束，运行时间: ', num2str(toc)]);
    s = best_path_ever; 
    %s=rec_data;
end

function path_value = cal_path_value1(pop)
    [px, ~] = size(pop);
    path_value = zeros(1, px); % 预分配
    for i = 1:px
        % 直接调用 cal_dis1 动态计算距离，忽略可能存在的多余列！
        path_value(1, i) = cal_dis1(pop{i, 1});
    end
end

function [path_value] = cal_path_value(pop)
    [n, ~] = size(pop);
    path_value = zeros(1, n);
    for i = 1 : n
        single_pop = pop{i, 1};
        %[m,~ ] = size(single_pop);
        path_value(1, i)=single_pop(1,4);%cal_dis1(single_pop);
%         for j = 1 : m - 1
%             %path_value(1, i) =path_value(1, i)+getDis(single_pop(j,:),single_pop(j+1,:));
%             path_value(1, i)=
%         end
    end
end

function s=cal_dis2(path)
    s=0;
    for i=1:length(path)-1
        dis=norm(path(i+1,1:2)-path(i,1:2),2);
        s=s+dis;
    end
end

function dis = cal_dis1(path)
    % 差分直接求解，0 循环！
    diff_vec = diff(path(:, 1:2), 1, 1); 
    dis = sum(sqrt(sum(diff_vec.^2, 2)));
end

function s=getInitPath(path,mindis)
    ot_path=[path(1,1:2),0];
    for i=1:length(path)-1
        dis=norm(path(i+1,1:2)-path(i,1:2),2);
        theta=mod(atan2(path(i+1,2)-path(i,2),path(i+1,1)-path(i,1)), 2*pi);
        if dis>mindis
            N = ceil(dis / mindis) + 1;
            % 生成插值点
            x_pts = linspace(path(i,1), path(i+1,1), N);
            y_pts = linspace(path(i,2), path(i+1,2), N);
            t_pts = ones(1,N)*theta;
            ot_path=[ot_path;x_pts(2:N)' y_pts(2:N)' t_pts(2:N)'];
        else
            ot_path=[ot_path;path(i+1,1:2),theta];
            continue;
        end
    end
    ot_path(1,3)=ot_path(2,3);
    s=ot_path;
end

function s=rrt_1(source,goal,op_dis,Imp,Delta)
    count=1;
    xL=size(Imp,2);%锟斤拷图x锟结长锟斤拷
    yL=size(Imp,1);%锟斤拷图y锟结长锟斤拷
    x_I=source(1,1); y_I=source(1,2);           % 锟斤拷锟矫筹拷始锟斤拷
    x_G=goal(1,1); y_G=goal(1,2);       % 锟斤拷锟斤拷目锟斤拷锟?
%     Thr=10;                 %锟斤拷锟斤拷目锟斤拷锟斤拷锟街?
%     Delta= 10;              % 锟斤拷锟斤拷锟斤拷展锟斤拷锟斤拷
    start_goal_dist = 1000000;
    path.pos(1).x = x_G;
    path.pos(1).y = y_G;
    knt=5000;
    total_dis=zeros(knt,1);
    %% 锟斤拷锟斤拷锟斤拷始锟斤拷
    T.v(1).x = x_I;         % T锟斤拷锟斤拷锟斤拷要锟斤拷锟斤拷锟斤拷锟斤拷v锟角节点，锟斤拷锟斤拷锟饺帮拷锟斤拷始锟斤拷锟斤拷氲絋锟斤拷锟斤拷锟斤拷
    T.v(1).y = y_I; 
    T.v(1).xPrev = x_I;     % 锟斤拷始锟节碉拷母锟斤拷诘锟斤拷锟饺伙拷锟斤拷浔撅拷锟?
    T.v(1).yPrev = y_I;
    T.v(1).dist=0;          %锟接革拷锟节点到锟矫节碉拷木锟斤拷耄拷锟斤拷锟斤拷取欧锟较撅拷锟斤拷
    T.v(1).indPrev = 0;     %
    %tic
    Time=tic;
    data_=[];
    s=[];
    for iter = 1:knt
        x_rand=[];
        if rand < 0.5 || start_goal_dist<1000
            x_rand(1) = xL*rand; 
            x_rand(2) = yL*rand;
        else
            x_rand=goal;
        end
        %%=======瀵绘壘x_near===========%%
        x_near=[];
        min_dist = 1000000;
        near_iter = 1;
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
        %plot(x_new(1), x_new(2), 'ro', 'MarkerSize',5, 'MarkerFaceColor','m');
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
        %plot([x_near(1),x_new(1)],[x_near(2),x_new(2)],'-r');
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

        %plot([x_near(1),x_new(1)],[x_near(2),x_new(2)],'-b','Linewidth', 0.5);
        %plot([x_near(1),x_new(1)],[x_near(2),x_new(2)],'color',[0.7, 0.7, 0.7],'Linewidth', 0.5);
        if norm(x_new - goal) < 30 || collisionChecking_deep(x_new,goal,Imp)%2*Thr
            daad=T.v(count).dist + norm(x_new - goal);
            if (T.v(count).dist + norm(x_new - goal)) < start_goal_dist
                start_goal_dist = (T.v(count).dist + norm(x_new - goal));
                total_dis(iter)=start_goal_dist;
                %disp([num2str(toc)]);
                tEnd = toc(Time); 
                xxxx=[start_goal_dist,tEnd];
                data_=[data_;xxxx];
                path.pos = [];
                if iter < knt
                    path.pos(1).x = x_G; path.pos(1).y = y_G;
                    path.pos(2).x = T.v(end).x; path.pos(2).y = T.v(end).y;
                    pathIndex = T.v(end).indPrev; % 锟秸碉拷锟斤拷锟铰凤拷锟?
                    j=0;
                    while 1
                        path.pos(j+3).x = T.v(pathIndex).x;
                        path.pos(j+3).y = T.v(pathIndex).y;
                        pathIndex = T.v(pathIndex).indPrev;
                        if pathIndex == 1
                            break
                        end
                        j=j+1;
                    end  % 锟斤拷锟秸碉拷锟斤拷莸锟斤拷锟斤拷
                    path.pos(end+1).x = x_I; path.pos(end).y = y_I; % 锟斤拷锟斤拷锟斤拷路锟斤拷
                else
                    disp('Error, no path found!');
                end
                if ~isempty(path.pos)
                    out_path=[];
                    for jj=1:length(path.pos)
                        ptt(1,1)=path.pos(jj).x;
                        ptt(1,2)=path.pos(jj).y;
                        out_path=[out_path;ptt];
                    end
                    s=flip(out_path);
                    break;
                end
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

function s = getInitGAPath(path, D, obs, p)
    [n, ~] = size(path);
    out = path(1, :);
    
    [rows, cols] = size(obs); % 提前获取地图尺寸，提高防御力
    
    for i = 2:n-1
        % 【核心修改 1】将逻辑翻转，让概率 p 成为“发生变异的概率”
        % 这样我们调用时传 0.8，就意味着有 80% 的节点会变动，大大增加种群多样性！
        if rand < p 
            pt = path(i, :);
            is_ok = false;
            
            % 【核心修改 2】给 while 加上最大 10 次重试限制，避免在狭窄通道里死循环！
            for attempt = 1:10
                num_rand = randi([-D, D], 1, 1);
                pt(1,1) = round(path(i,1) + num_rand * cos(path(i,3) + pi/2));
                pt(1,2) = round(path(i,2) + num_rand * sin(path(i,3) + pi/2));
                pt(1,3) = path(i,3);
                
                % 边界检查，彻底防止 map(x,y) 索引越界
                if pt(1,1) >= 1 && pt(1,1) <= cols && ...
                   pt(1,2) >= 1 && pt(1,2) <= rows && ...
                   obs(pt(1,2), pt(1,1)) ~= 0
                    out = [out; pt];
                    is_ok = true;
                    break; % 变异成功，跳出重试
                end
            end
            
            % 【保护机制】如果尝试了10次仍然撞墙，退一步，直接保留原节点，防止路径断开
            if ~is_ok
                out = [out; path(i, :)];
            end
        else
            out = [out; path(i, :)];
        end
    end
    s = [out; path(n, :)];
end

function s=getInitGAPath1(path,D,obs,p)
     [n,~]=size(path);
     out=path(1,:);
     for i=2:n-1
         if rand<p
            pt=path(i,:);
            out=[out;pt];
         else
             while 1
                 num_rand=randi([0,2*D],1,1)-D;
                 pt(1,1)=round(path(i,1)+num_rand*cos(path(i,3)+pi/2));
                 pt(1,2)=round(path(i,2)+num_rand*sin(path(i,3)+pi/2));
                 pt(1,3)=path(i,3);
                 if obs(pt(1,2),pt(1,1)) ~= 0
                     out=[out;pt];
                     break;
                 end
             end
         end
     end
     s=[out;path(n,:)];
end

function s=getInsertPoints(pt,pt1,D,obpath,mindis,refpath)
    theta=caculaterangle(pt,pt1);
    path=getStraightPath1( pt,pt1, 1,0 );
    lateral=traceerror( path, refpath );
    [n,~]=size(path);
    for i=1:n
        while 1
            if lateral(i)<2
                D=0;
            end
            num_rand=randi([0,2*D],1,1)-D;
            path1(i,1)=round(path(i,1)+num_rand*cos(theta+pi/2));
            path1(i,2)=round(path(i,2)+num_rand*sin(theta+pi/2));
            flag=0;
            for j=1:length(obpath)
                dis=getDis(path1(i,:),obpath(j,:));
                if dis>mindis
                    flag=1;
                    break;
                end
            end
            if flag==1 break;end
%             if map(path1(i,2),path1(i,1))==255 || map(path1(i,2),path1(i,1))==1
%                 break;
%             end
        end
    end
    path1(:,3)=0;
    s=path1;
end

function s=judgeCurve(p1,p2,p3)
    A=p3(1,2)-p1(1,2);
    B=p1(1,1)-p3(1,1);
    C=p1(1,2)*p3(1,1)-p1(1,1)*p3(1,2);
    D=A*p2(1,1)+B*p2(1,2)+C;
    if D>0
        s=1;%右侧
    elseif D<0
        s=2;%左侧
    else
        s=0;%左侧
    end
end

% 变异操作
% 函数说明
% 输入变量：pop：种群，pm：变异概率
% 输出变量：newpop变异以后的种群
function [new_pop] = mutation_shiche(pop, pm, obpath, knt,min_k,mindis)
[px, ~] = size(pop);
new_pop = {};
%knt=5;
for i = 1:px
    % 初始化最大迭代次数
    if i==min_k
        new_pop{i, 1} = pop{i, 1};
        continue;
    end
    max_iteration = 0;
    single_new_pop = pop{i, 1};
    [cnt,~] = size(single_new_pop);
    % single_new_pop_slice初始化
    %single_new_pop_slice = [];
    if(rand < pm)
        while max_iteration<50
%             if max_iteration >30
%                 rand_index=ceil(rand*cnt+60);
%                 if rand_index>cnt-knt 
%                     rand_index = cnt-knt;
%                 end
%             else
%                 rand_index=ceil(rand*cnt-knt);
%             end
            rand_index=ceil(rand*cnt-knt);
            if rand_index<1 rand_index=1;end
            p1=single_new_pop(rand_index,:);
            p2=single_new_pop(rand_index+knt,:);
            %if collisionChecking(p1(:,1:2),p2(:,1:2),G)
            if collisionChecking_shiche(p1(:,1:2),p2(:,1:2),obpath,mindis)
                %dis_l = sqrt((p2(1,2)- p1(1,2))^2 + (p2(1,1) - p1(1,1) )^2);
                a=1/knt;
                for j=1:(knt-1)
                    single_new_pop(rand_index+j,1) = round(p1(1,1) + a * j * (p2(1,1) - p1(1,1)));
                    single_new_pop(rand_index+j,2) = round(p1(1,2) + a * j * (p2(1,2) - p1(1,2)));
                end          
            end
            max_iteration=max_iteration+1;
        end
        new_pop{i, 1} = single_new_pop;
    else
        new_pop{i, 1} = pop{i, 1};
    end

end
end

function s=getPathCurve(Hf,Bf,ss_min)
    n=length(ss_min);
    cur=zeros(n,1);
    for i=1:n
       if i<Bf+1
           bd=0;
       else
           bd=Bf;
       end
       if (n-i)<Hf
           hd=n-i;
       else
           hd=Hf;
       end
       xx=getCurvature(ss_min((i-bd):(i+hd),1:2),2,ss_min(i,1:2));
       cur(i)=xx(2);
    end
    s=cur;
end

function s=getCurvature(fit_path,Dimen,c_point)
%%%%id为当前点
    m=length(fit_path);
    A=zeros(m,Dimen+1);
    for i=1:m
        for j=1:(Dimen+1)
            if j==(Dimen+1)
               A(i,j)=1;
            else
               A(i,j)=fit_path(i,1)^(Dimen+1-j);
            end
        end
    end
    C=fit_path(:,2);
    %D=inv(A'*A)*A'*C;
   % D=(A'*A)\A'*C;
    D=pinv(A'*A)*A'*C;
%     for i=length(D)
%         if abs(D(i))<1e-10
%             D(i)=0;
%         end
%     end
    %获取最近点
    near_point=getNearstPoint(D,c_point,fit_path);
    id=near_point(1,1);
    %cnt=length(D);
    val=mul(D,id);
    theta=atan(val);
    cur=2*D(1);
    cuu=sqrt((1+(2*D(1)*id+D(2))^2)^3);
    cur=cur/cuu;
    s=[theta,cur,D(1),D(2),D(3),near_point(1,1),near_point(1,2)];
%     for i=1:length(fit_path)
%         out_path(i,1)=fit_path(i,1);
%         out_path(i,2)=A(i,:)*D;
%         val=mul(D,fit_path(i,1));
%         out_path(i,3)=atan(val);
%     end
end

function s=getNearstPoint(D,p,path)
    cnt=0;
    x1=path(1,1);
    n=length(path);
    x2=path(n ,1);
    dis=sqrt((path(1,1)-path(2,1))^2+(path(1,2)-path(2,2))^2);
    if(dis<0.1)
        mul1=1;
    else
        mul1=round(dis*10);
    end
    x=linspace(x1,x2,n*mul1);
    n1=length(x);
    y=zeros(1,n1);
    for i=1:n1
        y(i)=D(1)*x(i)*x(i)+D(2)*x(i)+D(3);
    end
    dis1=999;
    for i=1:n1
        d=sqrt((x(i)-p(1,1))^2+(y(i)-p(1,2))^2);
        if d<dis1
            dis1=d;
            cnt=i;
        end    
    end
    % path=[x',y'];
    % plot(path(:,1),path(:,2),'r.');
    % hold on
    s=[x(cnt),y(cnt)];
end
function s=mul(d,x)
  val=0;
  n=length(d);
  for i=1:(n-1)
      val=val+(n-i)*d(i)*x^(n-i-1);
  end
  s=val;

end

function [new_pop] = selection1(pop, fit_value)
    [px, ~] = size(pop);
    
    % 防御性编程：避免适应度全为0导致总概率为0
    total_fit = sum(fit_value);
    if total_fit == 0
        fit_value = ones(px, 1); % 如果全0，随机赋予均等概率
        total_fit = px;
    end
    
    % 1. 计算累积和，并转置成 1行x列 的行向量（为后续矩阵广播做准备）
    edges = cumsum(fit_value)'; 
    
    % 2. 生成 px 个 0~total_fit 之间的随机数（不需要排序）
    r = rand(px, 1) * total_fit; 
    
    % 3. 【核心向量化替换】利用 MATLAB 的隐式扩展（广播）代替 discretize
    % r是 px行1列，edges是1行px列。相减后自动扩展成 px行px列的矩阵。
    % sum(..., 2) 对每一行求和，计算出的数字 +1 就是个体所在的区间索引。
    selected_indices = sum(r > edges, 2) + 1;
    
    % 4. 根据索引直接提取父代个体，生成新种群
    new_pop = pop(selected_indices);
end
% function [new_pop] = selection1(pop, fit_value)
%     [px, ~] = size(pop);
%     
%     % 【安全防御】：处理全0适应度或极小适应度的问题（防止 discretize 边界相等报错）
%     if sum(fit_value) == 0
%         fit_value = ones(px, 1); 
%     end
%     
%     % 1. 计算累积适应度（不再需要归一化到1，直接用原值做区间边界）
%     edges = cumsum(fit_value);
%     total_fit = edges(end);
%     
%     % 2. 生成 [0, total_fit] 区间内的 px 个随机数（无需排序）
%     r = rand(px, 1) * total_fit;
%     
%     % 3. 【核心提速】使用 discretize 将随机数映射到所在的累积区间 -> 返回个体索引
%     selected_indices = discretize(r, [0; edges]);
%     
%     % 4. 根据索引挑选出新的种群
%     new_pop = pop(selected_indices);
%     
%     % 5. 确保输出为列向量的元胞数组
%     new_pop = reshape(new_pop, [], 1); 
% end

function [new_pop] = selection2(pop, fit_value)
    [px, ~] = size(pop);
    
    % 计算累积概率（归一化到0~1区间）
    total_fit = sum(fit_value);
    p_fit_value = cumsum(fit_value) / total_fit;
    
    new_pop = cell(px, 1);
    
    % 直接循环遍历，不需要排序和复杂的双指针
    for i = 1:px
        r = rand(); % 生成一个 0~1 之间的随机数
        % 找到第一个大于 r 的概率位置
        idx = find(p_fit_value >= r, 1); 
        new_pop{i, 1} = pop{idx, 1};
    end
end

% 交叉变换（采用单点交叉）
% 输入变量：pop：父代种群（元胞数组，其中每个元素为 Nx2 的坐标矩阵），pc：交叉的概率
% 输出变量：newpop：交叉后的种群
function [new_pop] = crossover1(pop, pc)
    [px, ~] = size(pop);
    
    % 1. 预分配内存，不要使用 {} 动态扩容，极大提升运行速度
    new_pop = cell(px, 1); 
    
    % 2. 两两配对进行交叉
    for i = 1:2:px-1
        p1 = pop{i, 1};
        p2 = pop{i+1, 1};
        
        % 判断是否进行交叉
        if rand < pc
            % 3. 【核心安全修正】使用 intersect 的 'rows' 选项
            % 之前的 ismember 是按列比较，会导致坐标错乱。用 'rows' 严格按 [x, y] 点对匹配
            [~, idx_in_p1, idx_in_p2] = intersect(p1, p2, 'rows');
            m = length(idx_in_p1);
            
            % 必须有至少3个以上的公共节点才能进行交叉（保证不会切掉起点和终点）
            if m >= 3
                % 4. 随机挑选中间的一个公共交叉点（避免切掉起点/终点）
                r = randi([2, m-1]); 
                cut1 = idx_in_p1(r);
                cut2 = idx_in_p2(r);
                
                % 5. 执行单点交叉，拼接路径
                % 【维度安全】路径为 Nx2 矩阵，必须使用分号 ; 进行垂直拼接
                new_pop{i, 1}   = [p1(1:cut1, :); p2(cut2+1:end, :)];
                new_pop{i+1, 1} = [p2(1:cut2, :); p1(cut1+1:end, :)];
                
                % 交叉成功，直接跳过保留原始路径的步骤
                continue; 
            end
        end
        
        % 概率未命中、或找不到公共点的情况下，保留原个体
        new_pop{i, 1} = p1;
        new_pop{i+1, 1} = p2;
    end
    
    % 6. 如果种群数量是奇数，最后一个个体直接透传
    if mod(px, 2) == 1
        new_pop{px, 1} = pop{px, 1};
    end
end

function [new_pop] = mutation1(pop, pm, G, knt, min_k)
    [px, ~] = size(pop);
    new_pop = cell(px, 1); 
    
    for i = 1:px
        if i == min_k
            new_pop{i, 1} = pop{i, 1};
            continue;
        end
        
        single_new_pop = pop{i, 1};
        [cnt, ~] = size(single_new_pop);
        
        if rand < pm && cnt > knt
            max_start_idx = cnt - knt;
            attempt = 0; % 重置尝试次数
            
            % 【核心修改】变异的失败不直接跳过，而是允许尝试 5 次不同的截断点
            while attempt < 10
                rand_index = randi([1, max_start_idx]);
                p1 = single_new_pop(rand_index, :);
                p2 = single_new_pop(rand_index + knt, :);
                
                if collisionChecking_deep(p1, p2, G)
                    % 如果直线可行，立刻用 linspace 插值替换
                    x_pts = round(linspace(p1(1), p2(1), knt + 1));
                    y_pts = round(linspace(p1(2), p2(2), knt + 1));
                    single_new_pop(rand_index+1 : rand_index+knt-1, 1) = x_pts(2:knt);
                    single_new_pop(rand_index+1 : rand_index+knt-1, 2) = y_pts(2:knt);
                    
                    break; % 变异成功，跳出 while，保留此子代
                end
                attempt = attempt + 1;
            end
        end
        
        new_pop{i, 1} = single_new_pop;
    end
end
% 变异操作
% 输入变量：pop：种群（元胞数组），pm：变异概率，G：环境地图，knt：局部变异长度，min_k：最优个体索引
% 输出变量：new_pop：变异后的种群
function [new_pop] = mutation2(pop, pm, G, knt, min_k)
    [px, ~] = size(pop);
    % 1. 预分配内存（告别动态扩容，提速明显）
    new_pop = cell(px, 1); 
    
    for i = 1:px
        % 2. 精英保留策略：最佳个体不参与变异
        if i == min_k
            new_pop{i, 1} = pop{i, 1};
            continue;
        end
        
        single_new_pop = pop{i, 1};
        [cnt, ~] = size(single_new_pop);
        
        % 3. 判断是否触发变异
        if rand < pm
            % 检查路径长度是否足够做截断变异（边界防御）
            if cnt > knt
                % 安全地生成随机截断位置（消除了原代码的 ceil 计算，防止数组越界）
                max_start_idx = cnt - knt;
                rand_index = randi([1, max_start_idx]);
                
                p1 = single_new_pop(rand_index, :);
                p2 = single_new_pop(rand_index + knt, :);
                
                % 4. 【核心优化1】只需检测一次连线是否可行，废弃无意义的30次 while 循环！
                if collisionChecking_deep(p1, p2, G)
                    % 5. 【核心优化2】利用 linspace 矢量化生成插值点，彻底干掉内部的 for 循环！
                    x_pts = round(linspace(p1(1), p2(1), knt + 1));
                    y_pts = round(linspace(p1(2), p2(2), knt + 1));
                    
                    % 一次性赋值替换中间所有的点（切掉首尾的 p1 和 p2）
                    single_new_pop(rand_index+1 : rand_index+knt-1, 1) = x_pts(2:knt);
                    single_new_pop(rand_index+1 : rand_index+knt-1, 2) = y_pts(2:knt);
                end
            end
        end
        
        new_pop{i, 1} = single_new_pop;
    end
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
            end
        end
        %plot([x_near(1),x_new(1)],[x_near(2),x_new(2)],'color',[0.7, 0.7, 0.7],'Linewidth', 0.5);
        if norm(x_new - goal) < 30 || collisionChecking_deep(x_new,goal,Imp)%2*Thr
            if (T.v(count).dist + norm(x_new - goal)) < start_goal_dist
                start_goal_dist = (T.v(count).dist + norm(x_new - goal));
                total_dis(iter,1)=start_goal_dist;
                tEnd = toc(Time); 
                xxxx=[start_goal_dist,tEnd];
                data_=[data_;xxxx];
                %disp(['运行时间: ',num2str(toc)]);
                path.pos = [];
                if iter < 5000
                    path.pos(1).x = x_G; path.pos(1).y = y_G;
                    path.pos(2).x = T.v(end).x; path.pos(2).y = T.v(end).y;
                    pathIndex = T.v(end).indPrev; % 锟秸碉拷锟斤拷锟铰凤拷锟?
                    j=0;
                    while 1
                        path.pos(j+3).x = T.v(pathIndex).x;
                        path.pos(j+3).y = T.v(pathIndex).y;
                        pathIndex = T.v(pathIndex).indPrev;
                        if pathIndex == 1
                            break
                        end
                        j=j+1;
                    end  % 锟斤拷锟秸碉拷锟斤拷莸锟斤拷锟斤拷
                    path.pos(end+1).x = x_I; path.pos(end).y = y_I; % 锟斤拷锟斤拷锟斤拷路锟斤拷
                else
                    disp('Error, no path found!');
                end
                if ~isempty(path.pos)
                    out_path=[];
                    for jj=1:length(path.pos)
                        ptt(1,1)=path.pos(jj).x;
                        ptt(1,2)=path.pos(jj).y;
                        out_path=[out_path;ptt];
                    end
                    s=flip(out_path);
                    return;
                end
            end
        end
    end
end

function s=Q_rrt_informed_1(source,goal,op_dis,Imp)
    x_I = source(1,1); y_I = source(1,2);
    x_G = goal(1,1);   y_G = goal(1,2);
    xL = size(Imp,2);
    yL = size(Imp,1);
    count = 1;
    start_goal_dist = 1000000;
    path.pos(1).x = (source(1,1)+goal(1,1))/2;
    path.pos(1).y = (source(1,2)+goal(1,2))/2;
    knt = 5000;
    total_dis = zeros(knt,2);
    Thr = 10;                  % goal tolerance
    Delta = 10;                % extension step size
    Depth = 2;                 % Q-RRT* depth parameter (number of ancestors to consider)

    %% 初始化树 T
    T.v(1).x = x_I;
    T.v(1).y = y_I;
    T.v(1).xPrev = x_I;
    T.v(1).yPrev = y_I;
    T.v(1).dist = 0;
    T.v(1).indPrev = 0;
    Time = tic;
    data_ = [];
    for iter = 1:5000
        x_rand = [];
        %% ========= 采样 (Informed 采样策略) ========= %%
        if start_goal_dist < 1000000
            while 1
                x_rand(1) = xL * rand;
                x_rand(2) = yL * rand;
                if new_node(x_rand(1), x_rand(2), start_goal_dist, source, goal)
                    break;
                end
            end
        else
            if rand < 0.5
                x_rand(1) = xL * rand;
                x_rand(2) = yL * rand;
            else
                x_rand = goal;
            end
        end

        %% ========= 寻找最近邻 x_near ========= %%
        min_dist = 1000000;
        near_iter = 1;
        N = length(T.v);
        for j = 1:N
            x_near_tmp(1) = T.v(j).x;
            x_near_tmp(2) = T.v(j).y;
            dist = norm(x_rand - x_near_tmp);
            if min_dist > dist
                min_dist = dist;
                near_iter = j;
            end
        end
        x_near = [T.v(near_iter).x, T.v(near_iter).y];

        %% ========= 生成 x_new ========= %%
        near_to_rand = x_rand - x_near;
        normlized = near_to_rand / norm(near_to_rand) * Delta;
        x_new = x_near + normlized;
        %plot([x_near(1),x_new(1)],[x_near(2),x_new(2)],'color',[0.7,0.7,0.7], 'Linewidth', 0.5);
        %% ========= 碰撞检测 ========= %%
        if ~collisionChecking_deep1(x_near, x_new, Imp)
            continue;
        end

        %% ========= 1. 寻找 Near 集合 & 2. ChooseParent (Q-RRT* 版本) ========= %%
        % 收集所有在半径 50 内且与 x_new 无碰撞的节点索引（原逻辑保留）
        nearptr = [];
        nearcount = 0;
        near_dist = norm(x_new - x_near) + T.v(near_iter).dist;
        best_parent_idx = near_iter;   % 初始父节点为最近邻
        best_parent_cost = near_dist;

        % 对每个已有节点（除根节点外）进行检查
        for j = 2:N
            if j == near_iter
                continue;
            end
            x_j = [T.v(j).x, T.v(j).y];
            norm_dist = norm(x_new - x_j);
            if norm_dist < 50   % 半径阈值 (与原文一致)
                if collisionChecking_deep1(x_j, x_new, Imp)
                    nearcount = nearcount + 1;
                    nearptr(nearcount,1) = j;
                    % 先考虑该节点本身作为父节点
                    cost_via_j = T.v(j).dist + norm_dist;
                    if best_parent_cost > cost_via_j
                        best_parent_cost = cost_via_j;
                        best_parent_idx = j;
                    end
                    % 考虑该节点的祖先 (Depth 层)
                    ancList = getAncestors(T, j, Depth);
                    for k = 2:length(ancList)  % 跳过自身（已在上面考虑）
                        anc_idx = ancList(k);
                        anc_pos = [T.v(anc_idx).x, T.v(anc_idx).y];
                        cost_via_anc = T.v(anc_idx).dist + norm(x_new - anc_pos);
                        if best_parent_cost > cost_via_anc
                            % 检查路径无碰撞
                            if collisionChecking_deep1(anc_pos, x_new, Imp)
                                best_parent_cost = cost_via_anc;
                                best_parent_idx = anc_idx;
                            end
                        end
                    end
                end
            end
        end

        % 最终父节点
        x_parent = [T.v(best_parent_idx).x, T.v(best_parent_idx).y];

        %% ========= 将 x_new 加入树 ========= %%
        count = count + 1;
        T.v(count).x = x_new(1);
        T.v(count).y = x_new(2);
        T.v(count).xPrev = x_parent(1);
        T.v(count).yPrev = x_parent(2);
        T.v(count).dist = best_parent_cost;
        T.v(count).indPrev = best_parent_idx;

        %% ========= Rewire (Q-RRT* 版本) ========= %%
        % 获取 x_new 的祖先列表（包括自身）
        ancList_new = getAncestors(T, count, Depth);

        for k = 1:size(nearptr,1)
            idx_near = nearptr(k,1);
            x_near_pos = [T.v(idx_near).x, T.v(idx_near).y];
            current_cost = T.v(idx_near).dist;

            % 尝试用 x_new 及其每个祖先进行重连
            for a = 1:length(ancList_new)
                anc_idx = ancList_new(a);
                anc_pos = [T.v(anc_idx).x, T.v(anc_idx).y];
                new_cost = T.v(anc_idx).dist + norm(x_near_pos - anc_pos);
                if new_cost < current_cost
                    if collisionChecking_deep(anc_pos, x_near_pos, Imp)
                        % 更新父节点
                        T.v(idx_near).dist = new_cost;
                        T.v(idx_near).xPrev = anc_pos(1);
                        T.v(idx_near).yPrev = anc_pos(2);
                        T.v(idx_near).indPrev = anc_idx;
                        current_cost = new_cost;  % 更新为新的成本，以便后续祖先尝试
                        break;  % 找到更优即可跳出，因为祖先按距离递减，第一个有效的最优
                    end
                end
            end
        end

        %% ========= 检测是否到达目标 ========= %%
        if norm(x_new - goal) < 30 || collisionChecking_deep(x_new, goal, Imp)
            if (T.v(count).dist + norm(x_new - goal)) < start_goal_dist
                start_goal_dist = (T.v(count).dist + norm(x_new - goal));
                total_dis(iter,1) = start_goal_dist;
                tEnd = toc(Time);
                xxxx = [start_goal_dist, tEnd];
                data_ = [data_; xxxx];
                % 记录路径（用于绘图）
                path.pos = [];
                path.pos(1).x = x_G; path.pos(1).y = y_G;
                path.pos(2).x = T.v(end).x; path.pos(2).y = T.v(end).y;
                pathIndex = T.v(end).indPrev;
                j = 0;
                while 1
                    path.pos(j+3).x = T.v(pathIndex).x;
                    path.pos(j+3).y = T.v(pathIndex).y;
                    pathIndex = T.v(pathIndex).indPrev;
                    if pathIndex == 1 || pathIndex == 0
                        break
                    end
                    j = j + 1;
                end
                path.pos(end+1).x = x_I; path.pos(end).y = y_I;
            end
            if ~isempty(path.pos)
                    out_path=[];
                    for jj=1:length(path.pos)
                        ptt(1,1)=path.pos(jj).x;
                        ptt(1,2)=path.pos(jj).y;
                        out_path=[out_path;ptt];
                    end
                    s=flip(out_path);
                    return;
            end
            continue;
        end
    end
end

%% 辅助函数：获取节点 idx 的祖先列表（包括自身，向上 Depth 层）
function ancList = getAncestors(T, idx, Depth)
    ancList = idx;
    current = idx;
    for d = 1:Depth
        if T.v(current).indPrev ~= 0
            current = T.v(current).indPrev;
            ancList = [ancList, current];
        else
            break;
        end
    end
end

%% ================= 辅助函数 ================= %%
function feasible = collisionChecking_deep1(startPose, goalPose, map)
    x1 = round(startPose(1)); y1 = round(startPose(2));
    x2 = round(goalPose(1));  y2 = round(goalPose(2));
    [rows, cols] = size(map);
    if x1<1 || x1>cols || y1<1 || y1>rows || x2<1 || x2>cols || y2<1 || y2>rows
        feasible = false;
        return;
    end
    dx = abs(x2 - x1);
    dy = abs(y2 - y1);
    steps = max(dx, dy);
    if steps == 0
        feasible = (map(y1, x1) ~= 0);
        return;
    end
    x_points = round(linspace(x1, x2, steps+1));
    y_points = round(linspace(y1, y2, steps+1));
    idx = y_points + (x_points - 1) * rows;
    if any(map(idx) == 0)
        feasible = false;
    else
        feasible = true;
    end
end

% 用于 Informed 采样判断（原文已有，但未给出，此处补充）
function valid = new_node(x, y, c_best, source, goal)
    % 判断点 (x,y) 是否在 Informed 椭球体内
    % 椭球焦点为 source 和 goal，长轴为 c_best
    c_min = norm(source - goal);
    if c_best <= c_min
        valid = false;
        return;
    end
    center = (source + goal) / 2;
    % 将点变换到以 center 为原点，长轴为 x 轴的坐标系
    % 这里简化处理，仅用距离和判断（实际需旋转矩阵，但代码原样保留）
    dist_sum = norm([x,y] - source) + norm([x,y] - goal);
    if dist_sum <= c_best
        valid = true;
    else
        valid = false;
    end
end