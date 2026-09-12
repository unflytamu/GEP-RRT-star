
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Imp=im2bw(imread('ral_jincou.bmp')); 
source=[25 25];
goal=[384 690];
op_dis=759.2*1.01;
Delta= 10;    
mindis=15;
tic
%data=rrt_1(source,goal,op_dis,Imp,Delta);%与RRT*结合
data=informed_rrt_1(source,goal,op_dis,Imp);%与Informed RRT*结合
%data=Q_rrt_informed_1(source,goal,op_dis,Imp);%与Q_RRT*结合
init_path=getInitPath(data,mindis);
GA_path=GA_optmisation(init_path,op_dis,Imp,mindis,Delta);
disp([num2str(toc)]);
function s = GA_optmisation(init_path, op_dis, Imp, mindis, Delta)
    NP = 30;       % 种群数量
    max_gen = 30;  % 最大进化代数
    pc = 0.8;      % 交叉概率
    pm = 0.8;      % 变异概率
    
    new_pop1 = cell(NP, 1); 
    best_path_ever = init_path; 
    best_dis_ever = cal_dis1(init_path); 
    Time=tic;   
    for i = 1:NP
       if i == 1
           path = init_path;
       else
           path = getInitGAPath_new(init_path, Delta, Imp, 0.8);
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
    stagnation_count = 0; 
    best_dis_last = best_dis_ever; 

    for gen = 1 : max_gen
        path_value = cal_path_value1(new_pop1); 
        [min_dis, m] = min(path_value);         
        
        mean_path_value(gen) = mean(path_value);
        min_path_value(gen) = min_dis;
        if abs(min_dis - best_dis_last) < 0.5
            stagnation_count = stagnation_count + 1;
        else
            stagnation_count = 0;
        end
        best_dis_last = min_dis;
        if min_dis < best_dis_ever
            best_dis_ever = min_dis;
            best_path_ever = new_pop1{m, 1};
            if best_dis_ever < op_dis
                best_path_ever(:,3)=best_dis_ever;
                s = best_path_ever;
                return;
            end
        end 
        fit_value = path_value .^ -1;     
        new_pop2 = selection1(new_pop1, fit_value);
        new_pop2 = crossover1(new_pop2, pc); 
        if stagnation_count >=2 
            tmp_path = best_path_ever;
            [len, ~] = size(tmp_path);   
            if len > 20
                min_span = max(3, round(len * 0.10)); 
                max_span = max(8, round(len * 0.20)); 
                span = randi([min_span, max_span]);  
                max_start_idx = len - span - 1;
                if max_start_idx > 2
                    idx1 = randi([2, max_start_idx]);
                    idx2 = idx1 + span;                 
                    p1 = tmp_path(idx1, :);
                    p2 = tmp_path(idx2, :);                 
                    if collisionChecking_deep(p1, p2, Imp)
                        x_pts = round(linspace(p1(1), p2(1), idx2 - idx1 + 1));
                        y_pts = round(linspace(p1(2), p2(2), idx2 - idx1 + 1));
                        tmp_path(idx1:idx2, 1) = x_pts;
                        tmp_path(idx1:idx2, 2) = y_pts;
                        
                        new_len = cal_dis1(tmp_path);
                        if new_len < best_dis_ever
                            best_dis_ever = new_len;
                            best_path_ever = tmp_path;
                            if best_dis_ever<op_dis
                                best_path_ever(:,3)=best_dis_ever;
                                s = best_path_ever;
                                tEnd = toc(Time); 
                                return;
                            end
                            new_pop2{randi(NP), 1} = tmp_path;
                        end
                    end
                end
            end
        end
        new_pop2 = mutation1(new_pop2, pm, Imp, 3, m);%
        new_pop2{m, 1} = best_path_ever;
        new_pop1 = new_pop2;
    end
    best_path_ever(:,3)=cal_dis1(best_path_ever(:,1:2));
    s = best_path_ever; 
end

function path_value = cal_path_value1(pop)
    [px, ~] = size(pop);
    path_value = zeros(1, px);
    for i = 1:px
        path_value(1, i) = cal_dis1(pop{i, 1});
    end
end

function dis = cal_dis1(path)
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
    xL=size(Imp,2);
    yL=size(Imp,1);
    x_I=source(1,1); y_I=source(1,2);         
    x_G=goal(1,1); y_G=goal(1,2);     
    start_goal_dist = 1000000;
    path.pos(1).x = x_G;
    path.pos(1).y = y_G;
    knt=5000;
    total_dis=zeros(knt,1);
    T.v(1).x = x_I;        
    T.v(1).y = y_I; 
    T.v(1).xPrev = x_I;   
    T.v(1).yPrev = y_I;
    T.v(1).dist=0;        
    T.v(1).indPrev = 0;     
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
        x_new=[];
        near_to_rand = [x_rand(1)-x_near(1),x_rand(2)-x_near(2)];
        normlized = near_to_rand / norm(near_to_rand) * Delta;
        x_new = x_near + normlized;
        if ~collisionChecking_deep(x_near,x_new,Imp) 
           continue;
        end
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
        T.v(count).x = x_new(1);
        T.v(count).y = x_new(2); 
        T.v(count).xPrev = x_near(1);     
        T.v(count).yPrev = x_near(2);
        T.v(count).dist= norm(x_new - x_near) + T.v(near_iter).dist;          
        T.v(count).indPrev = near_iter;   
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
        if norm(x_new - goal) < 30 || collisionChecking_deep(x_new,goal,Imp)
            daad=T.v(count).dist + norm(x_new - goal);
            if (T.v(count).dist + norm(x_new - goal)) < start_goal_dist
                start_goal_dist = (T.v(count).dist + norm(x_new - goal));
                total_dis(iter)=start_goal_dist;
                tEnd = toc(Time); 
                xxxx=[start_goal_dist,tEnd];
                data_=[data_;xxxx];
                path.pos = [];
                if iter < knt
                    path.pos(1).x = x_G; path.pos(1).y = y_G;
                    path.pos(2).x = T.v(end).x; path.pos(2).y = T.v(end).y;
                    pathIndex = T.v(end).indPrev; 
                    j=0;
                    while 1
                        path.pos(j+3).x = T.v(pathIndex).x;
                        path.pos(j+3).y = T.v(pathIndex).y;
                        pathIndex = T.v(pathIndex).indPrev;
                        if pathIndex == 1
                            break
                        end
                        j=j+1;
                    end  
                    path.pos(end+1).x = x_I; path.pos(end).y = y_I; 
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
            continue;
        end
        if start_goal_dist<op_dis
            s=data_;
            break;
        end
    end
end

function feasible = collisionChecking_deep(startPose, goalPose, map)
    x1 = round(startPose(1)); y1 = round(startPose(2));
    x2 = round(goalPose(1)); y2 = round(goalPose(2));
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

function [new_pop] = crossover_ASAC_GPTT(pop, pc, d_th)
    [px, ~] = size(pop);
    new_pop =cell(px, 1); 
    for i = 1:2:px-1
        p1 = pop{i, 1};
        p2 = pop{i+1, 1};
        [n1, ~] = size(p1);
        [n2, ~] = size(p2);
        if n1 < 3 || n2 < 3
            return;
        end
        if rand < pc 
            [idx_in_p2, dists] = knnsearch(p2, p1);
            valid_idx = find(dists < d_th);
            if isempty(valid_idx)
                return;
            end  
            valid_idx = valid_idx(valid_idx > 1 & valid_idx < n1-1);
            if isempty(valid_idx)
                return;
            end 
            cross_idx_p1 = valid_idx(randi(length(valid_idx)));
            cross_idx_p2 = idx_in_p2(cross_idx_p1);
            if cross_idx_p2 <= 1 || cross_idx_p2 >= n2
                return;
            end
            new_pop{i, 1} = [p1(1:cross_idx_p1, :); p2(cross_idx_p2+1:end, :)];
            new_pop{i+1, 1} = [p2(1:cross_idx_p2, :); p1(cross_idx_p1+1:end, :)];
        end
        new_pop{i, 1} = p1;
        new_pop{i+1, 1} = p2;
    end
    if mod(px, 2) == 1
        new_pop{px, 1} = pop{px, 1};
    end
end

function [new_pop] = selection1(pop, fit_value)
    [px, ~] = size(pop);
    total_fit = sum(fit_value);
    if total_fit == 0
        fit_value = ones(px, 1); 
        total_fit = px;
    end
    
    edges = cumsum(fit_value)'; 
    
    r = rand(px, 1) * total_fit; 
   
    selected_indices = sum(r > edges, 2) + 1;
    new_pop = pop(selected_indices);
end

function [new_pop] = crossover1(pop, pc)
    [px, ~] = size(pop);
    new_pop = cell(px, 1);   
    for i = 1:2:px-1
        p1 = pop{i, 1};
        p2 = pop{i+1, 1};
        
        if rand < pc
            [~, idx_in_p1, idx_in_p2] = intersect(p1, p2, 'rows');
            m = length(idx_in_p1);
            if m >= 3
                r = randi([2, m-1]); 
                cut1 = idx_in_p1(r);
                cut2 = idx_in_p2(r);               
                new_pop{i, 1}   = [p1(1:cut1, :); p2(cut2+1:end, :)];
                new_pop{i+1, 1} = [p2(1:cut2, :); p1(cut1+1:end, :)];                
                continue; 
            end
        end
        new_pop{i, 1} = p1;
        new_pop{i+1, 1} = p2;
    end
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
            attempt = 0;           
            while attempt < 10
                rand_index = randi([1, max_start_idx]);
                p1 = single_new_pop(rand_index, :);
                p2 = single_new_pop(rand_index + knt, :);
                
                if collisionChecking_deep(p1, p2, G)
                    x_pts = round(linspace(p1(1), p2(1), knt + 1));
                    y_pts = round(linspace(p1(2), p2(2), knt + 1));
                    single_new_pop(rand_index+1 : rand_index+knt-1, 1) = x_pts(2:knt);
                    single_new_pop(rand_index+1 : rand_index+knt-1, 2) = y_pts(2:knt);
                    break; 
                end
                attempt = attempt + 1;
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
    xL=size(Imp,2);
    yL=size(Imp,1);
    count=1;
    start_goal_dist = 1000000;
    path.pos(1).x = (source(1,1)+goal(1,1))/2;
    path.pos(1).y = (source(1,2)+goal(1,2))/2;
    knt=5000;
    total_dis=zeros(knt,2);
    data_=[];
    Time=tic;
    for iter = 1:5000
        x_rand=[];
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
        x_new=[];
        near_to_rand = [x_rand(1)-x_near(1),x_rand(2)-x_near(2)];
        normlized = near_to_rand / norm(near_to_rand) * Delta;
        x_new = x_near + normlized;
        if ~collisionChecking_deep(x_near,x_new,Imp) 
           continue;
        end
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
        T.v(count).x = x_new(1);
        T.v(count).y = x_new(2); 
        T.v(count).xPrev = x_near(1);     
        T.v(count).yPrev = x_near(2);
        T.v(count).dist= norm(x_new - x_near) + T.v(near_iter).dist;          
        T.v(count).indPrev = near_iter;   
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
        if norm(x_new - goal) < 30 || collisionChecking_deep(x_new,goal,Imp)
            if (T.v(count).dist + norm(x_new - goal)) < start_goal_dist
                start_goal_dist = (T.v(count).dist + norm(x_new - goal));
                total_dis(iter,1)=start_goal_dist;
                tEnd = toc(Time); 
                xxxx=[start_goal_dist,tEnd];
                path.pos = [];
                if iter < 5000
                    path.pos(1).x = x_G; path.pos(1).y = y_G;
                    path.pos(2).x = T.v(end).x; path.pos(2).y = T.v(end).y;
                    pathIndex = T.v(end).indPrev;
                    j=0;
                    while 1
                        path.pos(j+3).x = T.v(pathIndex).x;
                        path.pos(j+3).y = T.v(pathIndex).y;
                        pathIndex = T.v(pathIndex).indPrev;
                        if pathIndex == 1
                            break
                        end
                        j=j+1;
                    end  
                    path.pos(end+1).x = x_I; path.pos(end).y = y_I;
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
    Thr = 10;            
    Delta = 10;         
    Depth = 2;                
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
        near_to_rand = x_rand - x_near;
        normlized = near_to_rand / norm(near_to_rand) * Delta;
        x_new = x_near + normlized;
        if ~collisionChecking_deep1(x_near, x_new, Imp)
            continue;
        end
        nearptr = [];
        nearcount = 0;
        near_dist = norm(x_new - x_near) + T.v(near_iter).dist;
        best_parent_idx = near_iter;   
        best_parent_cost = near_dist;
        for j = 2:N
            if j == near_iter
                continue;
            end
            x_j = [T.v(j).x, T.v(j).y];
            norm_dist = norm(x_new - x_j);
            if norm_dist < 50  
                if collisionChecking_deep1(x_j, x_new, Imp)
                    nearcount = nearcount + 1;
                    nearptr(nearcount,1) = j;
                    cost_via_j = T.v(j).dist + norm_dist;
                    if best_parent_cost > cost_via_j
                        best_parent_cost = cost_via_j;
                        best_parent_idx = j;
                    end
                    ancList = getAncestors(T, j, Depth);
                    for k = 2:length(ancList)  
                        anc_idx = ancList(k);
                        anc_pos = [T.v(anc_idx).x, T.v(anc_idx).y];
                        cost_via_anc = T.v(anc_idx).dist + norm(x_new - anc_pos);
                        if best_parent_cost > cost_via_anc
                            if collisionChecking_deep1(anc_pos, x_new, Imp)
                                best_parent_cost = cost_via_anc;
                                best_parent_idx = anc_idx;
                            end
                        end
                    end
                end
            end
        end
        x_parent = [T.v(best_parent_idx).x, T.v(best_parent_idx).y];
        count = count + 1;
        T.v(count).x = x_new(1);
        T.v(count).y = x_new(2);
        T.v(count).xPrev = x_parent(1);
        T.v(count).yPrev = x_parent(2);
        T.v(count).dist = best_parent_cost;
        T.v(count).indPrev = best_parent_idx;
        ancList_new = getAncestors(T, count, Depth);

        for k = 1:size(nearptr,1)
            idx_near = nearptr(k,1);
            x_near_pos = [T.v(idx_near).x, T.v(idx_near).y];
            current_cost = T.v(idx_near).dist;
            for a = 1:length(ancList_new)
                anc_idx = ancList_new(a);
                anc_pos = [T.v(anc_idx).x, T.v(anc_idx).y];
                new_cost = T.v(anc_idx).dist + norm(x_near_pos - anc_pos);
                if new_cost < current_cost
                    if collisionChecking_deep(anc_pos, x_near_pos, Imp)
                        T.v(idx_near).dist = new_cost;
                        T.v(idx_near).xPrev = anc_pos(1);
                        T.v(idx_near).yPrev = anc_pos(2);
                        T.v(idx_near).indPrev = anc_idx;
                        current_cost = new_cost;  
                        break; 
                    end
                end
            end
        end
        if norm(x_new - goal) < 30 || collisionChecking_deep(x_new, goal, Imp)
            if (T.v(count).dist + norm(x_new - goal)) < start_goal_dist
                start_goal_dist = (T.v(count).dist + norm(x_new - goal));
                total_dis(iter,1) = start_goal_dist;
                tEnd = toc(Time);
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
function valid = new_node(x, y, c_best, source, goal)
    c_min = norm(source - goal);
    if c_best <= c_min
        valid = false;
        return;
    end
    center = (source + goal) / 2;
    dist_sum = norm([x,y] - source) + norm([x,y] - goal);
    if dist_sum <= c_best
        valid = true;
    else
        valid = false;
    end
end

function s = getInitGAPath_new(path, D, obs, p)
    [n, ~] = size(path);
    out = path(1, :);
    last_shift = 0; 
    for i = 2:n-1
        if rand < p 
            pt = path(i, :);
            is_ok = false;
            dx = path(i+1, 1) - path(i, 1);
            dy = path(i+1, 2) - path(i, 2);
            len = sqrt(dx^2 + dy^2);
            if len < 1e-6
                out = [out; pt];
                continue;
            end
            normal_vec_x = -dy / len;
            normal_vec_y = dx / len;
            for attempt = 1:10
                target_shift = (rand * 2 - 1) * D;
                current_shift = 0.7 * last_shift + 0.3 * target_shift;
                
                pt_new_x = round(pt(1,1) + current_shift * normal_vec_x);
                pt_new_y = round(pt(1,2) + current_shift * normal_vec_y);
                if collisionChecking_deep([pt_new_x, pt_new_y], [pt_new_x, pt_new_y], obs)       
                    out = [out; [pt_new_x, pt_new_y, path(i,3)]];
                    last_shift = current_shift; 
                    is_ok = true;
                    break;
                end
            end
            if ~is_ok
                out = [out; path(i, :)];
                last_shift = 0;
            end
        else
            out = [out; path(i, :)];
            last_shift = 0; 
        end
    end
    s = [out; path(n, :)];
end
