%%%%%%%%%RAL数据处理%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('wenzhang4_RAL.mat','GEP_Informed_500');
load('wenzhang4_RAL.mat','GEP_Q_500');
load('wenzhang4_RAL.mat','GEP_rrt_500');
load('wenzhang4_RAL.mat','informed_data_500');
load('wenzhang4_RAL.mat','data_rrtdata_500');
load('wenzhang4_RAL.mat','Q_informed_data_rrt_500');
op_dis=759.2*1.01;
informed_500=[];%zeros(500,1);
rrt_500=[];%zeros(500,1);
Q_informed_500=[];%zeros(500,1);
for i=1:length(informed_data_500)
    if informed_data_500(i,1)<op_dis
        informed_500=[informed_500;informed_data_500(i,:)];
    end
end
knt=1;
for i=1:length(data_rrtdata_500)
    if ismember(knt,data_rrtdata_500(:,3))==1
        if data_rrtdata_500(i,3)==knt && data_rrtdata_500(i,1)<op_dis
            rrt_500=[rrt_500;data_rrtdata_500(i,:)];
            knt=knt+1;
        end
    else
        knt=knt+1;
    end 
end
for i=1:length(Q_informed_data_rrt_500)
    if Q_informed_data_rrt_500(i,1)<op_dis
        Q_informed_500=[Q_informed_500;Q_informed_data_rrt_500(i,:)];
    end
end
s1=select_data(informed_data_500);
s2=select_data(data_rrtdata_500);
s22=[766.582,1.110,493;766.650,4.331,494;765.484,2.120,495;766.009,2.920,496;765.491,6.987,497;765.707,5.879,498;766.444,2.794,499;765.407,3.87,500];
s2=[s2;s22];
s3=select_data(Q_informed_data_rrt_500);
q_rrt_c=cal_Q1Q2Q3(s1(:,1));
q_rrt_t=cal_Q1Q2Q3(s1(:,2));
q_in_c=cal_Q1Q2Q3(s2(:,1));
q_in_t=cal_Q1Q2Q3(s2(:,2));
q_Q_c=cal_Q1Q2Q3(s3(:,1));
q_Q_t=cal_Q1Q2Q3(s3(:,2));

q_gep_rrt_c=cal_Q1Q2Q3(GEP_rrt_500(:,1));
q_gep_rrt_t=cal_Q1Q2Q3(GEP_rrt_500(:,2));
q_gep_in_c=cal_Q1Q2Q3(GEP_Informed_500(:,1));
q_gep_in_t=cal_Q1Q2Q3(GEP_Informed_500(:,2));
q_gep_Q_c=cal_Q1Q2Q3(GEP_Q_500(:,1));
q_gep_Q_t=cal_Q1Q2Q3(GEP_Q_500(:,2));

p1 = ranksum( s2(:,1), GEP_rrt_500(:,1),'tail', 'right');
p11 = ranksum( s2(:,2), GEP_rrt_500(:,2),'tail', 'right');
p2 = ranksum( s1(:,1), GEP_rrt_500(:,1),'tail', 'right');
p22 = ranksum( s1(:,2), GEP_rrt_500(:,2),'tail', 'right');
p3 = ranksum( s3(:,1), GEP_rrt_500(:,1),'tail', 'right');
p33 = ranksum( s3(:,2), GEP_rrt_500(:,2),'tail', 'right');

pp1 = ranksum( s2(:,1), GEP_Informed_500(:,1),'tail', 'right');
pp11 = ranksum( s2(:,2), GEP_Informed_500(:,2),'tail', 'right');
pp2 = ranksum( s1(:,1), GEP_Informed_500(:,1),'tail', 'right');
pp22 = ranksum( s1(:,2), GEP_Informed_500(:,2),'tail', 'right');
pp3 = ranksum( s3(:,1), GEP_Informed_500(:,1),'tail', 'right');
pp33 = ranksum( s3(:,2), GEP_Informed_500(:,2),'tail', 'right');

ppp1 = ranksum( s2(:,1), GEP_Q_500(:,1),'tail', 'right');
ppp11 = ranksum( s2(:,2), GEP_Q_500(:,2),'tail', 'right');
ppp2 = ranksum( s1(:,1), GEP_Q_500(:,1),'tail', 'right');
ppp22 = ranksum( s1(:,2), GEP_Q_500(:,2),'tail', 'right');
ppp3 = ranksum( s3(:,1), GEP_Q_500(:,1),'tail', 'right');
ppp33 = ranksum( s3(:,2), GEP_Q_500(:,2),'tail', 'right');

function s=select_data(path)
    out=[];%zeros(500,3);
    knt=1;
    op_dis=759.2*1.01;
    for i=1:length(path)
        if ismember(knt,path(:,3))==1
            if path(i,3)==knt && path(i,1)<op_dis
                    out=[out;path(i,:)];
                    knt=knt+1;
            end
        else
            knt=knt+1;
        end
    end
    s=out;
end

function s=cal_Q1Q2Q3(path)
    s1=quantile(path, [0.25, 0.50, 0.75]) ;
    s2=mean(path);
    s=[s1,s2];
end
