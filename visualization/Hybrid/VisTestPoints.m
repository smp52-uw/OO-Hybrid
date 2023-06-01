%visualize points tested

output = argBasin_6_0_p2t_v11.output;
%output = argBasin_5_5_p2t_v9.output;
% %IT 1
% output = argBasin_6_0_p2t_v10_3it.output;
% figure(1)
% scatter3(output.Kwi_run{1,1},output.Ki_run{1,1},output.Kwa_run{1,1},'.')
% xlabel('Wind')
% ylabel('Solar')
% zlabel('Wave')
% title('IT 1')
% surv1 = output.surv{1,1};
% surv1(surv1< 0.96) = 0;
% surv1(surv1>0.993) = 0;
% Kd_old1 = output.Kd_run{1,1}(surv1 > 0);
% Ki_old1 = output.Ki_run{1,1}(surv1 > 0);
% Kwi_old1 = output.Kwi_run{1,1}(surv1 > 0);
% Kwa_old1 = output.Kwa_run{1,1}(surv1 > 0);
% Kc_old1 = output.Kc_run{1,1}(surv1 > 0);
% S_old1 = output.S_run{1,1}(surv1>0);
% 
% %IT 2 points
% Kd2 = output.Kd_run{1,2};
% Ki2 = output.Ki_run{1,2};
% Kwi2 = output.Kwi_run{1,2};
% Kwa2 = output.Kwa_run{1,2};
% Kc2 = output.Kc_run{1,2};
% S2 = output.S_run{1,2};
% 
% old_pts_test = zeros(1,length(Kd_old1));
% for i=1:length(Kd_old1)
%     old_pts_test(i) = find(Kd_old1(i) == Kd2 & Ki_old1(i) == Ki2 & Kwi_old1(i) == Kwi2 & Kwa_old1(i) == Kwa2 & Kc_old1(i) == Kc2 & S_old1(i) == S2);
% end
% 
% figure(2)
% scatter3(output.Kwi_run{1,2}(old_pts_test),output.Ki_run{1,2}(old_pts_test),output.Kwa_run{1,2}(old_pts_test),'r.')
% hold on
% scatter3(Kwi_old1,Ki_old1,Kwa_old1,'gx')
% xlabel('Wind')
% ylabel('Solar')
% zlabel('Wave')
% title('IT 2')
% 
% surv2 = output.surv{1,2};
% % surv2(surv2< 0.96) = 0;
% % surv2(surv2>0.993) = 0;
% Kd_old2 = output.Kd_run{1,2}(surv2 > 0);
% Ki_old2 = output.Ki_run{1,2}(surv2 > 0);
% Kwi_old2 = output.Kwi_run{1,2}(surv2 > 0);
% Kwa_old2 = output.Kwa_run{1,2}(surv2 > 0);
% Kc_old2 = output.Kc_run{1,2}(surv2 > 0);
% S_old2 = output.S_run{1,2}(surv2 >0);
% 
% %IT 3 points
% Kd3 = output.Kd_run{1,3};
% Ki3 = output.Ki_run{1,3};
% Kwi3 = output.Kwi_run{1,3};
% Kwa3 = output.Kwa_run{1,3};
% Kc3 = output.Kc_run{1,3};
% S3 = output.S_run{1,3};
% 
% % old_pts_test2 = zeros(1,length(Kd_old2));
% % for i=1:length(Kd_old2)
% %     old_pts_test2(i) = find(Kd_old2(i) == Kd3 & Ki_old2(i) == Ki3 & Kwi_old2(i) == Kwi3 & Kwa_old2(i) == Kwa3 & Kc_old2(i) == Kc3 & S_old2(i) == S3);
% % end
% 
% old_min = find(output.min.kWd{1,2} == Kd3 & output.min.kWi{1,2} == Ki3 & output.min.kWwi{1,2} == Kwi3 & output.min.kWwa{1,2} == Kwa3 & output.min.kWc{1,2} == Kc3 & output.min.Smax{1,2} == S3);
% 
% old_min_it2 = find(output.min.kWd{1,2} == Kd_old2 & output.min.kWi{1,2} == Ki_old2 & output.min.kWwi{1,2} == Kwi_old2 & output.min.kWwa{1,2} == Kwa_old2 & output.min.kWc{1,2} == Kc_old2 & output.min.Smax{1,2} == S_old2);

% id_wa3 = output.Kwa_run{1,3} == 0;
% cost = output.cost{1,3};
% surv = output.surv{1,3};
% min_cost = min(output.cost{1,3}(surv>=0.99));
% id_c3 = (cost< (min_cost*0.8) | cost > (min_cost*1.2));
% id_s3 = surv>0;
% id_t = id_wa3' & id_c3 & id_s3;
idi = output.Kwa_run{1,3} == 0;
figure(1)
scatter3(output.Kwi_run{1,3}(idi),output.Ki_run{1,3}(idi),output.Kwa_run{1,3}(idi),'.')
xlabel('Wind')
ylabel('Solar')
zlabel('Wave')
title('IT 3')
% 
% id_wa = output.Kwa_run{1,4} == 0;

figure(2)
idi4 = output.Kwa_run{1,4} == 0;
% scatter3(output.Kwi_run{1,4}(idi4),output.Ki_run{1,4}(idi4),output.Kwa_run{1,4}(idi4),'.')
scatter3(output.Kwi_run{1,4},output.Ki_run{1,4},output.Kwa_run{1,4},'.')
xlabel('Wind')
ylabel('Solar')
zlabel('Wave')
title('IT 4')