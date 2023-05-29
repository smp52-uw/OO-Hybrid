%visualize points tested

output = argBasin_5_5_p2t_v9.output;

figure(1)
scatter3(output.Kwi_run{1,3},output.Ki_run{1,3},output.Kwa_run{1,3},'.')
xlabel('Wind')
ylabel('Solar')
zlabel('Wave')
title('IT 3')

figure(2)
scatter3(output.Kwi_run{1,4},output.Ki_run{1,4},output.Kwa_run{1,4},'.')
xlabel('Wind')
ylabel('Solar')
zlabel('Wave')
title('IT 4')

figure(3)
scatter3(output.Kwi_run{1,5},output.Ki_run{1,5},output.Kwa_run{1,5},'.')
xlabel('Wind')
ylabel('Solar')
zlabel('Wave')
title('IT 5')