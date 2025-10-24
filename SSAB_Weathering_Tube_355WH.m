close all; clear;

square_row = [40; 50; 60; 70; 80; 90; 100; 120; 140; 150; 160; 180; 200; 220; 250; 300];
heavy_dim = [4.0; 5.0; 5.0; 5.0; 6.0; 6.0; 8.0; 10; 12.5; 12.5; 12.5; 12.5; 12.5;12.5; 12.5; 12.5];
light_dim = [2.0; 2.0; 2.0; 3.0; 3.0; 3.0; 3.0; 4.0;4.0;4.0;4.0;5.0; 5.0; 6.0;6.0;6.0]

model1h=fit(log(square_row),log(heavy_dim),'poly1')

model2h=fit(log(square_row(1:11)),log(heavy_dim(1:11)),'poly1')

model1l=fit(log(square_row),log(light_dim),'poly1')

model2l=fit(log(square_row(2:15)),log(light_dim(2:15)),'poly1')

area1l=4*exp(model1l.p2)*(square_row.^(model1l.p1+1))
weight1l=7800e-6*area1l

grayColor = [.7 .7 .7];

figure
hold on
patch([square_row' fliplr(square_row')], [light_dim' fliplr(heavy_dim')], [.8 .8 .8])
%plot(square_row,heavy_dim,'-r')
%plot(square_row,light_dim,'-k')
plot(square_row(1:11),exp(model2h.p2)*(square_row(1:11).^model2h.p1),':k')
plot(square_row,exp(model1l.p2)*(square_row.^model1l.p1),'--k')

xlabel('section width, mm')
ylabel('wall thickness, mm')

legend({'available sections','fit heavy t ~ h^0^.^8^9^7','fit light t ~ h^0^.^6^7^7'},...
    'Location','southeast') %'heavy section limit', 'light section limit',

hold off



