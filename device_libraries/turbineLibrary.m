% https://mainrenewableenergy.com/2018-wind-turbine-cost/ (lots of info
% sourced from this blog)
TurbPPI = 1.3; %2018->2024 PPI
% Nature Power 400 W turbine (home depot)
turbineLib(1).kW = .4;
turbineLib(1).cost = 431.54*TurbPPI;

% Coleman 400 W turbine (home depot)
turbineLib(2).kW = .4;
turbineLib(2).cost = 454.55*TurbPPI;

% Nature Power 2000 W turbine (home depot)
turbineLib(3).kW = 2;
turbineLib(3).cost = 2095.65*TurbPPI;

% Ramsond 1000 W turbine (home depot)
turbineLib(4).kW = 1;
turbineLib(4).cost = 909.99*TurbPPI;

% Ramsond 3000 W turbine (home depot)
turbineLib(5).kW = 3;
turbineLib(5).cost = 4356*TurbPPI;

% Nature's Generator 400 W turbine (home depot)
turbineLib(6).kW = .4;
turbineLib(6).cost = 549*TurbPPI;

% ISTA Breeze Heli 4.0 (ebay)
turbineLib(7).kW = 3.8;
turbineLib(7).cost = 4796.98*TurbPPI;

% ISTA Breeze Heli 2.0 (ebay)
turbineLib(8).kW = 2;
turbineLib(8).cost = 1640*TurbPPI;

% Pika Energy 1.7 kW (dwr + google + mainrenewableenergy)
turbineLib(9).kW = 1.7;
%turbine(9).cost = 5997;
turbineLib(9).cost = 6675*TurbPPI;

% Bergey Excel 10 (dwr + bergey + mainrenewableenergy)
turbineLib(10).kW = 10;
turbineLib(10).cost = 31770*TurbPPI; 

% Eveready Diversified Kestrel e400nb (dwr + 
turbineLib(11).kW = 3.5;
turbineLib(11).cost = nan;

% Kingspan Environmental KW6 (dwr + 
turbineLib(12).kW = 6;
turbineLib(12).cost = nan;

% Kingspan Environmental KW3ex (kingspan + 
turbineLib(13).kW = 2.5;
turbineLib(13).cost = nan;

% Kingspan Environmental KW3 (kingspan + 
turbineLib(14).kW = 2.5;
turbineLib(14).cost = nan;

% Lely Aircon BV LA (dwr + lely +
turbineLib(15).kW = 9.6;
turbineLib(15).cost = nan;

% Osiris 10 (dwr + mainrenewableenergy)
turbineLib(16).kW = 10;
turbineLib(16).cost = 27500*TurbPPI;

% Bergey Excel R (mainrenewableenergy + bergey)
turbineLib(17).kW = 7.5;
turbineLib(17).cost = 26870*TurbPPI;

% Bergey Excel 6 (mainrenewableenergy + bergey)
turbineLib(18).kW = 5.5;
turbineLib(18).cost = 21995*TurbPPI;

% Weaver 5 (mainrenewableenergy + weaver + homepower)
turbineLib(19).kW = 5.03;
turbineLib(19).cost = 38576*TurbPPI;
turbineLib(19).cost = nan;

% Weaver 2 (weaver + mainrenewableenergy)
turbineLib(20).kW = 2;
turbineLib(20).cost = nan;

% Ventera VT-10 (mainrenewableenergy + wingenpower)
turbineLib(21).kW = 10;
turbineLib(21).cost = 24000*TurbPPI;

% Luminous Whisper 500 (wish energy + mainrenewableenergy)
turbineLib(22).kW = 3.2;
turbineLib(22).cost = 9466*TurbPPI;

% Xzeres Skystream 3.7 (dwr + xzeres + mainrenewableenergy)
turbineLib(23).kW = 2.4;
turbineLib(23).cost = 3420*TurbPPI;

% Luminous Whisper 200 (southwest + mainrenewableenergy)
turbineLib(24).kW = 1;
turbineLib(24).cost = 3921*TurbPPI;

% Bergey Excel 1 (bergey + mainrenewableenergy)
turbineLib(25).kW = 1;
turbineLib(25).cost = 4595*TurbPPI;

% Sumec PWB01 30 48 (dwr + interek)
turbineLib(26).kW = 1.17;
turbineLib(26).cost = nan;

% Sumec PWB02 20 48 (dwr + interek)
turbineLib(27).kW = 1.7;
turbineLib(27).cost = nan;

% Sumec PWB03 44 250 (dwr + interek)
turbineLib(28).kW = 3.24;
turbineLib(28).cost = nan;

% Sumec PWB05 280 (dwr + interek)
turbineLib(29).kW = 5.01;
turbineLib(29).cost = nan;

%%Potential add on turbines for low power cost values
% %Shine Turbine (https://shineturbine.com/products/shinewindturbine)
% turbineLib(30) = 40/1000;
% turbineLib(30).cost = 399.99;
% 
% %ATO (https://www.ato.com/100w-wind-turbine?srsltid=AfmBOoo9fZn2Ohukte8TzJJ7BHwjpTHXr07ZeBZ99sZQMtD2wt9bx6-P)
% turbineLib(30) = 100/1000;
% turbineLib(30).cost = 293.69;

%outliers for linear regression
% turbine(26).kW = 9;
% turbine(26).cost = 5000;
% 
% turbine(27).kW = 9.5;
% turbine(27).cost = 40000;











