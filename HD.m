% File HD.m containing variables for the Hydrodynamic Data
% The data in this file is from the software Wamit, provided by Dominique
% Roddier

function [T,F,a,b] = HD

% Wave period, initially just test for one period
T = [4.24E+00;
4.95E+00;
5.66E+00;
6.36E+00;
7.07E+00;
7.78E+00;
8.49E+00;
9.19E+00;
9.55E+00;
9.90E+00;
1.03E+01;
1.06E+01;
1.13E+01;
1.20E+01;
1.27E+01;
1.34E+01;
1.41E+01;
1.56E+01;
1.70E+01;
1.84E+01;
2.12E+01;
2.55E+01;
3.11E+01;
3.68E+01;
4.24E+01];

% Wave exiting force (at wave amplitude 1) for the given wave period
F = [56046.8125;
76934.7285;
96051.515;
114007.9815;
131638.402;
149441.6455;
167436.7485;
185264.4675;
193953.564;
202418.608;
210598.472;
218460.2035;
233139.452;
246372.245;
258191.29;
268696.4;
278013.75;
293628.825;
306003.53;
315913.535;
330514.31;
344361.465;
355004.51;
361203.255;
365124.48];

% Corresponding added mass for the given wave period
a = [42596.9125;
43226.9125;
44068.4375;
44900.5;
45678.8375;
46427.975;
47166.025;
47882.05;
48222.05;
48547.525;
48858;
49148.5375;
49678.4875;
50146.9125;
50556.4625;
50922.65;
51250.7375;
51819.1;
52303.95;
52731.725;
53457.9875;
54346.475;
55299;
56083.1125;
56752.85];

% Corresponding Hydrodynamic damping for the given wave period
b = [5414.696938;
6426.680946;
6710.022235;
6635.40741;
6429.533594;
6176.755365;
5892.03195;
5578.435971;
5414.210689;
5247.947607;
5081.908153;
4918.026313;
4603.036372;
4311.327214;
4045.805829;
3806.430019;
3591.198384;
3223.388817;
2922.971759;
2674.178076;
2286.788872;
1881.529145;
1524.553823;
1282.811533;
1107.900985];