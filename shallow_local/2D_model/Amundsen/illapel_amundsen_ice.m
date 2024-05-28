%%% modified code from Nurbek Tazhimbetov's PhD thesis, 2022 %%%
%%% changes made by Nestor Walters, advised by Eric Dunham %%%

% illapel_amundsen_ice.m is the culminating code that runs the simulation
% it calls multiple other files, the most important being Fullplate.m

%%% Some notes:
%%% run experiments by adjusting bathymetry (bathy), ice thickness (thick),
%%% and bending stiffness (BendS) of the ice, instructions inline

%%% save_matrix = 'yes' and running the whole simulation outputs
% mat_ABq_bvar_tvar.mat, which has stored the A and B matrices and q0
% from Nurbek's original Adq/dt = Bq + forcing problem
% the filename for matrices can be changed in Fullplate.m

%%% note that every change to bathy, thick, or BendS will yield different
% matrices, but unless you change the matname variable in Fullplate.m
% they will overwrite the previously stored matrices

%%% when calling illapel_amundsen_ice.m in get_eigs, get_latlong, or
% plot_evecs, you can comment out everything after q = Fullplate([] ...
% to greatly reduce run time



% matrix name 

clear
clc
format long
close all

addpath('../')
addpath('../FullPlate')
addpath('../../sbplib')

load('BED.mat')
load('SURFACE.mat')
load('THICKNESS.mat')
load('MASK.mat')
load('points.mat')

cond_less = 1e3;

gravity = 9.8 / cond_less;
rho_ice = 920 / cond_less;
rho_water = 1000 / cond_less;
time_width = 30;
in_angle = 7 * pi / 4;
%time_in_sec = 60; %% added: shorten when only want to extract matrices A,B
time_in_sec = 43200; %% set this time for original simulation
x_base = 0;
y_base = 0;

boundary_type = 'mixed';
anime = 'no';
video_mode = 'no';
symmetry_check = 'no';
time_solver = 'implicit';
MMS = 'no';
save_matrix = 'no'; %% for the A,B matrices in Adq/dt = Bq + forcing
non_reflective = 'yes';
boundary_data = 'tsunami_data';
tsunami_data_construct.file = 'tsunami_data_amundsen.m'; %%edited
tsunami_data_construct.interpolated = 'illapel_amundsen_S_nr.mat';%%edited
save_in_binary = 'no'; %% edit to not save simulation
q_file_name = 'illapel_amundsen_ice_q15_binary.bin'; %% set to 'midrange'

order = 4;
m = 40; % 12
iLU_tol = 0.00000001; % 0.000001;
plot_details.left = 0;
plot_details.right = 275;
plot_details.top = -13;
plot_details.bottom = -220;
plot_details.caxis_low = -.15;
plot_details.caxis_up = .15;

% k_max = 1 / (3 * (1 / time_width));
k_max = 16; % 14.4
[k, iter] = alignedTimestep(k_max, time_in_sec);


BATHY = (SURFACE - BED - THICKNESS);
[Bn, Bm] = size(BATHY);
x_axis = linspace(0, Bm/2, Bm);
y_axis = linspace(-Bn/2, 0, Bn);
% % figure()
% % imagesc('XData', x_axis, 'YData', y_axis, 'CData', flip(BATHY))
% % colorbar
BATHY(BATHY <= 30) = 30;
BATHY = BATHY / cond_less;
[X, Y] = meshgrid(linspace(0, Bm/2, Bm), linspace(-Bn/2, 0, Bn));
%%% (modified) 1 of 2 experiments on constancy effects on eigenmodes
%%% to make bathymetry constant, comment out bathy = @(x,y)..."
%%% and uncomment conBath, bathy = @(x,y)*conBath
bathy = @(x, y) my_interpolation(x, y, X, Y, flip(BATHY));
% conBath = 0.6;
% bathy = @(x, y) 0.6;

nu_const = 0.3;
depth_const = 600 / cond_less; 
thick_const = 300 / cond_less;
bend_const = 27 * 1e15 / cond_less^6;

nu = @(x, y) nu_const * (1 + 0 * sin( - 1e-1 * x - 3 * 1e-1 * y) / 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run amundsen_mesh_setup.m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESTIMATING THE ICE THICKNESS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

THICKNESS(isnan(THICKNESS)) = 500;
THICKNESS = THICKNESS / cond_less;

%%% (modified) experiment 2 of 2 on constancy effects on eigenmodes
thick = @(x, y) interp2(X, Y, flip(THICKNESS), full(x), full(y), 'spline');
F = grid.evalOn(g_ice, thick);
%%% to make ice thickness constant, comment out "thick = @(x,y)..."
%%% and uncomment conThick, "F = 0*F + conThick"
% conThick = 0.6;
% F = 0*F + conThick;
F(isnan(F)) = .25;
F(F <= .25) = .25;
full_F = [F; zeros(n_water - n_ice, 1)]; % ONLY ICE

BendS = (full_F * cond_less / 100).^3 * 10^15 / cond_less^6;
NU = grid.evalOn(g_ice, nu);
full_NU = [NU; zeros(n_water - n_ice, 1)];


M_s = m;
M_r = g.grids{2}.logic.m(1);

h_x = 90 / (3 * M_s - 1);
lambda = 8 * 2 * 25 / 12;
wvn_K  = 2 * pi / lambda;
wvn_x  = wvn_K * cos(pi/4);
wvn_y  = wvn_K * sin(pi/4);
omega = wvn_K * sqrt(depth_const * (bend_const * wvn_K^4 + rho_water * gravity) / (rho_water + depth_const * wvn_K^2 * rho_ice * thick_const));
ppwv   = lambda / h_x;
period = 2 * pi / omega;

switch order
	case 2
        C_const = 0.09;
	case 4
        C_const = 0.03;
	case 6
        C_const = 0.003;
end

bc_ice_1_d3.boundary = def_ice.boundaryGroups.FR;
bc_ice_1_d3.type     = 'd3';
bc_ice_1_d3.data     = [];

bc_ice_1_d2.boundary  = def_ice.boundaryGroups.FR;
bc_ice_1_d2.type      = 'd2_modified';
bc_ice_1_d2.data      = [];

bc_ice_1_d1.boundary = def_ice.boundaryGroups.CL;
bc_ice_1_d1.type     = {'d1', def_ice.boundaryGroups.FR};
bc_ice_1_d1.data     = [];

bc_ice_1_e.boundary  = def_ice.boundaryGroups.CL;
bc_ice_1_e.type      = {'e', def_ice.boundaryGroups.FR};
bc_ice_1_e.data      = [];

bc_ice_2_free.boundary = def_ice.boundaryGroups.FR;
bc_ice_2_free.type = 'free_boundaries';
bc_ice_2_free.data = [];
 
bc_ice_2_clamped.boundary = def_ice.boundaryGroups.CL;
bc_ice_2_clamped.type     = {'clamped_boundaries', def_ice.boundaryGroups.CL, def_ice.boundaryGroups.FR};
bc_ice_2_clamped.data     = [];

bc_ice1 = {bc_ice_1_d3, bc_ice_1_d2, bc_ice_1_d1, bc_ice_1_e};
bc_ice2 = {bc_ice_2_free, bc_ice_2_clamped};

bc_water.boundary = def.boundaryGroups.all;
bc_water.type = 'N';
bc_water.data = [];
bc_water = {bc_water};

bc_water_NR.boundary = def.boundaryGroups.NR;
bc_water_NR.type     = 'N';
bc_water_NR.data     = [];
bc_water_NR = {bc_water_NR};

A = 0.01;
B = A * omega / wvn_K^2 / depth_const;

w0 = grid.evalOn(g, @(x, y) 0 * cos(wvn_x * x + wvn_y * y));
% w0 = grid.evalOn(g, @(x, y) exp(-.5 * ((x - 60).^2 / 2.5^2 + (y + 30).^2 / 2.5^2)));
phi0 = grid.evalOn(g, @(x, y) -0 * sin(wvn_x * x + wvn_y * y));

% clearvars -except bc_ice1 bc_ice2 bc_water bc_water_NR w0 phi0 order M_s k iter g rho_ice rho_water gravity full_NU BendS bathy full_F anime symmetry_check time_solver MMS save_matrix non_reflective g_ice time_width in_angle video_mode boundary_data tsunami_data_construct


% q = FullPlate([], [], [], ...
%            bc_ice1, bc_ice2, bc_water, [], bc_water_NR, ...
%            w0, phi0, order, M_s, k, iter, ...
%            g, rho_ice, rho_water, gravity, ...
%            full_NU, BendS, bathy, full_F, ...
%            anime, [], [], ...
%            symmetry_check, time_solver, ...
%            MMS, save_matrix, non_reflective, ...
%            g_ice, time_width, in_angle, video_mode, ...
%            boundary_data, tsunami_data_construct, ...
%            water_block_N, x_base, y_base, iLU_tol, ...
%            plot_details, q_file_name);
% 




