%% TASK 1
% Pau Cornudella Quer i Roger Jordi Martinez Pardell
clear; close all; clc;
load("model_data_clean.mat");


%% Aparatat 1
% 0)
K = KAAX; % Stiffness matrix
M = MAAX; % Mass matrix
n_nodes = size(node_coords,1); % N nodes
ndof    = 6*n_nodes; % N DOFs 
in_glo  = (1:ndof)'; % Global DOF index vector
support_nodes = case_control_sets.subcase_0_SET_2;   % 3 fixacions

% 1)
fixnodes = zeros(3*length(support_nodes),3);
k = 1;
for i = 1:length(support_nodes)
    for j = 1:3   % DOF 1,2,3 (ux, uy, uz)
        fixnodes(k,:) = [support_nodes(i) j 0];
        k = k + 1;
    end
end

% 2,3 i 4) index global dels DOF = dof + (node-1)*6
in_d = (fixnodes(:,1)-1)*6 + fixnodes(:,2);   % indexs globals dels DOF fixats
u_d  = fixnodes(:,3); % valors prescrits (tots 0)
in_n = setdiff(in_glo, in_d); % Free DOFs

% 5) Gravity F = M*g
g = [0; 0; -9.81e3; 0; 0; 0]; % Unitats! g = 9.81e3 mm/s^2
g_glo = repmat(g,n_nodes,1); %Replica l'array g a tots els nodes
F = M * g_glo;


% Theory_statics eq 6
% K_NN * u_N = F_N - K_ND*u_D
u_D = u_d;
K_ND = K(in_n, in_d);
K_NN = K(in_n, in_n);
K_DN = K(in_d, in_n);
K_DD = K(in_d, in_d);
F_N  = F(in_n);
u_N = K_NN \ (F_N-K_ND*u_D);
F_D = K_DD * u_D + K_DN * u_N; % ç
r_f = F_D - F(in_d); %F(in_d) = F externes

%Gravity sag
u_full = zeros(ndof, 1);
u_full(in_n) = u_N;
uz = u_full(3:6:end); % mm
uz0 = uz - mean(uz);
gravity_sag_nm_rms = sqrt(mean(uz0.^2)) * 1e6; % nm-RMS

%Comparació massa amb reaccions:
Fz = F(3:6:end);
massa_total = -sum(Fz) / 9.81e3;
r_z = r_f(3:3:end);
% Comparació (hauria de donar 0)
sum(r_z) + sum(Fz);

%Reaccions als suports
R_support_1 = r_f(1:3);
R_support_2 = r_f(4:6);
R_support_3 = r_f(7:9);

%Exportació al META
uz = u_full(3:6:end);
uz0 = uz - mean(uz);
u_grav = zeros(n_nodes,6);
u_grav(:,3) = uz0;
fillhdf('h5template.h5','META_1.h5',u_grav*1e10);
% Comprovació
D = h5read('META_1.h5','/NASTRAN/RESULT/NODAL/DISPLACEMENT');
max(abs(D.Z))

%% Apartat 2
n_modes = 100; % enunciat
Z = zeros(n_nodes,n_modes);

x = node_coords(:,1);
y = node_coords(:,2);
r = sqrt(x.^2 + y.^2);
R = max(r); % correcte = 75
x_norm = x/R; %normalized coordinates
y_norm = y/R; %normalized coordinates
R_norm = max(sqrt(x_norm.^2 + y_norm.^2));
    
for j = 1:n_modes
    [n, m, name, val, Z_func] = zernike_noll(x_norm, y_norm, j);
    Z(:, j) = val;  % Zernike values at each node
end

act = case_control_sets.subcase_0_SET_1;
n_act = numel(act);
U = zeros(n_nodes, n_act);

for i = 1:n_act
    F_act = zeros(ndof, 1);
    F_act((act(i)-1)*6+3) = 1; % Força unitària al DOF uz
    
    F_act_N = F_act(in_n);
    u_act_N = K_NN \ F_act_N; %dofs lliures
    u_act = zeros(ndof,1);
    u_act(in_n) = u_act_N;
    U(:, i) = u_act(3:6:end); % uz de tots els nodes
end

C = (Z \ U) * 1e6;   % de mm/N a nm/N (Influence matrix)
disp('Coeficients Zernike (modes 4-10) per actuador 13 [um/N]:')
disp(C(4:10, 13) *1e-3) % nm/N -> um/N


%% Apartat 3
% Least-squares method
% Aberracions Zt [nm-RMS]. Enunciat
Zt = zeros(n_modes, 1);
Zt(4)  = 1300;  % Defocus
Zt(5)  = -650;  % Astigmatism (Oblique)
Zt(6)  = 650;   % Astigmatism (Vertical)
Zt(7)  = -350;  % Coma (Vertical)
Zt(8)  = 350;   % Coma (Horizontal)
Zt(9)  = -180;  % Trefoil (Vertical)
Zt(10) = 180;   % Trefoil (Oblique)
Zt(11) = -120;  % Spherical (Primary)
Zt(12) = 70;    % Sec. Astigmatism (Vertical)
Zt(13) = -70;   % Sec. Astigmatism (Oblique)
Zt(14) = -45;   % Quadrafoil (Vertical)
Zt(15) = 45;    % Quadrafoil (Oblique) 

idx_opt = 4:n_modes; % No considerar modes 1, 2, 3
f_opt = C(idx_opt, :) \ Zt(idx_opt); 
disp('Forces actuadors 13 a 17 [N]:')
disp(f_opt(13:17))

% Residu RMS [nm]
residu_vec = C(idx_opt, :) * f_opt - Zt(idx_opt);
residual_surface = Z(:,idx_opt) * ((C(idx_opt,:) * f_opt) * 1e-6) - Z(:,idx_opt) * (Zt(idx_opt) * 1e-6);
total_residual_rms = sqrt(mean(residual_surface.^2)) * 1e6 % [nm]
residual_100_modes_rms = sqrt(mean(residu_vec.^2)) %[nm]
high_freq_residual = norm(Zt(101:end)) %[nm]

%% Apartat 4
nmodes_ap4 = 10;

% Unconstrained
[V_free, D_free] = eigs(K, M, nmodes_ap4, 'smallestabs');
lambda_free = diag(D_free); % lambda = w^2
freq_free = sqrt(abs(lambda_free)) / (2*pi);   % Hz
freq_free = sort(freq_free);

disp('Frequencies unconstrained [Hz]:')
disp(freq_free)

% Comprovació: haurien d'aparèixer 6 modes rígids aprox 0 Hz
n_rigid = sum(freq_free < 1e-3);
fprintf('Rigid body modes detectats: %d\n', n_rigid);

% Constrained (només DOF lliures)
M_NN = M(in_n, in_n);
[V_cons, D_cons] = eigs(K_NN, M_NN, nmodes_ap4, 'smallestabs');
lambda_cons = diag(D_cons);
freq_cons = sqrt(abs(lambda_cons)) / (2*pi);   % Hz

% Ordenar
[freq_cons, idx] = sort(freq_cons);
V_cons = V_cons(:, idx);

disp('Frequencies constrained [Hz]:')
disp(freq_cons)

% 6è mode restringit
mode6 = zeros(ndof,1);
mode6(in_n) = V_cons(:,6);

fprintf('Freqüència del 6è mode: %.6f Hz\n', freq_cons(6));

% Convertir a format nodal [X Y Z RX RY RZ]
u_mode6 = zeros(n_nodes,6);
u_mode6(:,1) = mode6(1:6:end);
u_mode6(:,2) = mode6(2:6:end);
u_mode6(:,3) = mode6(3:6:end);
u_mode6(:,4) = mode6(4:6:end);
u_mode6(:,5) = mode6(5:6:end);
u_mode6(:,6) = mode6(6:6:end);

% Normalització per visualització
u_mode6 = u_mode6 / max(abs(u_mode6(:)));

% Exportar a META
fillhdf('h5template.h5','mode6_constrained.h5',u_mode6);


%% Apartat 5 - Resposta dinàmica

f_range = 0:10:1000;   % Hz
zeta = 0.002;          % 0.2%

% Modes del sistema restringit
[V, D] = eigs(K_NN, M_NN, n_modes, 'smallestabs');
lambda = diag(D);
omega_n = sqrt(lambda);   % rad/s

% DOF vertical de l'actuador 13
dof_13_global = (act(13)-1)*6 + 3;
dof_13 = find(in_n == dof_13_global);

% Força de 1 N al mateix punt
F = zeros(length(in_n),1);
F(dof_13) = 1;

% Força modal
F_modal = V' * F;
H = zeros(size(f_range));

for k = 1:length(f_range)
    omega = 2*pi*f_range(k);

    % Resposta modal
    xi = F_modal ./ (omega_n.^2 - omega^2 + 1i*2*zeta*omega_n*omega);
    % Tornar a espai físic
    u = V * xi;
    % Compliment dinàmic
    H(k) = abs(u(dof_13)) * 1e3;   % um/N
end

Amp_dB = 20*log10(H / H(1)); %en dB

% Bandwidth +3 dB
bw_idx = find(Amp_dB >= 3, 1);
bandwidth_Hz = f_range(bw_idx);

fprintf('Max amplification: %.4f dB\n', max(Amp_dB));
fprintf('Static compliance H(0): %.6f um/N\n', H(1));
fprintf('Mechanical bandwidth (+3 dB): %.2f Hz\n', bandwidth_Hz);

