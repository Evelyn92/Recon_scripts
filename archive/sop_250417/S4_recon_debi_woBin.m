clear;clc; 
% =====================================================
% Author: Yiwei Jia
% Date: April 03
% ------------------------------------------------
% Recon on debi
% Update: just for test before running on HPC
% =====================================================

% change woBinning/Binning, change th75->th0, change eyeMask name!

%%
cs_recon=1;
subject_array = [2,3,4];
for subject_num = subject_array


    subject_num_C = subject_num;
    
    datasetDir = '/home/debi/yiwei/mreye_dataset/250417/';
    reconDir = '/home/debi/yiwei/recon_results/250417/';
    
    mask_note_list={'placeholder','libre_t1w', ...
        'libre_t2w', 'gre_original', 'gre_unishuffle'};


    mask_note = mask_note_list{subject_num};

    x0Dir = [reconDir, '/Sub00',num2str(subject_num),'/T1_LIBRE_woBinning/output/mask_',mask_note,'/'];
    xDir = [reconDir, '/Sub00',num2str(subject_num),'/T1_LIBRE_woBinning/output/mask_',mask_note,'/'];



    mDir = [reconDir, '/Sub00',num2str(subject_num),'/T1_LIBRE_woBinning/mitosius/mask_', mask_note, '/'];



    saveCDirList = {['/Sub00',num2str(subject_num),'/T1_LIBRE_Binning/C/'],['/Sub00',num2str(subject_num),'/T1_LIBRE_woBinning/C/']};


    CfileName = 'C.mat';
    saveCDir = [reconDir, saveCDirList{2}];
    CfilePath = fullfile(saveCDir, CfileName);


    %
    y   = bmMitosius_load(mDir, 'y'); 
    t   = bmMitosius_load(mDir, 't'); 
    ve  = bmMitosius_load(mDir, 've'); 
    
    disp('Mitosius has been loaded!')
    disp(mDir)
    %% compileScript()
    Matrix_size = 240;
    ReconFov = 480; %mm
    N_u     = [Matrix_size, Matrix_size, Matrix_size]; % Matrix size: Size of the Virtual cartesian grid in the fourier space (regridding)
    n_u     = [Matrix_size, Matrix_size, Matrix_size]; % Image size (output)
    dK_u    = [1, 1, 1]./ReconFov; % Spacing of the virtual cartesian grid
    nFr     = 1; 
    % best achivable resolution is 1/ N_u*dK_u If you have enough coverage
%%

    load(CfilePath); 
    C = bmImResize(C, [48, 48, 48], N_u);

%

x0 = cell(nFr, 1);
for i = 1:nFr
    x0{i} = bmMathilda(y{i}, t{i}, ve{i}, C, N_u, n_u, dK_u, [], [], [], []);
end
% isequal(x0_p, x0)
%
bmImage(x0);

%
if ~isfolder(x0Dir)
    % If it doesn't exist, create it
    mkdir(x0Dir);
    disp(['Directory created: ', x0Dir]);
else
    disp(['Directory already exists: ', x0Dir]);
end
x0Path = fullfile(x0Dir, 'x0.mat');
% Save the x0 to the .mat file
save(x0Path, 'x0', '-v7.3');
disp('x0 has been saved here:')
disp(x0Path)

%%
if cs_recon

[Gu, Gut] = bmTraj2SparseMat(t, ve, N_u, dK_u);
% bmSteva
deltaArray = 1;

% nIter = 30; % iterations before stopping
nIter = 15; %20, 30
witness_ind = [15,18];
delta = deltaArray(1);
% delta     = 0.1; %0.01, 0.1, 1
rho       = 10*delta;
nCGD      = 4;
ve_max    = 10*prod(dK_u(:));
if nFr<= 1

    x = bmSteva(  x0{1}, [], [], y{1}, ve{1}, C, Gu{1}, Gut{1}, n_u, ...
                                        delta, rho, nCGD, ve_max, ...
                                        nIter, ...
                                        bmWitnessInfo('steva_d1_r10_nCGD4', witness_ind));
else
    witness_info = bmWitnessInfo('tevaMorphosia_d0p1_r1_nCGD4', witness_ind);
    witness_info.save_witnessIm_flag = true;
    x = bmTevaMorphosia_chain(  x0, ...
                                [], [], ...
                                y, ve, C, ...
                                Gu, Gut, n_u, ...
                                [], [], ...
                                delta, rho, 'normal', ...
                                nCGD, ve_max, ...
                                nIter, ...
                                witness_info);
end

bmImage(x)




if ~isfolder(xDir)
    % If it doesn't exist, create it
    mkdir(xDir);
    disp(['Directory created: ', xDir]);
else
    disp(['Directory already exists: ', xDir]);
end

xPath = fullfile(xDir, sprintf('x_nIter%d_delta_%.3f.mat', nIter, delta));

% Save the x0 to the .mat file
save(xPath, 'x');
disp('x has been saved here:')
disp(xPath)

end

end