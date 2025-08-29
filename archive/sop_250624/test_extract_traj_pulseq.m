clear;
seqPath = '/Users/cag/Documents/forclone/pulseq4mreye/debug_0624/seq1_t1w_libre_debugTR6p2_rfdelay_gain_0.seq';

nSeg=22;
nShot=1000;
N = 480;
FoV   = 240; 
flagSelfNav = 1; 
nShot_off   = 14;


k_trj = extract_traj_pulseq(seqPath, nSeg, nShot, N);
k_trj = k_trj/max(abs(k_trj(:)))*0.5*N/FoV;

if flagSelfNav
k_trj(:,:, 1, :) = [];
end
k_trj(:,:, :, 1:nShot_off) = [];


k_trj_size = size(k_trj);
k_trj = reshape(k_trj, [3,N,k_trj_size(3)*k_trj_size(4)]);

t = bmTraj_fullRadial3_phyllotaxis_lineAssym2(N,nSeg, nShot, 1/FoV, flagSelfNav,nShot_off);
%%
check_traj_diff(k_trj, t, 1e-4); 

function k_trj = extract_traj_pulseq(seqPath, nSeg, nShot, N)
   
    
    seq=mr.Sequence();
    disp(['Reading seq file: ', seqPath]);
    seq.read(seqPath);
    
    kspace_traj = seq.calculateKspacePP();
    
    k_trj = reshape(kspace_traj, [3, N, nSeg, nShot]);

end




function isclose = check_traj_diff(t_monalisa, t_pulseq, tol) 
    
    fprintf('Setting tolerence: %f \n', tol);

    % Compute the difference
    diff = abs(squeeze(t_monalisa(:,:,:)) - squeeze(t_pulseq(:,:,:)));
    % Check if all differences are below tolerance
    if all(diff(:) < tol)
        disp('All elements are close within the specified tolerance.');
        isclose=1;
    else
        isclose=0;
        disp('Some elements differ beyond the tolerance.');    
        % Optionally, show where and how much
        [i, j] = find(diff >= tol);
        for idx = 1:length(i)
            fprintf('Mismatch at (%d, %d):  -t = %g, ref = %g, diff = %g\n', ...
                i(idx), j(idx), ...
                t_monalisa(1,i(idx),j(idx)), ...
                t_pulseq(1,i(idx),j(idx)), ...
                diff(i(idx),j(idx)));
        end
        
    end

end