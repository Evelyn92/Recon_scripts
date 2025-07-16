%% Ludovica
% Discard spokes to reach steady state or to just cut down the dataset
SpokesToDiscard=200; % for steady state
%SpokesToDiscard=length(k0_signals)-counter; % this is to only consider the breath-hold (manual, todo automated)
KSPACE=KSPACE(:,:,SpokesToDiscard:end,:,:);
Times=Times(:,:,SpokesToDiscard:end,:,:);

k=k(:,:,SpokesToDiscard:end,:,:);
w=w(:,:,SpokesToDiscard:end,:,:);

% Compress coils to speed up reconstruction keeping a minimum of 8
cKSPACE = reshape( KSPACE, size(KSPACE,1)*size(KSPACE,2)*size(KSPACE,3), size(KSPACE,4) );
[~,S,V] = svd(cKSPACE,'econ'); S=diag(S);
param.Nc=max(min(sum(param.UseCoils),8),knee_pt(S)); %*
cKSPACE = reshape( cKSPACE * V(:,1:param.Nc), size(KSPACE,1), size(KSPACE,2), size(KSPACE,3), param.Nc );

% KSPACE=cKSPACE;
% clear cKSPACE

% Compute coil sensitivities using enough spokes to meet nyquist
Mask=false(size(cKSPACE,2),1); Mask(1:round(param.Np*1))=true;
E = gpNUFFT((reshape(k(:,:,:,1),[size(kSPACE,1)*size(cKSPACE(:,:,Mask),2)])),...,col(w(:,:,Mask,1)).^-2,1.5,3,8,1*[param.Nx,param.Ny],[],true);
IC=(E'*reshape(kaiser(param.Nx,2)*kaiser(param.Ny,2)',[size(cKSPACE(:,:,Mask),1),size(cKSPACE(:,:,Mask),2),size(cKSPACE,4)]));
CSM=permute(ismrm_estimate_csm_walsh(IC),[1,2,4,3]);


%% Matteo
% Coil Compression
if coilCompression == 1
    D = reshape(kdata_raw, nx * ntviews, nc);
    [U, S, V] = svd(D, 'econ');
    singular_values = diag(S);
    total_variance = sum(singular_values .^ 2);
    explained_variance = singular_values(1:ncc) .^ 2;
    percentage_explained = (explained_variance / total_variance) * 100;
    kdata_raw = reshape(D * V(:, 1:ncc), nx, ntviews, ncc);
    % Plot the explained variance
    f=figure;
    f.Position = [100 100 500 800];
    plot(percentage_explained,'.-','LineWidth',2,'Color','r','MarkerSize',20)
    xlabel('Virtual Coil nr.')
    ylabel('Explained [%]')
    title(sprintf('Coil Compression, Explaination for ncc = %d is %.2f%%',ncc,sum(percentage_explained(1:ncc))))
else
    % coil combination to save memory (for now just coil selection) EP test
    selected_coils = [1,2,4,5,6:2:size(kdata_raw,3)];
    kdata_raw=kdata_raw(:,:,selected_coils);
end