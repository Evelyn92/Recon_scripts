%% Yiwei Jia July 14 2025
% The script is specifically used for analyzing Twix 

%%
clc; clear all;

raw_data = ['/home/debi/yiwei/mreye_dataset/250624/' ...
    '/meas_MID00388_FID34125_pulseqv15.dat'];
twix = mapVBVD_JB(raw_data);
seqHash_twix = twix{end}.hdr.Dicom.tSequenceVariant;
%%
if length(seqHash_twix)==32
    fprintf(['raw data contain pulseq-file signature ' seqHash_twix '\n']);
end
%%
seqFile = ['/home/debi/yiwei/forclone/pulseqmreye/debug_0624/' ...
    '/seq1_t1w_libre_debugTR6p2_rfdelay_gain_0.seq'];

% seq1_t1w_gre_debug_TR6p2wo_gzsp_cor_te_wrap_r
% seq1_t1w_libre_debugTR6p2_rfdelay_gain_0
% seq2_t1w_gre_debug_TR8p01_rfdelay_gdsp
% seq2_t1w_libre_debugTR8p01_rfdelay_gdsp_gain_0

% Read the entire content of the .seq file
raw = fileread(seqFile);

% Find the position of the '[SIGNATURE]' section
sigPos = strfind(raw, '[SIGNATURE]');

if ~isempty(sigPos)
    % Exclude the '[SIGNATURE]' section and the preceding newline
    raw = raw(1:sigPos(1)-2);
end

% Create a MessageDigest object for MD5
md = java.security.MessageDigest.getInstance('MD5');

% Update the digest with the content
md.update(uint8(raw));

% Compute the hash
hash = md.digest();

% Convert the hash to a hexadecimal string
hash_hex = dec2hex(typecast(hash, 'uint8'))';
hash_str = lower(reshape(hash_hex, 1, []));

% Display the recalculated hash
disp(['Recalculated MD5 hash: ', hash_str]);