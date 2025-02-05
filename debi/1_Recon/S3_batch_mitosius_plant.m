mode = '4_frs';
mode_list = {'discard', '4_frs'};
subject_num = 2;
reconDir = '/media/sinf/1,0 TB Disk/240922_480/';

if strcmp(mode, mode_list{1})
    
    if subject_num == 1
        c_path = [reconDir, '/Sub001/T1_LIBRE_Binning/C/C.mat'];
        cMask_path = [reconDir, '/Sub001/T1_LIBRE_Binning/other/discard_cMask_th2_12831.mat'];
        mDir = '/media/sinf/1,0 TB Disk/240922_480/Sub001/T1_LIBRE_Binning/mitosius/discard_12831/';
    elseif subject_num == 2
        c_path = [reconDir, '/Sub002/T1_LIBRE_Binning/C/C.mat'];
        cMask_path = [reconDir, '/Sub002/T1_LIBRE_Binning/other/discard_cMask_th2_11556.mat'];
        mDir = '/media/sinf/1,0 TB Disk/240922_480/Sub002/T1_LIBRE_Binning/mitosius/discard_11556/';
    elseif subject_num == 3
        c_path = [reconDir, '/Sub003/T1_LIBRE_Binning/C/C.mat'];
        cMask_path = [];
        mDir = '';
    else
        c_path = [reconDir, '/Sub004/T1_LIBRE_Binning/C/C.mat'];
        cMask_path = [reconDir, '/Sub004/T1_LIBRE_Binning/other/discard_cMask_th2_12912.mat'];
        mDir = '/media/sinf/1,0 TB Disk/240922_480/Sub004/T1_LIBRE_Binning/mitosius/discard_12912/';
    end
elseif strcmp(mode, mode_list{2})
    if subject_num == 1
        c_path = [reconDir, '/Sub001/T1_LIBRE_Binning/C/C.mat'];
        cMask_path = [reconDir, '/Sub001/T1_LIBRE_Binning/other/cMask_4fr.mat'];
        mDir = '/media/sinf/1,0 TB Disk/240922_480/Sub001/T1_LIBRE_Binning/mitosius/4fr/';
    elseif subject_num == 2
        c_path = [reconDir, '/Sub002/T1_LIBRE_Binning/C/C.mat'];
        cMask_path = [reconDir, '/Sub002/T1_LIBRE_Binning/other/cMask_4fr.mat'];
        mDir = '/media/sinf/1,0 TB Disk/240922_480/Sub002/T1_LIBRE_Binning/mitosius/4fr/';
    elseif subject_num == 3
        c_path = [reconDir, '/Sub003/T1_LIBRE_Binning/C/C.mat'];
        cMask_path = [];
        mDir = '';
    else
        c_path = [reconDir, '/Sub004/T1_LIBRE_Binning/C/C.mat'];
        cMask_path = [reconDir, '/Sub004/T1_LIBRE_Binning/other/cMask_4fr.mat'];
        mDir = '/media/sinf/1,0 TB Disk/240922_480/Sub004/T1_LIBRE_Binning/mitosius/4fr/';
    end
    

end

mitosius_func(subject_num, c_path, cMask_path, mDir);