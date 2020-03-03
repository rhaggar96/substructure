from modules import *


inputdir_all = ('/run/media/ppxrh2/166AA4B87A2DD3B7/NewMDCLUSTER_data/'
        + 'NewMDCLUSTER_0002_reruns/NewMDCLUSTER_%04d/')
outputdir = ('/run/media/ppxrh2/166AA4B87A2DD3B7/NewMDCLUSTER_data/'
        + 'NewMDCLUSTER_0002_reruns/reduced_cluster_info/')



def get_ahf_h_fname(dir, snap_inq):
    """ Get AHF halos file name for correct snapshot """
    
    snap_inq = '%03d' % snap_inq
    for file in os.listdir(dir):
        if file.endswith('_halos'):
            fnames = os.path.join(dir, file)
            if fnames[-20:-17] == snap_inq or fnames[-21:-18] == snap_inq:
                fname = fnames
    return fname

def extract_halo_positions(c_no, suffs=['subs', 'ms', 'npars', 'xs', 'ys', 
        'zs', 'vx', 'vy', 'vz', 'rvirs', 'mgas', 'mstars']):
    """ Take halo catalogues, and extract positions, virial radii, 
    masses and subhalo status to separate files """
    
    inputdir = inputdir_all % c_no

    for i in range(len(suffs)):
        if not os.path.exists(outputdir+suffs[i]+'/'):
            os.mkdir(outputdir+suffs[i]+'/')
    
    path_list = []
    for i in range(len(suffs)):
        path_list += [outputdir + suffs[i] + '/CLUSTER_%04d_'%c_no + suffs[i]]
        if os.path.exists(path_list[i]):
            os.remove(path_list[i])

    hf_list = []
    for i in range(len(suffs)):
        hf_list += [h5py.File(path_list[i], 'w')]
    int_rng = [0,1,2,10,11]
    flo_rng = range(3, 10)

    for i in range(129):
        print(i)
        if (c_no==10 and i==100) or (c_no==228 and i==110):
            data_in = np.zeros(0)
            
        else:
            inputname = get_ahf_h_fname(inputdir, i)
            keys_list = ['hostHalo(2)', 'numSubStruct(3)', 'Mvir(4)', 
                    'npart(5)', 'Xc(6)', 'Yc(7)', 'Zc(8)', 'VXc(9)', 'VYc(10)',
                    'cNFW(43)', 'Epot_gas(63)']
            ###DM-ONLY CASE: C_NO=2000###
            if c_no==2000:
                keys_list = ['hostHalo(2)', 'Mvir(4)', 'npart(5)', 'Xc(6)', 
                'Yc(7)', 'Zc(8)', 'VXc(9)', 'VYc(10)', 'VZc(11)', 'Rvir(12)']
            
            data_in = pd.read_csv(inputname, sep='\t', usecols=keys_list)
            subs = np.asarray(list(data_in.index))
            data_in = np.asarray(data_in)

            if len(data_in) > 0:
                data_subs_in = np.zeros((len(data_in), len(suffs)))
                if c_no!=2000:
                    data_subs_in[:, 1:] = data_in
                    data_subs_in[:, 0] = subs[:, 1]
                else:
                    data_subs_in[:, :10] = data_in
                
                data_int = np.array(data_subs_in[:, int_rng], dtype='int')
                data_flo = np.array(data_subs_in[:, flo_rng])
            
        
        if len(data_in) == 0:
            for j in range(len(int_rng)):
                hf_list[int_rng[j]].create_dataset('%03d' % i, 
                        data=np.array(np.zeros(0), dtype='int'))
            for j in range(len(flo_rng)):
                hf_list[flo_rng[j]].create_dataset('%03d' % i, 
                        data=np.array(np.zeros(0), dtype='float'))

        else:
            for j in range(len(int_rng)):
                hf_list[int_rng[j]].create_dataset('%03d' % i, 
                        data=data_int[:, j])
            for j in range(len(flo_rng)):
                hf_list[flo_rng[j]].create_dataset('%03d' % i, 
                        data=data_flo[:, j])

    for hf in hf_list:
        hf.close()
    
    return None 




extract_halo_positions(2000)
extract_halo_positions(2001)
