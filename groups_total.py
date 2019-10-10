from modules import *


def find_groups(c, m_lim=0., ms_lim=0., rat_lim=1.):
    """ Find all the objects which pass within R200 of the cluster, having 
    been outside of R200 at a previous snapshot """

    loaddir = ('/run/media/ppxrh2/166AA4B87A2DD3B7/NewMDCLUSTER_data/reduced_'
            'cluster_info/')
    xs_in  = h5py.File(loaddir+'xs/CLUSTER_%04d_xs' % c)
    xs_in  = h5py.File(loaddir+'xs/CLUSTER_%04d_xs' % c)
    ys_in  = h5py.File(loaddir+'ys/CLUSTER_%04d_ys' % c)
    zs_in  = h5py.File(loaddir+'zs/CLUSTER_%04d_zs' % c)
    ms_in  = h5py.File(loaddir+'ms/CLUSTER_%04d_ms' % c)
    mstar_in  = h5py.File(loaddir+'mstars/CLUSTER_%04d_mstars' % c)
    r200_in = h5py.File(loaddir+'rvirs/CLUSTER_%04d_rvirs' % c)

    halo_ids = ld_arr('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
            'MergerTreeAHF_General_Tree_Comp/NewMDCLUSTER_%04d/snap_128/' % c
            + 'CLUSTER_%04d.txt' % c, dtype='int')[c_ids[c-1]]

    halo_ids = halo_ids[halo_ids > 0]
    keys = np.char.mod(u'%03d', halo_ids / mod)


    j = keys[0]
    xs = np.array(xs_in[j])
    ys = np.array(ys_in[j])
    zs = np.array(zs_in[j])
    ms = np.array(ms_in[j])+1.
    ids = np.array(np.arange(len(xs)), dtype='int')+int(j)*mod+1
    mstars = np.array(mstar_in[j])
    r200s = np.array(r200_in[j])
    
    h_id = halo_ids[0] - (int(j)*mod+1)
    clust = [xs[h_id], ys[h_id], zs[h_id], r200s[h_id]]
    ratio = mstars / ms

    m_bool = (ms >= m_lim) * (mstars >= ms_lim) * (ratio <= rat_lim)
    m_bool[h_id] = False
    xs, ys, zs, ids = xs[m_bool], ys[m_bool], zs[m_bool], ids[m_bool]
        
    xs, ys, zs = xs-clust[0], ys-clust[1], zs-clust[2]
    rs = (xs**2. + ys**2. + zs**2.)**0.5 / clust[3]
    ids_out = ids[rs > 1.]

    ids_infall_t = np.array(np.zeros(0), dtype='int')
    
    for i in range(1, len(keys)):
        j = keys[i]
        #print j
        rs_p = rs

        xs = np.array(xs_in[j])
        ys = np.array(ys_in[j])
        zs = np.array(zs_in[j])
        ms = np.array(ms_in[j])+1.
        mstars = np.array(mstar_in[j])
        r200s = np.array(r200_in[j])
        ids = np.array(np.arange(len(xs)), dtype='int')+int(j)*mod+1
        ids_main_p_in = ld_arr('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTree'
                'AHF/MergerTreeAHF_HDF5_MainProg/CLUSTER_%04d_CRMratio2_Main'
                'Prog/CLUSTER_%04d_prog_output.snapshot_' % (c, c) + j
                + '.AHF.MainProg.txt', dtype='int')
        j = int(j)
        ids_main_p = np.array(np.zeros(ids_main_p_in[-1, 0]-j*mod), 
                dtype='int')
        ids_main_p[ids_main_p_in[:, 0]-(j*mod+1)] = ids_main_p_in[:, 1]

        h_id = halo_ids[i] - (int(j)*mod+1)
        clust = [xs[h_id], ys[h_id], zs[h_id], r200s[h_id]]
        ratio = mstars / ms

        m_bool = (ms >= m_lim) * (mstars >= ms_lim) * (ratio <= rat_lim)
        m_bool[h_id] = False
        xs, ys, zs, ids = xs[m_bool], ys[m_bool], zs[m_bool], ids[m_bool]
        
        xs, ys, zs = xs-clust[0], ys-clust[1], zs-clust[2]
        rs = (xs**2. + ys**2. + zs**2.)**0.5 / clust[3]
        

        ids_infall = ids[rs < 1.]
        ids_prev = ids_main_p[ids_infall-(j*mod+1)]
        #ids_prev: ids of where haloes in r200 previously were located 
        #ids_out: haloes in ALL previous snapshots which were outside r200 
        
        for k in range(len(ids_prev)):
            if len(np.where(ids_out==ids_prev[k])[0])==0:
                ids_prev[k] = 0
        
        #ids_infall: ids of every halo which has JUST fallen into the r200
        ids_infall = ids_infall[ids_prev>0]
        ids_infall_t = np.append(ids_infall_t, ids_infall)

        ids_out = np.append(ids_out, ids[rs > 1.])
    
    total_tree = np.array(pd.read_csv('/run/media/ppxrh2/166AA4B87A2DD3B7/'
            'MergerTreeAHF/MergerTreeAHF_ASCII/MergerTree_GadgetX-NewMD'
            'CLUSTER_%04d.txt-CRMratio2' % c, sep='\s+', skiprows=2, 
            usecols=[0], dtype='str')[:-1], dtype='int')[:, 0]

    hf = h5py.File('/run/media/ppxrh2/166AA4B87A2DD3B7/Other_Results/complete'
            '_infall_groups_list//cluster_%04d_groups.hdf5' % c, 'w')
    i_nonrep = 0
    for i in range(len(ids_infall_t)):
        #print i
        id_temp = ids_infall_t[i]
        ids_infall_t[i] = 0
        
        branch = find_branch(np.array([id_temp]), total_tree)
        n_infall = 1
        for j in branch[1:]:
            repeat = np.where(ids_infall_t == j)[0]
            if len(repeat) == 1:
                n_infall += 1
                ids_infall_t[repeat[0]] = 0
        if len(branch) > 0:
            #label dataset with index, and number of infalls
            hf.create_dataset('%04d_%02d' % (i_nonrep, n_infall), data=branch)
            i_nonrep += 1
    
    hf.close()

    return None


def find_branch(id_a, id_list):
    id = id_a[-1]
    tree_search = np.where(id_list==id)[0]
    if len(tree_search)==0:
        result = np.array(np.zeros(0), dtype='int')
    elif id/mod == id_list[tree_search[0]-1]/mod:
        result = id_a
    elif len(tree_search)==1:
        result = id_a
    else:
        id_a = np.append(id_a, id_list[tree_search[0]-1])
        result = find_branch(id_a, id_list)
    
    return result



#for c_id in range(1, 325):
#    print c_id
#    find_groups(c_id, 10.**10.5, 10.**9.5, 0.3)


