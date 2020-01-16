from modules import *



run_infall_finder = False
run_infaller_branches = False


all_infalling_objects = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'Infalling_Groups/MergerTreeAHF_Infalling_Re-written/all_infalling'
        + '_objects/')
all_infalling_branches = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'Infalling_Groups/MergerTreeAHF_Infalling_Re-written/all_infaller'
        + '_branches/')
all_infalling_bounded = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'Infalling_Groups/MergerTreeAHF_Infalling_Re-written/all_bounded_at'
        + '_infall_ids/')
halo_data = ('/run/media/ppxrh2/166AA4B87A2DD3B7/NewMDCLUSTER_data/reduced_'
        + 'cluster_info/')



def find_infalling(c, outdir):
    """ Find all the objects which pass within R200 of the cluster, 
    having been outside of R200 at a previous snapshot. Does not 
    account for double infalls, sub-subhaloes, etc. """

    loaddir = ('/run/media/ppxrh2/166AA4B87A2DD3B7/NewMDCLUSTER_data/reduced_'
            'cluster_info/')
    xs_in  = h5py.File(loaddir+'xs/CLUSTER_%04d_xs' % c)
    ys_in  = h5py.File(loaddir+'ys/CLUSTER_%04d_ys' % c)
    zs_in  = h5py.File(loaddir+'zs/CLUSTER_%04d_zs' % c)
    r200_in = h5py.File(loaddir+'rvirs/CLUSTER_%04d_rvirs' % c)

    halo_ids = ld_arr('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
            'MergerTreeAHF_General_Tree_Comp/NewMDCLUSTER_%04d/snap_128/' % c
            + 'CLUSTER_%04d.txt' % c, dtype='int')[c_ids[c-1]]

    halo_ids = halo_ids[halo_ids > 0]
    keys = np.char.mod(u'%03d', halo_ids // mod)


    j = keys[0]
    xs = np.array(xs_in[j])
    ys = np.array(ys_in[j])
    zs = np.array(zs_in[j])
    ids = np.array(np.arange(len(xs)), dtype='int')+int(j)*mod+1
    r200s = np.array(r200_in[j])
    
    h_id = halo_ids[0] - (int(j)*mod+1)
    clust = [xs[h_id], ys[h_id], zs[h_id], r200s[h_id]]
        
    #ids_out: id of every object outside the cluster in the first snapshot
    xs, ys, zs = xs-clust[0], ys-clust[1], zs-clust[2]
    rs = (xs**2. + ys**2. + zs**2.)**0.5 / clust[3]
    ids_out = ids[rs > 1.]


    ids_infall_t = np.array(np.zeros(0), dtype='int')
    
    for i in range(1, len(keys)):
        j = keys[i]

        xs = np.array(xs_in[j])
        ys = np.array(ys_in[j])
        zs = np.array(zs_in[j])
        r200s = np.array(r200_in[j])
        ids = np.array(np.arange(len(xs)), dtype='int')+int(j)*mod+1
        ids_main_p_in = ld_arr('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTree'
                'AHF/MergerTreeAHF_HDF5_MainProg/CLUSTER_%04d_CRMratio2_Main'
                'Prog/CLUSTER_%04d_prog_output.snapshot_' % (c, c) + j
                + '.AHF.MainProg.txt', dtype='int')
        j = int(j)
        ids_main_p = np.array(np.zeros(len(xs)), dtype='int')
        ids_main_p[ids_main_p_in[:, 0]-(j*mod+1)] = ids_main_p_in[:, 1]
        
        h_id = halo_ids[i] - (int(j)*mod+1)
        clust = [xs[h_id], ys[h_id], zs[h_id], r200s[h_id]]

        xs, ys, zs = xs-clust[0], ys-clust[1], zs-clust[2]
        rs = (xs**2. + ys**2. + zs**2.)**0.5 / clust[3]


        #ids_infall: ids of haloes that are now in r200
        #ids_prev: ids of where haloes in r200 previously were located 
        #ids_out: haloes in ALL previous snapshots which were outside r200 
        ids_infall = ids[rs < 1.]
        ids_prev = ids_main_p[ids_infall-(j*mod+1)]

        ids_infall = ids_infall[ids_prev > 0]
        ids_prev = ids_prev[ids_prev > 0]
        for k in range(len(ids_prev)):
            if len(np.where(ids_out==ids_prev[k])[0])==0:
                ids_prev[k] = 0


        #ids_infall: ids of every halo which has JUST fallen into the r200
        ids_infall = ids_infall[ids_prev>0]
        ids_infall_t = np.append(ids_infall_t, ids_infall)

        ids_out = np.append(ids_out, ids[rs > 1.])

    pd.DataFrame(np.char.mod('%15d', ids_infall_t)).to_csv(outdir + 
            '/CLUSTER_%04d_all_infallers.txt' % c, 
            header=['  all_infallers'], index=None)
    

    return None


def find_infaller_evolution(c, indir, outdir):
    """ Finds the subsequent evolution of every object that has fallen 
    into the cluster. Objects with multiple infalls are eliminated, 
    such that only the ID at first infall is recorded. The total number 
    of infalls is stored """

    total_tree = np.array(pd.read_csv('/run/media/ppxrh2/166AA4B87A2DD3B7/'
            'MergerTreeAHF/MergerTreeAHF_ASCII/MergerTree_GadgetX-NewMD'
            'CLUSTER_%04d.txt-CRMratio2' % c, sep='\s+', skiprows=2, 
            usecols=[0], dtype='str')[:-1], dtype='int')[:, 0]
    
    ids_infall_t = ld_arr(indir + '/CLUSTER_%04d_all_infallers.txt' % c, 
            dtype='int')[:, 0]

    hf = h5py.File(outdir + 'cluster_%04d_groups.hdf5' % c, 'w')
    infall_number = np.array(np.ones(len(ids_infall_t)), dtype='int')
    i_nonrep = 0
    for i in range(len(ids_infall_t)):
        branch = find_branch(np.array([ids_infall_t[i]]),total_tree)[0]
        for j in branch[1:]:
            repeat = np.where(ids_infall_t == j)[0]
            if len(repeat) == 1:
                infall_number[repeat] = infall_number[i] + 1
        if len(branch) > 0:
            #label dataset with index, and number of infalls
            hf.create_dataset('%04d_%02d' % (i_nonrep, infall_number[i]), data=branch)
            i_nonrep += 1
    
    hf.close()

    return None








def find_bound_groups(c, loaddir, datadir, outdir):
    """ Take the full list of all infalling objects (regardless of 
    mass, etc.) and find which objects are bound to them at infall """
    
    infalls = h5py.File(loaddir + 'cluster_%04d_groups.hdf5' % c, 'r')
    keys = np.array(list(infalls.keys()), dtype='str')

    if len(keys)>0:
        #splt key strings at underscore, into two keys
        n_infalls = np.array(np.core.defchararray.partition(keys, '_')[:, 2],
                dtype='int') #infall number for each object
    else:
        n_infalls = np.zeros(0)
    
    #ids of everything infalling
    infall_ids_t = np.array([np.array(infalls[keys[i]], dtype='int')[0] 
            for i in np.arange(len(keys))])
    
    infall_snaps = infall_ids_t // mod
    infall_ids = infall_ids_t - (infall_snaps*mod+1) #zero index
    #snapshots at which objects infall, and number of objects at each
    infall_snaps = np.array(list(Counter(infall_snaps).keys()))
    infall_snaps_ct = np.array(list(Counter(infall_snaps).values()))


    xs = h5py.File(datadir + 'xs/CLUSTER_%04d_xs' % c, 'r')
    ys = h5py.File(datadir + 'ys/CLUSTER_%04d_ys' % c, 'r')
    zs = h5py.File(datadir + 'zs/CLUSTER_%04d_zs' % c, 'r')
    vx = h5py.File(datadir + 'vx/CLUSTER_%04d_vx' % c, 'r')
    vy = h5py.File(datadir + 'vy/CLUSTER_%04d_vy' % c, 'r')
    vz = h5py.File(datadir + 'vz/CLUSTER_%04d_vz' % c, 'r')
    ms = h5py.File(datadir + 'ms/CLUSTER_%04d_ms' % c, 'r')
    mstars = h5py.File(datadir + 'mstars/CLUSTER_%04d_mstars' % c, 'r')
    r200s = h5py.File(datadir + 'rvirs/CLUSTER_%04d_rvirs' % c, 'r')

    print(alphabet)

    rs_out, vs_out = np.zeros(0), np.zeros(0) 
    ids_out = np.array(np.zeros(0), dtype='int')
    hosts_out = np.array(np.zeros(0), dtype='int')
    counter = 0
    bound_out = h5py.File(outdir+'cluster_%04d_members.hdf5' % c, 'w')
    for s in range(len(infall_snaps)):
        snap = '%03d' % infall_snaps[s]
        xs_s, ys_s, zs_s, r200s_s = xs[snap], ys[snap], zs[snap], r200s[snap]
        vx_s, vy_s, vz_s = vx[snap], vy[snap], vz[snap]
        ms_s = np.array(ms[snap], dtype='int') 
        mstars_s = np.array(mstars[snap], dtype='int')
        ratio = (mstars_s+1) / ms_s
        ids_s = np.array(np.arange(len(xs_s)), dtype='int')

        boolean = (ms_s >= m_lim) * (mstars_s >= ms_lim) * (ratio <= rat_lim)
        xs_s, ys_s, zs_s = xs_s[boolean], ys_s[boolean], zs_s[boolean]
        vx_s, vy_s, vz_s = vx_s[boolean], vy_s[boolean], vz_s[boolean]
        ms_s, mstars_s = ms_s[boolean], mstars_s[boolean]
        ids_s, r200s_s = ids_s[boolean], r200s_s[boolean]
        
        for i in range(infall_snaps_ct[s]):
            id_n = infall_ids[counter]
            id_n = np.where(ids_s == id_n)[0][0]
            rs_si = ((xs_s-xs_s[id_n])**2. + (ys_s-ys_s[id_n])**2.
                    + (zs_s-zs_s[id_n])**2.)**0.5
            vs_si = ((vx_s-vx_s[id_n])**2. + (vy_s-vy_s[id_n])**2.
                    + (vz_s-vz_s[id_n])**2.)**0.5
            ms_si = (ms_s < ms_s[id_n])
            ids_si, rs_si, vs_si = ids_s[ms_si], rs_si[ms_si], vs_si[ms_si]
            boundres = bound(vs_si, rs_si, ms_s[id_n], r200s_s[id_n])
            boundres_bool = boundres[0]
            if len(np.where(ids_out==infall_ids_t[counter])[0]) > 0:
                ids_bound = np.array(np.zeros(0), dtype='int')
                boundres_bool = np.zeros(len(boundres_bool))<-1.
            rs_out = np.append(rs_out, boundres[1][boundres_bool])
            vs_out = np.append(vs_out, boundres[2][boundres_bool])
            ids_bound = infall_snaps[s]*mod+1+ids_si[boundres_bool]
            ids_out = np.append(ids_out, ids_bound)
            hosts_out = np.append(hosts_out, np.array(
                    [infall_ids_t[counter]]*len(ids_bound)))
            bound_out.create_dataset('%16d_%03d' % (infall_ids_t[counter], 
                    len(ids_bound)), data=ids_bound)
            counter += 1
    bound_out.close()
    
    return rs_out, vs_out, ids_out, hosts_out














def find_branch(id_a, id_list, keep_sing=False):
    """ Finds the full evolution of a halo, starting with the earliest 
    id """
    id = id_a[-1]
    final_state = 0
    tree_search = np.where(id_list==id)[0]
    if len(tree_search)==0:
        result = np.array([id]*keep_sing, dtype='int')
        final_state = 1 #no progenitor - should not be the case for these
        #only the case when id_a is '0', for second infall
    elif id//mod == id_list[tree_search[0]-1]//mod:
        result = id_a
        final_state = 2 #absorbed by a larger halo - not main progenitor
    elif len(tree_search)==1:
        result = id_a
        final_state = 1 #lost by halo finder
    else:
        id_a = np.append(id_a, id_list[tree_search[0]-1])
        result = find_branch(id_a, id_list)[0]
    
    return result, final_state



if run_infall_finder==True:
    for i in range(1, 325):
        print(i)
        find_infalling(i, all_infalling_objects)


if run_infaller_branches==True:
    for i in range(1, 325):
        print(i)
        find_infaller_evolution(i, all_infalling_objects, 
                all_infalling_branches)



find_bound_groups(1, all_infalling_branches, halo_data, all_infalling_bounded)