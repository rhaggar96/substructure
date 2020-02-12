from modules import *


run_infall_finder = False
run_infaller_branches = False
run_bounded_at_infall = False
run_all_grp_member_branches = False
run_all_memb_data = False
run_find_cluster_data = False
run_all_memb_data_rel_clus = False


all_infalling_objects = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'Infalling_Groups/MergerTreeAHF_Infalling_Re-written/all_infalling'
        + '_objects/')
all_infalling_branches = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'Infalling_Groups/MergerTreeAHF_Infalling_Re-written/all_infaller'
        + '_branches/')
all_infalling_bounded = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'Infalling_Groups/MergerTreeAHF_Infalling_Re-written/all_bounded_at'
        + '_infall_ids/')
all_grp_member_branches = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'Infalling_Groups/MergerTreeAHF_Infalling_Re-written/all_group_'
        + 'member_branches/')
all_grp_member_data = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'Infalling_Groups/MergerTreeAHF_Infalling_Re-written/'
        + 'all_members_halo_data/absolute/')
all_relative_clus_data = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'Infalling_Groups/MergerTreeAHF_Infalling_Re-written/'
        + 'all_members_halo_data/rel_to_cluster/')
halo_data = ('/run/media/ppxrh2/166AA4B87A2DD3B7/NewMDCLUSTER_data/reduced_'
        + 'cluster_info/')
cluster_halo_data = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/Infall'
        + 'ing_Groups/MergerTreeAHF_Infalling_Re-written/cluster_halo_data/')
backsplash_track_data = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'MergerTreeAHF_General_Tree_Comp/')



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

    hf = h5py.File(outdir + 'CLUSTER_%04d_groups.hdf5' % c, 'w')
    infall_number = np.array(np.ones(len(ids_infall_t)), dtype='int')
    i_nonrep = 0
    for i in range(len(ids_infall_t)):
        branch = find_branch(np.array([ids_infall_t[i]]),total_tree)
        fate = branch[1]
        branch = branch[0]
        for j in branch[1:]:
            repeat = np.where(ids_infall_t == j)[0]
            if len(repeat) == 1:
                infall_number[repeat] = infall_number[i] + 1
        if len(branch) > 0:
            #label dataset with index, final state amd number of infalls
            hf.create_dataset('%04d_%01d_%02d' % (i_nonrep, fate, 
                    infall_number[i]), data=branch)
            i_nonrep += 1
    
    hf.close()

    return None



def find_bound_groups(c, loaddir, datadir, outdir):
    """ Take the full list of all infalling objects (regardless of 
    mass, etc.) and find which objects are bound to them at infall. 
    Outputs lists labelled with an index, number of infalls of the 
    host, and the number of members in the group (including the host). 
    List begins with host ID, followed by IDs of any additional bound 
    objects, infalling or not. No restictions on repeated objects. """
    
    infalls = h5py.File(loaddir + 'CLUSTER_%04d_groups.hdf5' % c, 'r')
    keys = np.array(list(infalls.keys()), dtype='str')

    #ids of everything infalling
    infall_ids_t = np.array([np.array(infalls[keys[i]], dtype='int')[0] 
            for i in np.arange(len(keys))])
    
    infall_snaps = infall_ids_t // mod
    infall_ids = infall_ids_t - (infall_snaps*mod+1) #zero index
    #infall_snaps: snapshots at which objects infall
    #infall_snaps_ct: number of objects at each snapshot
    infall_snaps_ct = np.array(list(Counter(infall_snaps).values()))
    infall_snaps = np.array(list(Counter(infall_snaps).keys()))

    #loading relevant data for finding boundness
    xs = h5py.File(datadir + 'xs/CLUSTER_%04d_xs' % c, 'r')
    ys = h5py.File(datadir + 'ys/CLUSTER_%04d_ys' % c, 'r')
    zs = h5py.File(datadir + 'zs/CLUSTER_%04d_zs' % c, 'r')
    vx = h5py.File(datadir + 'vx/CLUSTER_%04d_vx' % c, 'r')
    vy = h5py.File(datadir + 'vy/CLUSTER_%04d_vy' % c, 'r')
    vz = h5py.File(datadir + 'vz/CLUSTER_%04d_vz' % c, 'r')
    ms = h5py.File(datadir + 'ms/CLUSTER_%04d_ms' % c, 'r')
    r200s = h5py.File(datadir + 'rvirs/CLUSTER_%04d_rvirs' % c, 'r')

    counter = 0
    #file to write lists of bound objects to
    bound_out = h5py.File(outdir+'CLUSTER_%04d_bound_members.hdf5' % c, 'w')
    

    for s in range(len(infall_snaps)):
        snap = '%03d' % infall_snaps[s]
        xs_s, ys_s, zs_s, r200s_s = xs[snap], ys[snap], zs[snap], r200s[snap]
        vx_s, vy_s, vz_s = vx[snap], vy[snap], vz[snap]
        ms_s = np.array(ms[snap], dtype='int') 
        ids_s = np.array(np.arange(len(xs_s)), dtype='int')
        
        for i in range(infall_snaps_ct[s]):
            id_n = infall_ids[counter]
            #relative position and velocity of each halo wrt h
            rs_si = ((xs_s-xs_s[id_n])**2. + (ys_s-ys_s[id_n])**2.
                    + (zs_s-zs_s[id_n])**2.)**0.5
            vs_si = ((vx_s-vx_s[id_n])**2. + (vy_s-vy_s[id_n])**2.
                    + (vz_s-vz_s[id_n])**2.)**0.5

            #everything bounded to main infaller
            boundres = bound(vs_si, rs_si, ms_s[id_n], r200s_s[id_n])[0]

            ids_bound = infall_snaps[s]*mod+1+ids_s[boundres]
            #adding bound objects to infaller
            ids_bound = np.append(np.array([infall_ids_t[counter]]), 
                    ids_bound)

            bound_out.create_dataset(keys[counter]+'_%03d' % len(ids_bound), 
                    data=ids_bound)
            counter += 1
    bound_out.close()
    
    return None



def find_member_trees(c, loaddir, mainbranchdir, outdir):
    """ Finds the full branches of each bound member of an infalling 
    group """

    group_objs = h5py.File(loaddir+'CLUSTER_%04d_bound_members.hdf5' % c, 'r')
    infaller_branch = h5py.File(mainbranchdir+'CLUSTER_%04d_groups.hdf5' % c,
            'r')
    
    total_tree = np.array(pd.read_csv('/run/media/ppxrh2/166AA4B87A2DD3B7/'
            'MergerTreeAHF/MergerTreeAHF_ASCII/MergerTree_GadgetX-NewMD'
            'CLUSTER_%04d.txt-CRMratio2' % c, sep='\s+', skiprows=2, 
            usecols=[0], dtype='str')[:-1], dtype='int')[:, 0]

    #keys including id, fate, infall no, no of members
    keys_init = np.array(list(group_objs.keys()))
    if len(keys_init) > 0:
        keys = np.array(np.zeros((len(keys_init), 4)), dtype='str')
        keys[:, [0,1,2]] = np.core.defchararray.partition(keys_init, '_')
        keys[:, [1,2,3]] = np.core.defchararray.partition(
                keys[:, 2], '_')
        keys[:, [2,3]] = np.core.defchararray.partition(
                keys[:, 3], '_')[:, [0, 2]]

    hf = h5py.File(outdir + 'CLUSTER_%04d_grp_memb_branches.hdf5' % c, 'w')
    for i in range(len(keys_init)):
        #first, save branch of host object
        hf.create_dataset(keys_init[i] + '/000_' + keys[i,1], 
                data=infaller_branch[keys[i,0]+'_'+keys[i,1]+'_'+keys[i,2]])
        if keys[i, 3] != '001':
            for j in range(1, int(keys[i, 3])):
                #for each member, find branch and save
                grp_memb_id = np.array(list(group_objs[keys_init[i]]))[j]
                memb_branch = find_branch(np.array([grp_memb_id]), 
                        total_tree, True)
                memb_fate = memb_branch[1]
                memb_branch = memb_branch[0]
                #save in sub-directory, with id number for each member
                hf.create_dataset(keys_init[i] + '/%03d_%01d' % (j, memb_fate),
                        data=memb_branch)

    hf.close()

    return None


def find_object_data(c, loaddir, datadir, outdir):
    """ Finds masses, positions, velocities, etc. of group members 
    after their infall. 

    Note: create_dataset lines are VERY variable in speed. Setting 
    libver='latest' in h5py.File can speed this up, but sacrifice 
    backwards compatibility """
    
    make_files_in_dir(outdir, 
            ['xs','ys','zs','vx','vy','vz','ms','mstar','r200'])
  
    #load trees, and halo data
    trees = h5py.File(loaddir + 'CLUSTER_%04d_grp_memb_branches.hdf5' % c, 'r')
    xs_data = conv_hdf5_to_list(h5py.File(datadir + 'xs/CLUSTER_%04d_xs' % c))
    ys_data = conv_hdf5_to_list(h5py.File(datadir + 'ys/CLUSTER_%04d_ys' % c))
    zs_data = conv_hdf5_to_list(h5py.File(datadir + 'zs/CLUSTER_%04d_zs' % c))
    vx_data = conv_hdf5_to_list(h5py.File(datadir + 'vx/CLUSTER_%04d_vx' % c))
    vy_data = conv_hdf5_to_list(h5py.File(datadir + 'vy/CLUSTER_%04d_vy' % c))
    vz_data = conv_hdf5_to_list(h5py.File(datadir + 'vz/CLUSTER_%04d_vz' % c))
    ms_data = conv_hdf5_to_list(h5py.File(datadir + 'ms/CLUSTER_%04d_ms' % c))
    mstar_data = conv_hdf5_to_list(h5py.File(datadir + 
            'mstars/CLUSTER_%04d_mstars' % c))
    r2_data = conv_hdf5_to_list(h5py.File(datadir + 
            'rvirs/CLUSTER_%04d_rvirs' % c))

    keys = np.array(list(trees.keys()))

    #output files
    hf_xs = h5py.File(outdir + 'xs/CLUSTER_%04d_grp_memb_xs.hdf5' % c, 'w')
    hf_ys = h5py.File(outdir + 'ys/CLUSTER_%04d_grp_memb_ys.hdf5' % c, 'w')
    hf_zs = h5py.File(outdir + 'zs/CLUSTER_%04d_grp_memb_zs.hdf5' % c, 'w')
    hf_vx = h5py.File(outdir + 'vx/CLUSTER_%04d_grp_memb_vx.hdf5' % c, 'w')
    hf_vy = h5py.File(outdir + 'vy/CLUSTER_%04d_grp_memb_vy.hdf5' % c, 'w')
    hf_vz = h5py.File(outdir + 'vz/CLUSTER_%04d_grp_memb_vz.hdf5' % c, 'w')
    hf_ms = h5py.File(outdir + 'ms/CLUSTER_%04d_grp_memb_ms.hdf5' % c, 'w')
    hf_mstar = h5py.File(outdir + 'mstar/CLUSTER_%04d_grp_memb_mstar.hdf5' % c,
            'w')
    hf_r2 = h5py.File(outdir + 'r200/CLUSTER_%04d_grp_memb_r200.hdf5' % c, 'w')

    for i in range(len(keys)):
        key = keys[i]
        keys2 = np.array(list(trees[key].keys()))
        for j in range(len(keys2)):
            key_tot = key + '/' + keys2[j] #key and 'subkey' for halo
            ids = np.array(trees[key_tot])
            snap_ids = ids // mod
            halo_ids = ids - (snap_ids*mod+1)
            l = len(ids)
            xs, ys, zs = np.zeros(l), np.zeros(l), np.zeros(l)
            vx, vy, vz, r2 = np.zeros(l), np.zeros(l), np.zeros(l), np.zeros(l)
            ms = np.array(np.zeros(l), dtype='int')
            mstar = np.array(np.zeros(l), dtype='int')
            for k in range(l):
                #taking data from snapshots/haloes for each grp member 
                xs[k] = xs_data[snap_ids[k]][halo_ids[k]]
                ys[k] = ys_data[snap_ids[k]][halo_ids[k]]
                zs[k] = zs_data[snap_ids[k]][halo_ids[k]]
                vx[k] = vx_data[snap_ids[k]][halo_ids[k]]
                vy[k] = vy_data[snap_ids[k]][halo_ids[k]]
                vz[k] = vz_data[snap_ids[k]][halo_ids[k]]
                ms[k] = ms_data[snap_ids[k]][halo_ids[k]]
                mstar[k] = mstar_data[snap_ids[k]][halo_ids[k]]
                r2[k] = r2_data[snap_ids[k]][halo_ids[k]]
            hf_xs.create_dataset(key_tot, data=xs)
            hf_ys.create_dataset(key_tot, data=ys)
            hf_zs.create_dataset(key_tot, data=zs)
            hf_vx.create_dataset(key_tot, data=vx)
            hf_vy.create_dataset(key_tot, data=vy)
            hf_vz.create_dataset(key_tot, data=vz)
            hf_ms.create_dataset(key_tot, data=ms)
            hf_mstar.create_dataset(key_tot, data=mstar)
            hf_r2.create_dataset(key_tot, data=r2)
            
    hf_xs.close()
    hf_ys.close()
    hf_zs.close()
    hf_vx.close()
    hf_vy.close()
    hf_vz.close()
    hf_ms.close()
    hf_mstar.close()
    hf_r2.close()
   
    return None


def find_cluster_data(clus_range, outdir, oldtrackdir):
    """ Finds masses, positions, velocities, etc. of cluster haloes 
    over their full history. """

    #output files
    hf_id = h5py.File(outdir + 'ALL_CLUSTER_host_id.hdf5', 'w')
    hf_xs = h5py.File(outdir + 'ALL_CLUSTER_host_xs.hdf5', 'w')
    hf_ys = h5py.File(outdir + 'ALL_CLUSTER_host_ys.hdf5', 'w')
    hf_zs = h5py.File(outdir + 'ALL_CLUSTER_host_zs.hdf5', 'w')
    hf_vx = h5py.File(outdir + 'ALL_CLUSTER_host_vx.hdf5', 'w')
    hf_vy = h5py.File(outdir + 'ALL_CLUSTER_host_vy.hdf5', 'w')
    hf_vz = h5py.File(outdir + 'ALL_CLUSTER_host_vz.hdf5', 'w')
    hf_ms = h5py.File(outdir + 'ALL_CLUSTER_host_ms.hdf5', 'w')
    hf_mstar = h5py.File(outdir + 'ALL_CLUSTER_host_mstar.hdf5', 'w')
    hf_r2 = h5py.File(outdir + 'ALL_CLUSTER_host_r200.hdf5', 'w')

    for c in clus_range:
        print(c)
        h_id = c_ids[c-1]
        c_str = '%04d' % c
        bs_track_dir = oldtrackdir + ('NewMDCLUSTER_%04d/snap_128/'
                + 'CLUSTER_%04d') % (c, c)
        ids = ld_arr(bs_track_dir + '.txt', dtype='int')[h_id]
        boo = ids > 0
        ids = ids[boo]
        xs = ld_arr(bs_track_dir + '_xs.txt')[h_id][boo]
        ys = ld_arr(bs_track_dir + '_ys.txt')[h_id][boo]
        zs = ld_arr(bs_track_dir + '_zs.txt')[h_id][boo]
        vx = ld_arr(bs_track_dir + '_vx.txt')[h_id][boo]
        vy = ld_arr(bs_track_dir + '_vy.txt')[h_id][boo]
        vz = ld_arr(bs_track_dir + '_vz.txt')[h_id][boo]
        ms = ld_arr(bs_track_dir + '_ms.txt', dtype='int')[h_id][boo]
        mstar = ld_arr(bs_track_dir + '_mstars.txt', dtype='int')[h_id][boo]
        r2 = ld_arr(bs_track_dir + '_rvirs.txt')[h_id][boo]

        hf_id.create_dataset(c_str, data=ids)
        hf_xs.create_dataset(c_str, data=xs)
        hf_ys.create_dataset(c_str, data=ys)
        hf_zs.create_dataset(c_str, data=zs)
        hf_vx.create_dataset(c_str, data=vx)
        hf_vy.create_dataset(c_str, data=vy)
        hf_vz.create_dataset(c_str, data=vz)
        hf_ms.create_dataset(c_str, data=ms)
        hf_mstar.create_dataset(c_str, data=mstar)
        hf_r2.create_dataset(c_str, data=r2)

    hf_id.close()
    hf_xs.close()
    hf_ys.close()
    hf_zs.close()
    hf_vx.close()
    hf_vy.close()
    hf_vz.close()
    hf_ms.close()
    hf_mstar.close()
    hf_r2.close()

    return None




def find_clust_relative_positions(c, loaddir, membids, cdatadir, outdir):
    """ Find the position, velocity of all infalling objects, relative 
    to the cluster, and their masses at the applicable snapshots """

    make_files_in_dir(outdir, 
            ['xs','ys','zs','vx','vy','vz','ms','mstar','r200'])

    clus_idin = h5py.File(cdatadir + 'ALL_CLUSTER_host_id.hdf5', 'r')['%04d'%c]
    clus_xsin = h5py.File(cdatadir + 'ALL_CLUSTER_host_xs.hdf5', 'r')['%04d'%c]
    clus_ysin = h5py.File(cdatadir + 'ALL_CLUSTER_host_ys.hdf5', 'r')['%04d'%c]
    clus_zsin = h5py.File(cdatadir + 'ALL_CLUSTER_host_zs.hdf5', 'r')['%04d'%c]
    clus_r2in = h5py.File(cdatadir+'ALL_CLUSTER_host_r200.hdf5', 'r')['%04d'%c]
    clus_vxin = h5py.File(cdatadir + 'ALL_CLUSTER_host_vx.hdf5', 'r')['%04d'%c]
    clus_vyin = h5py.File(cdatadir + 'ALL_CLUSTER_host_vy.hdf5', 'r')['%04d'%c]
    clus_vzin = h5py.File(cdatadir + 'ALL_CLUSTER_host_vz.hdf5', 'r')['%04d'%c]

    clus_id = np.array(clus_idin)//mod
    clus_xs, clus_ys = np.zeros(129), np.zeros(129)
    clus_zs, clus_r2 = np.zeros(129), np.zeros(129)
    clus_vx, clus_vy = np.zeros(129), np.zeros(129)
    clus_vz = np.zeros(129)


    clus_xs[clus_id] = np.array(clus_xsin)
    clus_ys[clus_id] = np.array(clus_ysin)
    clus_zs[clus_id] = np.array(clus_zsin)
    clus_r2[clus_id] = np.array(clus_r2in)
    clus_vx[clus_id] = np.array(clus_vxin)
    clus_vy[clus_id] = np.array(clus_vyin)
    clus_vz[clus_id] = np.array(clus_vzin)

    halo_id = h5py.File(membids + 'CLUSTER_%04d_grp_memb_branches.hdf5'%c, 'r')
    halo_xs = h5py.File(loaddir + 'xs/CLUSTER_%04d_grp_memb_xs.hdf5' % c, 'r')
    halo_ys = h5py.File(loaddir + 'ys/CLUSTER_%04d_grp_memb_ys.hdf5' % c, 'r')
    halo_zs = h5py.File(loaddir + 'zs/CLUSTER_%04d_grp_memb_zs.hdf5' % c, 'r')
    halo_vx = h5py.File(loaddir + 'vx/CLUSTER_%04d_grp_memb_vx.hdf5' % c, 'r')
    halo_vy = h5py.File(loaddir + 'vy/CLUSTER_%04d_grp_memb_vy.hdf5' % c, 'r')
    halo_vz = h5py.File(loaddir + 'vz/CLUSTER_%04d_grp_memb_vz.hdf5' % c, 'r')
    halo_ms = h5py.File(loaddir + 'ms/CLUSTER_%04d_grp_memb_ms.hdf5' % c, 'r')
    halo_mstar = h5py.File(loaddir + 
            'mstar/CLUSTER_%04d_grp_memb_mstar.hdf5' % c, 'r')
    halo_r2 = h5py.File(loaddir + 
            'r200/CLUSTER_%04d_grp_memb_r200.hdf5' % c, 'r')

    keys = np.array(list(halo_id.keys()))


    hf_xs = h5py.File(outdir + 'xs/CLUSTER_%04d_xs_reltoCLUS.hdf5'%c, 'w')
    hf_ys = h5py.File(outdir + 'ys/CLUSTER_%04d_ys_reltoCLUS.hdf5'%c, 'w')
    hf_zs = h5py.File(outdir + 'zs/CLUSTER_%04d_zs_reltoCLUS.hdf5'%c, 'w')
    hf_vx = h5py.File(outdir + 'vx/CLUSTER_%04d_vx_reltoCLUS.hdf5'%c, 'w')
    hf_vy = h5py.File(outdir + 'vy/CLUSTER_%04d_vy_reltoCLUS.hdf5'%c, 'w')
    hf_vz = h5py.File(outdir + 'vz/CLUSTER_%04d_vz_reltoCLUS.hdf5'%c, 'w')
    hf_ms = h5py.File(outdir + 'ms/CLUSTER_%04d_ms_reltoCLUS.hdf5'%c, 'w')
    hf_mstar = h5py.File(outdir + 
            'mstar/CLUSTER_%04d_mstar_reltoCLUS.hdf5'%c, 'w')
    hf_r2 = h5py.File(outdir + 'r200/CLUSTER_%04d_r200_reltoCLUS.hdf5'%c, 'w')

    for i in range(len(keys)):
        key = keys[i]
        keys2 = np.array(list(halo_id[key].keys()))
        for j in range(len(keys2)):
            key_tot = key + '/' + keys2[j] #key and 'subkey' for halo
            ids = np.array(halo_id[key_tot])
            snap_ids = ids // mod

            hal_xs, hal_ys = np.zeros(129), np.zeros(129)
            hal_zs, hal_vx = np.zeros(129), np.zeros(129)
            hal_vy, hal_vz = np.zeros(129), np.zeros(129)
            hal_ms, hal_r2 = np.zeros(129), np.zeros(129)
            hal_mstar = np.zeros(129)
            
            hal_xs[snap_ids] = np.array(halo_xs[key_tot])
            hal_ys[snap_ids] = np.array(halo_ys[key_tot])
            hal_zs[snap_ids] = np.array(halo_zs[key_tot])
            hal_vx[snap_ids] = np.array(halo_vx[key_tot])
            hal_vy[snap_ids] = np.array(halo_vy[key_tot])
            hal_vz[snap_ids] = np.array(halo_vz[key_tot])
            hal_ms[snap_ids] = np.array(halo_ms[key_tot])
            hal_mstar[snap_ids] = np.array(halo_mstar[key_tot])
            hal_r2[snap_ids] = np.array(halo_r2[key_tot])
            
            boo = np.where((hal_xs > 0.) * (clus_xs > 0.))[0]
            hal_xs, hal_ys, hal_zs = hal_xs[boo], hal_ys[boo], hal_zs[boo]
            hal_vx, hal_vy, hal_vz = hal_vx[boo], hal_vy[boo], hal_vz[boo]
            hal_ms, hal_mstar, hal_r2 = hal_ms[boo],hal_mstar[boo],hal_r2[boo]
            clu_xs, clu_ys = clus_xs[boo], clus_ys[boo]
            clu_zs, clu_r2 = clus_zs[boo], clus_r2[boo]
            clu_vx, clu_vy = clus_vx[boo], clus_vy[boo]
            clu_vz = clus_vz[boo]

            rel_xs = (hal_xs-clu_xs) / clu_r2
            rel_ys = (hal_ys-clu_ys) / clu_r2
            rel_zs = (hal_zs-clu_zs) / clu_r2
            rel_vx = (hal_vx-clu_vx)
            rel_vy = (hal_vy-clu_vy)
            rel_vz = (hal_vz-clu_vz)

            hf_xs.create_dataset(key_tot, data=rel_xs)
            hf_ys.create_dataset(key_tot, data=rel_ys)
            hf_zs.create_dataset(key_tot, data=rel_zs)
            hf_vx.create_dataset(key_tot, data=rel_vx)
            hf_vy.create_dataset(key_tot, data=rel_vy)
            hf_vz.create_dataset(key_tot, data=rel_vz)
            hf_ms.create_dataset(key_tot, data=hal_ms)
            hf_mstar.create_dataset(key_tot, data=hal_mstar)
            hf_r2.create_dataset(key_tot, data=hal_r2)


    hf_xs.close()
    hf_ys.close()
    hf_zs.close()
    hf_vx.close()
    hf_vy.close()
    hf_vz.close()
    hf_ms.close()
    hf_mstar.close()
    hf_r2.close()


    return None












def find_branch(id_a, id_list, keep_sing=False):
    """ Finds the full evolution of a halo, starting with the earliest 
    id """
    id = id_a[-1]
    tree_search = np.where(id_list==id)[0]
    final_state = 0 #survives to z=0
    if len(tree_search)==0:
        result = np.array([id]*keep_sing, dtype='int')
        final_state = 1 #no progenitor - should not be the case for these
        #only the case when id_a is '0', for second infall
    elif id//mod == id_list[tree_search[0]-1]//mod:
        result = id_a
        final_state = 2 #absorbed by a larger halo - not main progenitor
    elif len(tree_search)==1 and id//mod > id_list[tree_search[0]+1]//mod:
        result = id_a
        final_state = 1 #lost by halo finder
    else:
        id_a = np.append(id_a, id_list[tree_search[0]-1])
        result, final_state = find_branch(id_a, id_list)
    
    return result, final_state


def bound(vrels, rrel, m, r200, incgroup=False):
    """ Determine whether a list of haloes are bound to their host """
    rrels = rrel
    host = rrels<0.001
    rrels[host] = 0.001
    ke = 0.5 * (1000.*vrels)**2.
    ge = ((1.3271244*10.**20.) * m) / (3.0857*10.**19.)
    vir = (ke - (ge / rrels)) < (ge / (-2.5*r200))
    v_cr = 0.001*(1.2*ge/r200)**0.5
    if incgroup==False:
        vir[host]=False
    return vir, rrels/r200, vrels/v_cr


def conv_hdf5_to_list(hdf_file):
    """ Convert a 2D HDF5 file to a list of 1D arrays """
    keys = np.array(list(hdf_file.keys()))
    out_list = []
    for i in range(len(keys)):
        out_list += [np.array(hdf_file[keys[i]])]
    return out_list

def make_files_in_dir(outdir, dirlist):
    """ Make paths for output data if not exist """
    if not os.path.exists(outdir + dirlist[0]):
        for i in range(len(dirlist)):
            os.mkdir(outdir + dirlist[i])
    return None



crange = range(1, 325)

if run_infall_finder==True:
    print('run_infall_finder')
    for clus_no in crange:
        print(clus_no)
        find_infalling(clus_no, all_infalling_objects)


if run_infaller_branches==True:
    print('run_infaller_branches')
    for clus_no in crange:
        print(clus_no)
        find_infaller_evolution(clus_no, all_infalling_objects, 
                all_infalling_branches)

if run_bounded_at_infall==True:
    print('run_bounded_at_infall')
    for clus_no in crange:
        print(clus_no)
        find_bound_groups(clus_no, all_infalling_branches, halo_data, 
                all_infalling_bounded)        

if run_all_grp_member_branches==True:
    print('run_all_grp_member_branches')
    for clus_no in crange:
        print(clus_no)
        find_member_trees(clus_no, all_infalling_bounded,
                all_infalling_branches, all_grp_member_branches)

if run_all_memb_data==True:
    print('run_all_memb_data')
    for clus_no in crange:
        print(clus_no)
        t0 = time.time()
        find_object_data(clus_no, all_grp_member_branches, halo_data,
                all_grp_member_data)

if run_find_cluster_data==True:
    print(run_find_cluster_data)
    find_cluster_data(crange, cluster_halo_data, backsplash_track_data)

if run_all_memb_data_rel_clus==True:
    print('run_all_memb_data_rel_clus')
    for clus_no in crange:
        print(clus_no)
        t0 = time.time()
        find_clust_relative_positions(clus_no, all_grp_member_data, 
                all_grp_member_branches, cluster_halo_data, 
                all_relative_clus_data)







"""
hx = h5py.File('selected_infaller_xs.hdf5', 'w')
hy = h5py.File('selected_infaller_ys.hdf5', 'w')
hz = h5py.File('selected_infaller_zs.hdf5', 'w')
#hvx = h5py.File('selected_infaller_vx.hdf5', 'w')
#hvy = h5py.File('selected_infaller_vy.hdf5', 'w')
#hvz = h5py.File('selected_infaller_vz.hdf5', 'w')
hr = h5py.File('selected_infaller_rs.hdf5', 'w')
crange = np.array(np.loadtxt('/home/ppxrh2/Documents/test_pollux/TheThreeHundred/playground/rhaggar/G3X_data/G3X_300_selected_sample_257.txt'), dtype='int')
counter = 0
for c_no in crange:
    print(c_no)

    data_xs = h5py.File(all_relative_clus_data + 'xs/CLUSTER_%04d_xs_reltoCLUS.hdf5' % c_no, 'r')
    data_ys = h5py.File(all_relative_clus_data + 'ys/CLUSTER_%04d_ys_reltoCLUS.hdf5' % c_no, 'r')
    data_zs = h5py.File(all_relative_clus_data + 'zs/CLUSTER_%04d_zs_reltoCLUS.hdf5' % c_no, 'r')
    data_vx = h5py.File(all_relative_clus_data + 'vx/CLUSTER_%04d_vx_reltoCLUS.hdf5' % c_no, 'r')
    data_vy = h5py.File(all_relative_clus_data + 'vy/CLUSTER_%04d_vy_reltoCLUS.hdf5' % c_no, 'r')
    data_vz = h5py.File(all_relative_clus_data + 'vz/CLUSTER_%04d_vz_reltoCLUS.hdf5' % c_no, 'r')
    data_ms = h5py.File(all_relative_clus_data + 'ms/CLUSTER_%04d_ms_reltoCLUS.hdf5' % c_no, 'r')
    data_mstar = h5py.File(all_relative_clus_data + 'mstar/CLUSTER_%04d_mstar_reltoCLUS.hdf5' % c_no, 'r')

    keys = np.array(list(data_xs.keys()))
    keys_l = []
    for i in range(len(keys)):
        if len(list(data_xs[keys[i]].keys())) > 1:
            keys_l += [keys[i]]
    keys_l = np.array(keys_l)

    for i in range(len(keys_l)):
        keys_l_i = keys_l[i]

        data_xs_n = data_xs[keys_l_i]
        data_ys_n = data_ys[keys_l_i]
        data_zs_n = data_zs[keys_l_i]
        data_vx_n = data_vx[keys_l_i]
        data_vy_n = data_vy[keys_l_i]
        data_vz_n = data_vz[keys_l_i]
        data_ms_n = data_ms[keys_l_i]
        data_mstar_n = data_mstar[keys_l_i]

        keys = np.array(list(data_xs_n.keys()))
        X = np.array(data_xs_n[keys[0]])[0]
        Y = np.array(data_ys_n[keys[0]])[0]
        Z = np.array(data_zs_n[keys[0]])[0]
        R = (X**2. + Y**2. + Z**2.)**0.5
        theta = np.arctan2(Y, X)
        phi = np.arccos(Z/R)

        M = np.array(data_ms_n[keys[0]])[0]
        MS = np.array(data_mstar_n[keys[0]])[0]
        RAT = MS/M
        mass_bool = (M > 10.**10.5) * (MS > 10.**9.5) * (RAT < 0.3)
        if mass_bool==False:
            keys = np.zeros(0)

        for j in range(len(keys)):
            xs = np.array(data_xs_n[keys[j]])
            ys = np.array(data_ys_n[keys[j]])
            zs = np.array(data_zs_n[keys[j]])
            vx = np.array(data_vx_n[keys[j]])
            vy = np.array(data_vy_n[keys[j]])
            vz = np.array(data_vz_n[keys[j]])
            ms = np.array(data_ms_n[keys[j]])
            mstar = np.array(data_mstar_n[keys[j]])
            rat = mstar / ms
            rs = (xs**2. + ys**2. + zs**2.)**0.5
            
            #ts = np.arctan2(ys, xs) - theta
            #ps = np.arccos(zs/rs) + (np.pi/2) - phi
            #xs_n = (1./R) * (xs*X + ys*Y + zs*Z*np.cos(ts-theta))
            #ys_n = (1./R) * (ys*X - xs*Y + zs*Z*np.sin(ts-theta))

            #xs_n = rs * np.cos(ts)*np.sin(ps)
            #ys_n = rs * np.sin(ts)*np.sin(ps)
            #zs_n = rs * np.cos(ps)

            s_theta, c_theta = np.sin(theta), np.cos(theta)
            s_phi, c_phi = np.sin(phi), np.cos(phi)

            xs_n = (zs*c_phi) + s_phi*((xs*c_theta)+(ys*s_theta))
            ys_n = (ys*c_theta) - (xs*s_theta)
            zs_n = (zs*s_phi) - c_phi*((xs*c_theta)+(ys*s_theta))

            vx_n = (vz*c_phi) + s_phi*((vx*c_theta)+(vy*s_theta))
            vy_n = (vy*c_theta) - (vx*s_theta)
            vz_n = (vz*s_phi) - c_phi*((vx*c_theta)+(vy*s_theta))

            ms_n = copy.deepcopy(ms)
            mstar_n = copy.deepcopy(mstar)
            rat_n = copy.deepcopy(rat)


            rs_cut_lim = 0
            rs_cut = np.where(rs[1:] > rs[:-1])[0]
            if len(rs_cut)>0:
                rs_cut = rs_cut[0]
                rs_n = rs[rs_cut:]
                rs_cut_n = np.where(rs_n[1:] < rs_n[:-1])[0]
                if len(rs_cut_n)>0:
                    rs_cut_lim = rs_cut + rs_cut_n[0] + 1
                    xs_n = xs_n[:rs_cut_lim]
                    ys_n = ys_n[:rs_cut_lim]
                    zs_n = zs_n[:rs_cut_lim]
                    vx_n = vx[:rs_cut_lim]
                    vy_n = vy[:rs_cut_lim]
                    vz_n = vz[:rs_cut_lim]
                    ms_n = ms[:rs_cut_lim]
                    mstar_n = mstar[:rs_cut_lim]
                    rat_n = rat[:rs_cut_lim]
            mass_bool = (ms_n > 10.**10.5) * (mstar_n > 10.**9.5) * (rat_n < 0.3)
            xs_n, ys_n, zs_n = xs_n[mass_bool], ys_n[mass_bool], zs_n[mass_bool]
            rs_n = (ys_n**2. + zs_n**2.)**0.5
            rs_n[ys_n<0.] = -1.*rs_n[ys_n<0.]
            #if j==0:
            #    flip_dir = 1. #1 for keep, -1 for reverse direction
            #    if len(ys_n)>1:
            #        if ys_n[0] > ys_n[1]:
            #        #if vy_n[0] < 0.:
            #            flip_dir = -1.
            #rs_n *= flip_dir


            if len(xs_n) > 1:
                hx.create_dataset('%07d' % counter, data=xs_n)
                hy.create_dataset('%07d' % counter, data=ys_n)
                hz.create_dataset('%07d' % counter, data=zs_n)
                hr.create_dataset('%07d' % counter, data=rs_n)
                counter += 1
print(counter)
hx.close()
hy.close()
hz.close()
hr.close()
"""




plt.figure(figsize=(7, 6))
hx = h5py.File('selected_infaller_xs.hdf5', 'r')
hy = h5py.File('selected_infaller_rs.hdf5', 'r')
#hy = h5py.File('selected_infaller_zs.hdf5', 'r')
keys = np.array(list(hx.keys()))
xs_all, ys_all = np.zeros(0), np.zeros(0)
main=0
for i in range(len(keys)):
    xs_all = np.append(xs_all, np.array(hy[keys[i]]))
    ys_all = np.append(ys_all, -1.*np.array(hx[keys[i]]))
    #plt.plot(-1.*np.array(hx[keys[i]]), hy[keys[i]], alpha=0.5, c='b')
    #if main==0:
    #    plt.plot(-1.*np.array(hx[keys[i]]), hy[keys[i]], alpha=1, linewidth=4., c='b')
    #main=1
plt.xlim(-1.5, 2.)
plt.ylim(-1.5, 1.5)
plt.xlabel(r'$x/R_{200}$')
plt.ylabel(r'$y/R_{200}$')
#plt.tight_layout()
#plt.savefig('example_002.png', dpi=400)
#plt.show()

#xbins = np.linspace(-1.5, 1.5, 30)
#ybins = (np.arange(36)-15)/10.
res = 10
xbins = np.linspace(-1.5, 1.5, 6*res)
ybins = (np.arange((7*res)+1)-(3*res))/(2*res)
xvals = np.zeros((len(ybins)-1, 3))
xvals2 = np.zeros((len(ybins)-1, 3))
xvals3 = np.zeros((len(ybins)-1, 3))
xvals4 = np.zeros((len(ybins)-1, 3))
#xhist = np.zeros((len(xbins)-1, len(ybins)-1))
for i in range(len(ybins)-1):
    where = np.where((ys_all > ybins[i]) * (ys_all < ybins[i+1]))[0]
    xvals[i] = find_stdev(xs_all[where], 0.5)
    xvals2[i] = find_stdev(xs_all[where], 0.2)
    xvals3[i] = find_stdev(xs_all[where], 0.05)
    xvals4[i] = find_stdev(xs_all[where], 0.8)
    #xhist[:, i] = np.histogram(xs_all[where], xbins)[0]
    #xhist[:, i] = np.log10(xhist[:, i] + 1.*(xhist[:, i]<1))

xhist = plt.hist2d(ys_all, xs_all, bins=(len(ybins)-1, len(xbins)-1), range=[[-1.5, 2.], [-1.5, 1.5]])[0]

xhist[:, (len(xbins)//2)-1] *= 2.
##########NNOTE: NEEDS NORMALISING TO VOLUME OF ANNULI: BIAS AGAINST R=0 CURRENTLY#############
plt.clf()
xhist = np.transpose(np.log10(xhist + 1.*(xhist<1)))

xbins += 0.5*(xbins[1]-xbins[0])
ybins += 0.5*(ybins[1]-ybins[0])
#plt.plot(ybins[5:-1], xvals[5:,0], c='k', linewidth=3.)
plt.plot(ybins[res:-1], xvals[res:,0]+xvals[res:,1], c='k', linewidth=2.)
plt.plot(ybins[res:-1], xvals[res:,0]-xvals[res:,2], c='k', linewidth=2.)
plt.plot(ybins[res:-1], xvals2[res:,0]+xvals2[res:,1], c='k', linewidth=1.)
plt.plot(ybins[res:-1], xvals2[res:,0]-xvals2[res:,2], c='k', linewidth=1.)
plt.plot(ybins[res:-1], xvals3[res:,0]+xvals3[res:,1], c='k', linewidth=1., linestyle=':')
plt.plot(ybins[res:-1], xvals3[res:,0]-xvals3[res:,2], c='k', linewidth=1., linestyle=':')
#plt.plot(ybins[5:-1], xvals4[5:,0]+xvals4[5:,1], c='k', linewidth=1.)
#plt.plot(ybins[5:-1], xvals4[5:,0]-xvals4[5:,2], c='k', linewidth=1.)

def find_angles(ys, xs):
    ys[:,1] = ys[:,0]+ys[:,1]
    ys[:,2] = ys[:,0]-ys[:,2]
    rs1 = (xs**2. + ys[:,1]**2.)**0.5
    rs2 = (xs**2. + ys[:,2]**2.)**0.5
    rs1, rs2 = rs1[1:], rs2[1:]
    int1 = np.where(rs1<1.)[0][-1]+1
    int2 = np.where(rs2<1.)[0][-1]+1
    xs1 = xs[int1:int1+2]
    xs2 = xs[int2:int2+2]
    ys1 = ys[int1:int1+2, 1]
    ys2 = ys[int2:int2+2, 2]
    ang1 = find_where_r1(xs1, ys1)
    ang2 = find_where_r1(xs2, ys2)
    return ang1+ang2

def find_where_r1(xs, ys):
    """ Find point on line [x0, x1], [y0, y1] where x^2+y^2=1 """
    m = (ys[1]-ys[0]) / (xs[1]-xs[0])
    c = ys[0] - m*xs[0]
    A = m**2. + 1.
    B = 2.*m*c
    C = c**2. - 1.
    x = ((-1.*B) + (B**2. - 4.*A*C)**0.5) / (2.*A)
    y = (1. - x**2.)**0.5
    angle = np.arctan(y/x)*180./np.pi
    return angle

#print(find_angles(xvals4[5:], ybins[5:-1]))
#print(find_angles(xvals[5:], ybins[5:-1]))
#print(find_angles(xvals2[5:], ybins[5:-1]))
#print(find_angles(xvals3[5:], ybins[5:-1]))




rng = np.arange(101)*2.*np.pi/100.
plt.plot(np.sin(rng), np.cos(rng), linewidth=1., c='k', linestyle='--', zorder=999999)
plt.scatter(0., 0., s=120., marker='x', c='r', zorder=1000000)
plt.scatter(-1.15, 0., s=250., marker='>', c='r', zorder=1000000)
plt.plot([-1.4, -1.2], [0., 0.], linewidth=5., c='r', zorder=1000000)

plt.imshow(xhist, extent=(-1.5, 2., -1.5, 1.5), origin='lower', cmap='viridis')
Y, X = np.meshgrid(xbins[:-1], ybins[:-1])
plt.xlabel(r'$x/R_{200}$')
plt.ylabel(r'$y/R_{200}$')
plt.tight_layout()
#plt.savefig('hist2d_rs_poster.png', dpi=600, facecolor='#333333')
#plt.savefig('hist2d_rs.png', dpi=400)

plt.show()


hx.close()
hy.close()






"""
#3d plotting stuff
#plt.figure(figsize=(7, 6))
fig = plt.figure(figsize=(10,10))
ax = Axes3D(fig)
ax.set_aspect('equal')
ax.xaxis.pane.fill =False
ax.yaxis.pane.fill =False
ax.zaxis.pane.fill =False
ax.grid(False)
ax.set_axis_off()
ax.view_init(50,90)

u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)
ax.plot_wireframe(x, y, z, color='k', alpha=0.3, zorder=1000) 

hx = h5py.File('selected_infaller_xs.hdf5', 'r')
hy = h5py.File('selected_infaller_ys.hdf5', 'r')
hz = h5py.File('selected_infaller_zs.hdf5', 'r')
keys = np.array(list(hx.keys()))
xs_all, ys_all = np.zeros(0), np.zeros(0)
main=0


for i in range(len(keys))[-48:-26]:
    xs_all = np.append(xs_all, np.array(hy[keys[i]]))
    ys_all = np.append(ys_all, -1.*np.array(hx[keys[i]]))
    plt.plot(np.array(hx[keys[i]]), hy[keys[i]], hz[keys[i]], alpha=0.8)#, c='b')
    ax.scatter(np.array(hx[keys[i]])[-1], hy[keys[i]][-1], hz[keys[i]][-1], alpha=0.8, s=30., linewidth=0.)#, c='b')
    if main==0:
        plt.plot(np.array(hx[keys[i]]), hy[keys[i]], hz[keys[i]], alpha=1, linewidth=3., c='b', zorder=999)
        ax.scatter(np.array(hx[keys[i]])[-1], hy[keys[i]][-1], hz[keys[i]][-1], alpha=1, c='b', s=50., zorder=999)
    main=1
#plt.xlim(-1.5, 2.)
#plt.ylim(-1.5, 1.5)
plt.scatter(0., 0., s=100., marker='x', c='r', zorder=999)
plt.scatter(1.1, 0., s=100., marker='>', c='r', zorder=999)
ax.set_xlim(-2.8, 0.7)
ax.set_ylim(-1.5, 2.)
ax.set_zlim(-1.5, 2.)
plt.xlabel(r'$x/R_{200}$')
plt.ylabel(r'$y/R_{200}$')
#plt.tight_layout()
plt.savefig('example_3d.png', dpi=400)
plt.show()
print(alphabet)
"""