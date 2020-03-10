from modules import *


run_infall_finder = False
run_infaller_branches = False
run_bounded_at_infall = False
run_all_grp_member_branches = False
run_all_memb_data = False
run_find_cluster_data = False
run_all_memb_data_rel_clus = False
run_all_memb_data_rel_grp = False

"""
crange = [2]#range(1, 325)

base_folder = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'Infalling_Groups/MergerTreeAHF_Infalling_Re-written/')
halo_data = ('/run/media/ppxrh2/166AA4B87A2DD3B7/NewMDCLUSTER_data/reduced_'
        + 'cluster_info/')
cluster_halo_data = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/Infall'
        + 'ing_Groups/MergerTreeAHF_Infalling_Re-written/cluster_halo_data/')
backsplash_track_data = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'MergerTreeAHF_General_Tree_Comp/')
main_prog_data = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'MergerTreeAHF_HDF5_MainProg/')
ascii_prog_data = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
        + 'MergerTreeAHF_ASCII/')

"""

base_folder = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/Infalling_Gr'
        + 'oups/MergerTreeAHF_Infalling_Re-written/NewMDCLUSTER_0002_reruns/')
halo_data = (base_folder + 'reduced_cluster_info/')
cluster_halo_data = (base_folder + 'cluster_halo_data/')
backsplash_track_data = (base_folder + 'full_tracking/')
main_prog_data = (base_folder + 'main_progs/')
ascii_prog_data = base_folder
#"""

all_infalling_objects = base_folder + 'all_infalling_objects/'
all_infalling_branches = base_folder + 'all_infaller_branches/'
all_infalling_bounded = base_folder + 'all_bounded_at_infall_ids/'
all_grp_member_branches = base_folder + 'all_group_member_branches/'
all_grp_member_data = base_folder + 'all_members_halo_data/'
all_relative_clus_data = base_folder + 'all_members_halo_data/'
all_relative_grp_data = base_folder + 'all_members_halo_data/'



def find_infalling(c, outdir, loaddir):
    """ Find all the objects which pass within R200 of the cluster, 
    having been outside of R200 at a previous snapshot. Does not 
    account for double infalls, sub-subhaloes, etc. """
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    xs_in  = h5py.File(loaddir+'xs/CLUSTER_%04d_xs' % c)
    ys_in  = h5py.File(loaddir+'ys/CLUSTER_%04d_ys' % c)
    zs_in  = h5py.File(loaddir+'zs/CLUSTER_%04d_zs' % c)
    r200_in = h5py.File(loaddir+'rvirs/CLUSTER_%04d_rvirs' % c)

    halo_ids = ld_arr(backsplash_track_data + ('NewMDCLUSTER_%04d/snap_128/'
            + 'CLUSTER_%04d.txt') % (c, c), dtype='int')
    halo_ids_c = halo_ids[c_ids[c-1]]
    halo_ids_c = halo_ids_c[halo_ids_c > 0]
    keys = np.char.mod(u'%03d', halo_ids_c // mod)


    j = keys[0]
    xs = np.array(xs_in[j])
    ys = np.array(ys_in[j])
    zs = np.array(zs_in[j])
    ids = np.array(np.arange(len(xs)), dtype='int')+int(j)*mod+1
    r200s = np.array(r200_in[j])
    
    h_id = halo_ids_c[0] - (int(j)*mod+1)
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
        ids_main_p_in = ld_arr(main_prog_data + ('CLUSTER_%04d_CRMratio2_Main'
                + 'Prog/CLUSTER_%04d_prog_output.snapshot_' + j 
                + '.AHF.MainProg.txt') % (c, c), dtype='int')
        j = int(j)
        ids_main_p = np.array(np.zeros(len(xs)), dtype='int')
        ids_main_p[ids_main_p_in[:, 0]-(j*mod+1)] = ids_main_p_in[:, 1]
        
        h_id = halo_ids_c[i] - (int(j)*mod+1)
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

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    total_tree = np.array(pd.read_csv(ascii_prog_data + ('MergerTree_GadgetX-'
            + 'NewMDCLUSTER_%04d.txt-CRMratio2') % c, sep='\s+', skiprows=2, 
            usecols=[0], dtype='str')[:-1], dtype='int')[:, 0]
    
    ids_infall_t = ld_arr(indir + 'CLUSTER_%04d_all_infallers.txt' % c, 
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
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)

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

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    group_objs = h5py.File(loaddir+'CLUSTER_%04d_bound_members.hdf5' % c, 'r')
    infaller_branch = h5py.File(mainbranchdir+'CLUSTER_%04d_groups.hdf5' % c,
            'r')
    
    total_tree = np.array(pd.read_csv(ascii_prog_data + ('MergerTree_GadgetX-'
            + 'NewMDCLUSTER_%04d.txt-CRMratio2') % c, sep='\s+', skiprows=2, 
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
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outdir += 'absolute/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    make_files_in_dir(outdir, 
            ['xs','ys','zs','vx','vy','vz','ms','mstars','mgas','r200'])
  
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
    mgas_data = conv_hdf5_to_list(h5py.File(datadir + 
            'mgas/CLUSTER_%04d_mgas' % c))
    r2_data = conv_hdf5_to_list(h5py.File(datadir + 
            'r200/CLUSTER_%04d_r200' % c))

    keys = np.array(list(trees.keys()))

    #output files
    hf_xs = h5py.File(outdir + 'xs/CLUSTER_%04d_grp_memb_xs.hdf5' % c, 'w')
    hf_ys = h5py.File(outdir + 'ys/CLUSTER_%04d_grp_memb_ys.hdf5' % c, 'w')
    hf_zs = h5py.File(outdir + 'zs/CLUSTER_%04d_grp_memb_zs.hdf5' % c, 'w')
    hf_vx = h5py.File(outdir + 'vx/CLUSTER_%04d_grp_memb_vx.hdf5' % c, 'w')
    hf_vy = h5py.File(outdir + 'vy/CLUSTER_%04d_grp_memb_vy.hdf5' % c, 'w')
    hf_vz = h5py.File(outdir + 'vz/CLUSTER_%04d_grp_memb_vz.hdf5' % c, 'w')
    hf_ms = h5py.File(outdir + 'ms/CLUSTER_%04d_grp_memb_ms.hdf5' % c, 'w')
    hf_mstar = h5py.File(outdir + 
            'mstars/CLUSTER_%04d_grp_memb_mstars.hdf5' % c, 'w')
    hf_mgas = h5py.File(outdir + 'mgas/CLUSTER_%04d_grp_memb_mgas.hdf5'%c, 'w')
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
            mgas = np.array(np.zeros(l), dtype='int')
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
                mgas[k] = mgas_data[snap_ids[k]][halo_ids[k]]
                r2[k] = r2_data[snap_ids[k]][halo_ids[k]]
            hf_xs.create_dataset(key_tot, data=xs)
            hf_ys.create_dataset(key_tot, data=ys)
            hf_zs.create_dataset(key_tot, data=zs)
            hf_vx.create_dataset(key_tot, data=vx)
            hf_vy.create_dataset(key_tot, data=vy)
            hf_vz.create_dataset(key_tot, data=vz)
            hf_ms.create_dataset(key_tot, data=ms)
            hf_mstar.create_dataset(key_tot, data=mstar)
            hf_mgas.create_dataset(key_tot, data=mgas)
            hf_r2.create_dataset(key_tot, data=r2)
            
    hf_xs.close()
    hf_ys.close()
    hf_zs.close()
    hf_vx.close()
    hf_vy.close()
    hf_vz.close()
    hf_ms.close()
    hf_mstar.close()
    hf_mgas.close()
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
        r2 = ld_arr(bs_track_dir + '_r200.txt')[h_id][boo]

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

    loaddir += 'absolute/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outdir += 'rel_to_cluster/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    make_files_in_dir(outdir, 
            ['id','xs','ys','zs','vx','vy','vz','ms','mstars','mgas','r200'])

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
            'mstars/CLUSTER_%04d_grp_memb_mstars.hdf5' % c, 'r')
    halo_mgas = h5py.File(loaddir + 
            'mgas/CLUSTER_%04d_grp_memb_mgas.hdf5' % c, 'r')
    halo_r2 = h5py.File(loaddir + 
            'r200/CLUSTER_%04d_grp_memb_r200.hdf5' % c, 'r')

    keys = np.array(list(halo_id.keys()))


    hf_id = h5py.File(outdir + 'id/CLUSTER_%04d_id_reltoCLUS.hdf5'%c, 'w')
    hf_xs = h5py.File(outdir + 'xs/CLUSTER_%04d_xs_reltoCLUS.hdf5'%c, 'w')
    hf_ys = h5py.File(outdir + 'ys/CLUSTER_%04d_ys_reltoCLUS.hdf5'%c, 'w')
    hf_zs = h5py.File(outdir + 'zs/CLUSTER_%04d_zs_reltoCLUS.hdf5'%c, 'w')
    hf_vx = h5py.File(outdir + 'vx/CLUSTER_%04d_vx_reltoCLUS.hdf5'%c, 'w')
    hf_vy = h5py.File(outdir + 'vy/CLUSTER_%04d_vy_reltoCLUS.hdf5'%c, 'w')
    hf_vz = h5py.File(outdir + 'vz/CLUSTER_%04d_vz_reltoCLUS.hdf5'%c, 'w')
    hf_ms = h5py.File(outdir + 'ms/CLUSTER_%04d_ms_reltoCLUS.hdf5'%c, 'w')
    hf_mstar = h5py.File(outdir + 
            'mstars/CLUSTER_%04d_mstars_reltoCLUS.hdf5'%c, 'w')
    hf_mgas = h5py.File(outdir + 
            'mgas/CLUSTER_%04d_mgas_reltoCLUS.hdf5'%c, 'w')
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
            hal_r2 = np.zeros(129)
            hal_ms = np.array(np.zeros(129), dtype='int')
            hal_mstar = np.array(np.zeros(129), dtype='int')
            hal_mgas = np.array(np.zeros(129), dtype='int')
            hal_ids = np.array(np.zeros(129), dtype='int')
            
            hal_xs[snap_ids] = np.array(halo_xs[key_tot])
            hal_ys[snap_ids] = np.array(halo_ys[key_tot])
            hal_zs[snap_ids] = np.array(halo_zs[key_tot])
            hal_vx[snap_ids] = np.array(halo_vx[key_tot])
            hal_vy[snap_ids] = np.array(halo_vy[key_tot])
            hal_vz[snap_ids] = np.array(halo_vz[key_tot])
            hal_ms[snap_ids] = np.array(halo_ms[key_tot])
            hal_mstar[snap_ids] = np.array(halo_mstar[key_tot])
            hal_mgas[snap_ids] = np.array(halo_mgas[key_tot])
            hal_r2[snap_ids] = np.array(halo_r2[key_tot])
            hal_ids[snap_ids] = ids
            
            boo = np.where((hal_xs > 0.) * (clus_xs > 0.))[0]
            hal_xs, hal_ys, hal_zs = hal_xs[boo], hal_ys[boo], hal_zs[boo]
            hal_vx, hal_vy, hal_vz = hal_vx[boo], hal_vy[boo], hal_vz[boo]
            hal_ms, hal_mstar, hal_r2 = hal_ms[boo],hal_mstar[boo],hal_r2[boo]
            hal_mgas, ids = hal_mgas[boo], hal_ids[boo]
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
            hf_mgas.create_dataset(key_tot, data=hal_mgas)
            hf_r2.create_dataset(key_tot, data=hal_r2)
            hf_id.create_dataset(key_tot, data=ids)


    hf_xs.close()
    hf_ys.close()
    hf_zs.close()
    hf_vx.close()
    hf_vy.close()
    hf_vz.close()
    hf_ms.close()
    hf_mstar.close()
    hf_mgas.close()
    hf_r2.close()
    hf_id.close()


    return None



def find_grp_relative_positions(c, loaddir, membids, cdatadir, outdir):
    """ Find the position, velocity of all infalling objects, relative 
    to their host group, and their masses at the applicable snapshots """

    loaddir += 'absolute/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outdir += 'rel_to_grp/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    make_files_in_dir(outdir, 
            ['id','xs','ys','zs','vx','vy','vz','ms','mstars','mgas','r200'])

    clus_idin = h5py.File(cdatadir + 'ALL_CLUSTER_host_id.hdf5', 'r')['%04d'%c]
    clus_idin = np.array(clus_idin)//mod
    clus_id = np.zeros(129)
    clus_id[clus_idin] = np.array(clus_idin)

    halo_id = h5py.File(membids + 'CLUSTER_%04d_grp_memb_branches.hdf5'%c, 'r')
    halo_xs = h5py.File(loaddir + 'xs/CLUSTER_%04d_grp_memb_xs.hdf5' % c, 'r')
    halo_ys = h5py.File(loaddir + 'ys/CLUSTER_%04d_grp_memb_ys.hdf5' % c, 'r')
    halo_zs = h5py.File(loaddir + 'zs/CLUSTER_%04d_grp_memb_zs.hdf5' % c, 'r')
    halo_vx = h5py.File(loaddir + 'vx/CLUSTER_%04d_grp_memb_vx.hdf5' % c, 'r')
    halo_vy = h5py.File(loaddir + 'vy/CLUSTER_%04d_grp_memb_vy.hdf5' % c, 'r')
    halo_vz = h5py.File(loaddir + 'vz/CLUSTER_%04d_grp_memb_vz.hdf5' % c, 'r')
    halo_ms = h5py.File(loaddir + 'ms/CLUSTER_%04d_grp_memb_ms.hdf5' % c, 'r')
    halo_mstar = h5py.File(loaddir + 
            'mstars/CLUSTER_%04d_grp_memb_mstars.hdf5' % c, 'r')
    halo_mgas = h5py.File(loaddir + 
            'mgas/CLUSTER_%04d_grp_memb_mgas.hdf5' % c, 'r')
    halo_r2 = h5py.File(loaddir + 
            'r200/CLUSTER_%04d_grp_memb_r200.hdf5' % c, 'r')

    keys = np.array(list(halo_id.keys()))


    hf_id = h5py.File(outdir + 'id/CLUSTER_%04d_id_reltoGRP.hdf5'%c, 'w')
    hf_xs = h5py.File(outdir + 'xs/CLUSTER_%04d_xs_reltoGRP.hdf5'%c, 'w')
    hf_ys = h5py.File(outdir + 'ys/CLUSTER_%04d_ys_reltoGRP.hdf5'%c, 'w')
    hf_zs = h5py.File(outdir + 'zs/CLUSTER_%04d_zs_reltoGRP.hdf5'%c, 'w')
    hf_vx = h5py.File(outdir + 'vx/CLUSTER_%04d_vx_reltoGRP.hdf5'%c, 'w')
    hf_vy = h5py.File(outdir + 'vy/CLUSTER_%04d_vy_reltoGRP.hdf5'%c, 'w')
    hf_vz = h5py.File(outdir + 'vz/CLUSTER_%04d_vz_reltoGRP.hdf5'%c, 'w')
    hf_ms = h5py.File(outdir + 'ms/CLUSTER_%04d_ms_reltoGRP.hdf5'%c, 'w')
    hf_mstar = h5py.File(outdir + 
            'mstars/CLUSTER_%04d_mstars_reltoGRP.hdf5'%c, 'w')
    hf_mgas = h5py.File(outdir + 
            'mgas/CLUSTER_%04d_mgas_reltoGRP.hdf5'%c, 'w')
    hf_r2 = h5py.File(outdir + 'r200/CLUSTER_%04d_r200_reltoGRP.hdf5'%c, 'w')

    for i in range(len(keys)):
        key = keys[i]
        keys2 = np.array(list(halo_id[key].keys()))

        grp_idin = np.array(halo_id[key][keys2[0]])//mod
        grp_xs, grp_ys = np.zeros(129), np.zeros(129)
        grp_zs, grp_r2 = np.zeros(129), np.zeros(129)
        grp_vx, grp_vy = np.zeros(129), np.zeros(129)
        grp_vz = np.zeros(129)
        grp_id = np.array(np.zeros(129), dtype='int')

        grp_xs[grp_idin] = np.array(halo_xs[key][keys2[0]])
        grp_ys[grp_idin] = np.array(halo_ys[key][keys2[0]])
        grp_zs[grp_idin] = np.array(halo_zs[key][keys2[0]])
        grp_vx[grp_idin] = np.array(halo_vx[key][keys2[0]])
        grp_vy[grp_idin] = np.array(halo_vy[key][keys2[0]])
        grp_vz[grp_idin] = np.array(halo_vz[key][keys2[0]])
        grp_r2[grp_idin] = np.array(halo_r2[key][keys2[0]])
        grp_id[grp_idin] = np.array(halo_id[key][keys2[0]])

        for j in range(len(keys2)):
            key_tot = key + '/' + keys2[j] #key and 'subkey' for halo
            ids = np.array(halo_id[key_tot])
            snap_ids = ids // mod

            hal_xs, hal_ys = np.zeros(129), np.zeros(129)
            hal_zs, hal_vx = np.zeros(129), np.zeros(129)
            hal_vy, hal_vz = np.zeros(129), np.zeros(129)
            hal_r2 = np.zeros(129)
            hal_ms = np.array(np.zeros(129), dtype='int')
            hal_mstar = np.array(np.zeros(129), dtype='int')
            hal_mgas = np.array(np.zeros(129), dtype='int')
            hal_ids = np.array(np.zeros(129), dtype='int')
            
            hal_xs[snap_ids] = np.array(halo_xs[key_tot])
            hal_ys[snap_ids] = np.array(halo_ys[key_tot])
            hal_zs[snap_ids] = np.array(halo_zs[key_tot])
            hal_vx[snap_ids] = np.array(halo_vx[key_tot])
            hal_vy[snap_ids] = np.array(halo_vy[key_tot])
            hal_vz[snap_ids] = np.array(halo_vz[key_tot])
            hal_ms[snap_ids] = np.array(halo_ms[key_tot])
            hal_mstar[snap_ids] = np.array(halo_mstar[key_tot])
            hal_mgas[snap_ids] = np.array(halo_mgas[key_tot])
            hal_r2[snap_ids] = np.array(halo_r2[key_tot])
            hal_ids[snap_ids] = ids
            if (grp_id[grp_id>0][-1]//mod < hal_ids[hal_ids>0][-1]//mod)==True:
                key_tot = key_tot[:-1]+'3'
                #label haloes where group is lost before halo as '3'
            
            boo = np.where((hal_xs > 0.) * (grp_id > 0))[0]
            hal_xs, hal_ys, hal_zs = hal_xs[boo], hal_ys[boo], hal_zs[boo]
            hal_vx, hal_vy, hal_vz = hal_vx[boo], hal_vy[boo], hal_vz[boo]
            hal_ms, hal_mstar, hal_r2 = hal_ms[boo],hal_mstar[boo],hal_r2[boo]
            hal_mgas, ids = hal_mgas[boo], hal_ids[boo]
            gr_xs, gr_ys = grp_xs[boo], grp_ys[boo]
            gr_zs, gr_r2 = grp_zs[boo], grp_r2[boo]
            gr_vx, gr_vy = grp_vx[boo], grp_vy[boo]
            gr_vz = grp_vz[boo]

            rel_xs = (hal_xs-gr_xs) / gr_r2
            rel_ys = (hal_ys-gr_ys) / gr_r2
            rel_zs = (hal_zs-gr_zs) / gr_r2
            rel_vx = (hal_vx-gr_vx)
            rel_vy = (hal_vy-gr_vy)
            rel_vz = (hal_vz-gr_vz)

            hf_xs.create_dataset(key_tot, data=rel_xs)
            hf_ys.create_dataset(key_tot, data=rel_ys)
            hf_zs.create_dataset(key_tot, data=rel_zs)
            hf_vx.create_dataset(key_tot, data=rel_vx)
            hf_vy.create_dataset(key_tot, data=rel_vy)
            hf_vz.create_dataset(key_tot, data=rel_vz)
            hf_ms.create_dataset(key_tot, data=hal_ms)
            hf_mstar.create_dataset(key_tot, data=hal_mstar)
            hf_mgas.create_dataset(key_tot, data=hal_mgas)
            hf_r2.create_dataset(key_tot, data=hal_r2)
            hf_id.create_dataset(key_tot, data=ids)


    hf_xs.close()
    hf_ys.close()
    hf_zs.close()
    hf_vx.close()
    hf_vy.close()
    hf_vz.close()
    hf_ms.close()
    hf_mstar.close()
    hf_mgas.close()
    hf_r2.close()
    hf_id.close()


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

def v_crit(m200, r200):
    ge = ((1.3271244*10.**20.) * m200) / (3.0857*10.**19.)
    v_cr = 0.001*(1.2*ge/r200)**0.5
    return v_cr


def conv_hdf5_to_list(hdf_file):
    """ Convert a 2D HDF5 file to a list of 1D arrays """
    keys = np.array(list(hdf_file.keys()))
    out_list = []
    for i in range(len(keys)):
        out_list += [np.array(hdf_file[keys[i]])]
    return out_list

def make_files_in_dir(outdir, dirlist):
    """ Make paths for output data if not exist """
    for i in range(len(dirlist)):
        if not os.path.exists(outdir + dirlist[i]):
            os.mkdir(outdir + dirlist[i])
    return None




###FOR CLUSTER 2 re-runs:
#c_ids_new = np.array(np.zeros(2002), dtype='int')
#c_ids_new[crange] = c_ids[1]
#c_ids = c_ids_new



if run_infall_finder==True:
    print('run_infall_finder')
    for clus_no in crange:
        print(clus_no)
        find_infalling(clus_no, all_infalling_objects, halo_data)


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
    print('run_find_cluster_data')
    find_cluster_data(crange, cluster_halo_data, backsplash_track_data)

if run_all_memb_data_rel_clus==True:
    print('run_all_memb_data_rel_clus')
    for clus_no in crange:
        print(clus_no)
        t0 = time.time()
        find_clust_relative_positions(clus_no, all_grp_member_data, 
                all_grp_member_branches, cluster_halo_data, 
                all_relative_clus_data)

if run_all_memb_data_rel_grp==True:
    print('run_all_memb_data_rel_grp')
    for clus_no in crange:
        print(clus_no)
        t0 = time.time()
        find_grp_relative_positions(clus_no, all_grp_member_data, 
                all_grp_member_branches, cluster_halo_data, 
                all_relative_grp_data)




c_no = 2000
rs = np.zeros(0)
vs = np.zeros(0)

id = h5py.File(all_grp_member_data + 'rel_to_grp/id/CLUSTER_%04d_id_reltoGRP.hdf5' % c_no, 'r')
id = h5py.File(all_grp_member_branches + 'CLUSTER_%04d_grp_memb_branches.hdf5' % c_no, 'r')
xs = h5py.File(all_grp_member_data + 'rel_to_grp/xs/CLUSTER_%04d_xs_reltoGRP.hdf5' % c_no, 'r')
ys = h5py.File(all_grp_member_data + 'rel_to_grp/ys/CLUSTER_%04d_ys_reltoGRP.hdf5' % c_no, 'r')
zs = h5py.File(all_grp_member_data + 'rel_to_grp/zs/CLUSTER_%04d_zs_reltoGRP.hdf5' % c_no, 'r')
vx = h5py.File(all_grp_member_data + 'rel_to_grp/vx/CLUSTER_%04d_vx_reltoGRP.hdf5' % c_no, 'r')
vy = h5py.File(all_grp_member_data + 'rel_to_grp/vy/CLUSTER_%04d_vy_reltoGRP.hdf5' % c_no, 'r')
vz = h5py.File(all_grp_member_data + 'rel_to_grp/vz/CLUSTER_%04d_vz_reltoGRP.hdf5' % c_no, 'r')
ms = h5py.File(all_grp_member_data + 'rel_to_grp/ms/CLUSTER_%04d_ms_reltoGRP.hdf5' % c_no, 'r')
mstars = h5py.File(all_grp_member_data + 'rel_to_grp/mstars/CLUSTER_%04d_mstars_reltoGRP.hdf5' % c_no, 'r')
r2 = h5py.File(all_grp_member_data + 'rel_to_grp/r200/CLUSTER_%04d_r200_reltoGRP.hdf5' % c_no, 'r')

keys = np.array(list(id.keys()))
keys1 = keys[-1]
print(keys1)
keys2 = np.array(list(id[keys[-1]].keys()))[-1]
print(np.array(id[keys1][keys2]))



keys = np.array(list(xs.keys()))
for i in range(len(keys)):
    keysi = np.array(list(xs[keys[i]].keys()))
    xsi = xs[keys[i]]
    ysi = ys[keys[i]]
    zsi = zs[keys[i]]
    vxi = vx[keys[i]]
    vyi = vy[keys[i]]
    vzi = vz[keys[i]]
    msi = ms[keys[i]]
    mti = mstars[keys[i]]
    r2i = r2[keys[i]]
    ms_grp = np.array(msi[keysi[0]])[0]
    #mt_grp = np.array(mti[keysi[0]])[0]
    #gal_sel = (ms_grp > 10.**10.5) * (mt_grp > 10.**9.5) * (
    #        (mt_grp/ms_grp)<0.3)
    gal_sel = (ms_grp > 10.**10.5)
    if gal_sel == True:
        #print(keys[i])
        for j in range(1, len(keysi)):
            ms_arr = np.array(msi[keysi[j]])[0]
            #mt_arr = np.array(mti[keysi[j]])[0]
            #gal_sel = (ms_arr > 10.**10.5) * (mt_arr > 10.**9.5) * (
            #        (mt_arr/ms_arr)<0.3) * (ms_grp > ms_arr)
            gal_sel = (ms_arr > 10.**10.5) * (ms_grp > ms_arr)

            if gal_sel == True:
                r = (np.array(xsi[keysi[j]])[0]**2. + 
                np.array(ysi[keysi[j]])[0]**2. + 
                np.array(zsi[keysi[j]])[0]**2.)**0.5

                v = (np.array(vxi[keysi[j]])[0]**2. + 
                np.array(vyi[keysi[j]])[0]**2. + 
                np.array(vzi[keysi[j]])[0]**2.)**0.5
                critical = v_crit(np.array(msi[keysi[0]])[0], 
                        np.array(r2i[keysi[0]])[0])
                rs = np.append(rs, r)
                vs = np.append(vs, v/critical)



def bound_crit(x):
    return ((5./(3.*x)) - (2./3.))**0.5



plt.figure(figsize=(7, 5.5))
plt.plot(np.arange(1,251)/100., bound_crit(np.arange(1,251)/100.), linewidth=2., c='r')
plt.scatter(rs, vs, s=30.)
#plt.title('Original softening length', size=18)
plt.title('DM-only', size=18)
#plt.title('Reduced softening length', size=18)
plt.xlim(0., 3.)
plt.ylim(0., 2.5)
plt.xlabel(r'$r/R_{200,\rm{host}}$')
plt.ylabel(r'$v/v_{\rm{crit}}$')
plt.tight_layout()
#plt.xlim(0., 3.)
#plt.ylim(0., 3.)
#plt.savefig('cluster_0002.png', dpi=200)
plt.savefig('cluster_0002_dm.png', dpi=200)
#plt.savefig('cluster_0002_red-soft.png', dpi=200)
plt.show()