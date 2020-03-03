from modules import *



reduced_cinfo = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/Infalling_'
        + 'Groups/MergerTreeAHF_Infalling_Re-written/NewMDCLUSTER_0002_reruns'
        + '/reduced_cluster_info/xs/CLUSTER_%04d_xs')

mainproginfo = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/Infalling_'
        + 'Groups/MergerTreeAHF_Infalling_Re-written/NewMDCLUSTER_0002_reruns'
        + '/main_progs/CLUSTER_%04d_CRMratio2_MainProg/CLUSTER_%04d_prog'
        + '_output.snapshot_%03d.AHF.MainProg.txt')

savedir = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/Infalling_Groups'
        + '/MergerTreeAHF_Infalling_Re-written/NewMDCLUSTER_0002_reruns/'
        + 'full_tracking/NewMDCLUSTER_%04d/snap_128/')


def find_prog_tree(c, present_snap, redcinfo, mainpinfo, outdir, 
        save=False):
    """ Find the progenitors of all haloes in a given snapshot, for a 
    particular cluster. Saves the tree, and rvir, x, y, z values """
    
    outdir = outdir % c
    if not os.path.exists(outdir[:-10]):
        os.mkdir(outdir[:-10])
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    lookback_snaps = present_snap+20 - np.arange(20, present_snap)

    hr = h5py.File(redcinfo % c)
    hf = hr['%03d' % present_snap]
    length = len(hf)
    final_list = np.arange(1, length+1) + present_snap*mod

    full_tree = np.zeros((length, len(lookback_snaps)+1), dtype='int')
    maxold = len(full_tree)
    full_tree[:, -1] = final_list #merger tree array
    nonzero = np.where(full_tree[:, -1] > 0)[0]
    ns = full_tree[:, -1] - (present_snap*mod+1) #indeces of present haloes
    

    if ((c == 10 and present_snap == 100) or (c == 228 and 
            present_snap == 110)): #exceptions for missing snaps
        full_tree = np.zeros((1, len(lookback_snaps)+1), dtype='int')
    else:
        for i in lookback_snaps:
            tree_i = ld_arr(mainpinfo % (c, c, i), dtype='int')

            #exception for missing snapshots:
            if ((c == 10 and i == 100) or (c == 228 and i == 110)):
                tree_i = np.asarray(np.zeros((2, 2)), dtype='int')

            max = np.max([tree_i[-1, 0]-(i*mod), len(full_tree), maxold, 
                    np.max(full_tree[:, i-20]) - (i*mod)]) #safe length limit

            maxold = tree_i[-1, 0]-(i*mod)
            ext_tree = np.array(np.zeros((max, 2)), dtype='int')

            tree_i_ids = tree_i[:, 0] - (i*mod+1)
            #or for missing snapshots:
            if ((c == 10 and i == 100) or (c == 228 and i == 110)):
                tree_i_ids = tree_i[:, 0]
            ext_tree[tree_i_ids] = tree_i

            #snapshots of progenitor haloes, and writing to merger tree
            snapnews = np.asarray(ext_tree[ns, 1], dtype='int') // mod
            sn_mask = snapnews>0
            full_tree[nonzero[sn_mask], snapnews[sn_mask]-20] = (
                    ext_tree[ns, 1][sn_mask])
    
            nonzero = np.where(full_tree[:, i-21] > 0)[0]
            ns = full_tree[nonzero, i-21] - ((i-1)*mod+1)

    if save==True: #save results

        full_tree_str = np.char.mod('%15d', full_tree)
        headers = np.char.mod('%15d', np.arange(20, present_snap+1))
        pd.DataFrame(full_tree_str).to_csv(outdir+'CLUSTER_%04d.txt' % c,
                header=headers, index=None, sep='\t')
	
    return full_tree


find_prog_tree(2000, 128, reduced_cinfo, mainproginfo, savedir, True)
find_prog_tree(2001, 128, reduced_cinfo, mainproginfo, savedir, True)


