from modules import *

""" Find main progenitors from a large ASCII file """

savedir = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/Infalling_Groups'
        + '/MergerTreeAHF_Infalling_Re-written/NewMDCLUSTER_0002_reruns/'
        + 'main_progs/')
loaddir = ('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/Infalling_Groups'
        + '/MergerTreeAHF_Infalling_Re-written/NewMDCLUSTER_0002_reruns/')


def find_progenitors(c, indir, outdir, snap_rng=range(20, 129)):
    tree = pd.read_csv(indir
            +'MergerTree_GadgetX-NewMDCLUSTER_%04d.txt-CRMratio2' % c)

    outdirc = outdir+'CLUSTER_%04d_CRMratio2_MainProg' % c
    if not os.path.exists(outdirc):
        os.mkdir(outdirc)

    tree = np.array(tree, dtype='str')[2:-1][:, 0]
    tree = np.core.defchararray.partition(tree, ' ')
    jumps = tree[:, 2]
    tree = np.array(tree[:, 0], dtype='int')
    tree_snaps = tree // mod

    for i in range(len(snap_rng)):
        whr = np.where((tree_snaps==snap_rng[-1-i]) * (jumps != ''))[0]
        mainprogs = np.array(np.zeros((len(whr), 2)), dtype='int')
        mainprogs[:, 0] = tree[whr]
        mainprogs[:, 1] = tree[whr+1]
        sorting = np.argsort(mainprogs[:, 0])
        mainprogs = mainprogs[sorting]
        if len(mainprogs)>0:
            mainprogs = np.char.mod('%15d', mainprogs)
        else:
            mainprogs = np.array([['', '']])
        
        pd.DataFrame(mainprogs).to_csv(outdirc+
                '/CLUSTER_%04d_prog_output.snapshot_%03d.AHF.MainProg.txt' % (
                c, snap_rng[-1-i]), sep='\t', index=None, 
                header=['           halo', '           prog'])

    return None



#find_progenitors(2000, loaddir, savedir)
#find_progenitors(2001, loaddir, savedir)