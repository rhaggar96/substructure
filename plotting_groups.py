from modules import *

clis = np.loadtxt('/home/ppxrh2/Documents/test_pollux/TheThreeHundred/'
            'playground/rhaggar/G3X_data/G3X_300_selected_sample_257.txt')
clis = np.array(clis, dtype='int')

def plot_paths(crange):
    """ Finds the paths of galaxies in an infalling group  """

    

    """
    masses, r200, mstar = np.zeros(0), np.zeros(0), np.zeros(0)
    for c in crange:
        print c
        data_in = ld_arr('/run/media/ppxrh2/166AA4B87A2DD3B7/NewMDCLUSTER_data/'
                'NewMDCLUSTER/NewMDCLUSTER_%04d/GadgetX-NewMDCLUSTER_%04d'
                '.snap_128.z0.000.AHF_halos' % (c, c))
        nlim = np.sum(data_in[:, 1] > 99)
        data_in = data_in[:nlim]
        masses = np.append(masses, data_in[:, 0])
        r200 = np.append(r200, data_in[:, 8])
        mstar = np.append(mstar, data_in[:, 61])
    np.savetxt('data_temp_masses.txt', masses)
    np.savetxt('data_temp_r200.txt', r200)
    np.savetxt('data_temp_mstar.txt', mstar)
    """
    masses = ld_arr('data_temp_masses.txt')
    r200 = ld_arr('data_temp_r200.txt')
    mstar = ld_arr('data_temp_mstar.txt')


    r200 = r200[mstar > 10.**9.5]
    masses = masses[mstar > 10.**9.5]
    mstar = mstar[mstar > 10.**9.5]

    #mstar = mstar[r200<3.]
    #masses = masses[r200<3.]
    #r200 = r200[r200<3.]

    stell_f = mstar / masses

    plt.figure()
    plt.scatter(r200, stell_f, s=1.)
    plt.loglog()
    plt.ylim(0.00001, 1.)
    plt.show()


    flim = 0.3
    a = np.sum(stell_f>flim)
    b = len(stell_f)
    print float(a)/float(b)

    mstar = np.log10(mstar+1.)
    masses = np.log10(masses)

    masses = masses[mstar>1.]
    mstar = mstar[mstar>1.]

    plt.figure()
    plotting = plt.hist2d(masses, mstar, bins=200, 
            range=[[9., 16.], [7., 14.]], norm=SymLogNorm(linthresh=2.))[0]
    plt.imshow(np.transpose(plotting), norm=SymLogNorm(linthresh=2.), 
            extent = [9., 16., 7., 14.], origin='lower', 
            interpolation = 'gaussian')
    plt.colorbar()
    plt.show()
    
    masses = 10.**masses
    mstar = 10.**mstar

    plt.figure(figsize=(7, 5))
    plt.scatter(masses, mstar, s=0.2)
    plt.plot([10.**9., 10.**16.], [10.**9., 10.**16.], c='k')
    plt.plot([10.**9., 10.**16.], [flim*10.**9., flim*10.**16.], c='k')
    plt.xlim(10.**9., 10.**16.)
    plt.ylim(10.**7., 10.**14.)
    plt.loglog()
    #plt.savefig('m_comp.pdf')
    #plt.savefig('m_comp.png', dpi=500)
    plt.show()
    """
    plt.figure()
    plt.hist(stell_f)
    plt.semilogy()
    #plt.show()

    smallest = np.argsort(np.array(r200, dtype='float'))[:5]
    """
    

    return None

plot_paths(clis)
"""
            
            ids_t = ids_i[np.where(bs[0]==True)[0], inf_snap[i]]
            print ids_i[inf_id[i]][-1]
            xs_t = xs_i[np.where(bs[0]==True)[0], inf_snap[i]:] - xs_i[c_h,inf_snap[i]:]
            ys_t = ys_i[np.where(bs[0]==True)[0], inf_snap[i]:] - ys_i[c_h,inf_snap[i]:]
            zs_t = zs_i[np.where(bs[0]==True)[0], inf_snap[i]:] - zs_i[c_h,inf_snap[i]:]
            
            halo_id = ids_i[inf_id[i], -1]
            x_host = xs_i[inf_id[i], inf_snap[i]:] - xs_i[c_h,inf_snap[i]:]
            y_host = ys_i[inf_id[i], inf_snap[i]:] - ys_i[c_h,inf_snap[i]:]
            z_host = zs_i[inf_id[i], inf_snap[i]:] - zs_i[c_h,inf_snap[i]:]
            
            



            if inf_snap[i] < 100 and len(xs_t) > 4:
                plt.close()
                fig = plt.figure(figsize=(16, 8))
                ax = fig.add_subplot(121, projection='3d')
                ax.grid(False)
                viewang = [15, 181]
                #viewang = [20, 135]
                ax.view_init(viewang[0], viewang[1])
                for nh in range(len(xs_t)):
                    nh_bool = xs_t[nh]>-100000.
                    ax.plot(xs_t[nh][nh_bool], ys_t[nh][nh_bool], zs_t[nh][nh_bool], linewidth=1.)
                    if nh==0:
                        ax.scatter(xs_t[nh, 0], ys_t[nh, 0], zs_t[nh, 0], s=10., alpha=0.)
                    else:
                        ax.scatter(xs_t[nh, 0], ys_t[nh, 0], zs_t[nh, 0], s=10.)
                u = np.linspace(0, 2*pi, 100)
                v = np.linspace(0, pi, 100)
                x = np.outer(np.cos(u), np.sin(v))
                y = np.outer(np.sin(u), np.sin(v))
                z = np.outer(np.ones(np.size(u)), np.cos(v))
                ax.plot_surface(rvirs_i[c_h,inf_snap[i]] * x, rvirs_i[c_h,inf_snap[i]] * y, rvirs_i[c_h,inf_snap[i]] * z, color='k', alpha=0.15)
                ax.plot_surface(x_host[0]+rvirs_i[inf_id[i],inf_snap[i]] * x, 
                        y_host[0]+rvirs_i[inf_id[i],inf_snap[i]] * y, 
                        z_host[0]+rvirs_i[inf_id[i],inf_snap[i]] * z, color='b', alpha=0.2)
                ax.set_xlim(-2950, 2950)
                ax.set_ylim(-2950, 2950)
                ax.set_zlim(-2950, 2950)
                ax.set_xlabel(r'$x$ / kpc/$h$')
                ax.set_ylabel(r'$y$ / kpc/$h$')
                ax.set_zlabel(r'$z$ / kpc/$h$')
                ax.xaxis.labelpad=40
                ax.yaxis.labelpad=10
                ax.zaxis.labelpad=20
                ax.tick_params(axis='x', pad=15)
                ax.tick_params(axis='z', pad=10)
                #plt.show()
                #plt.close()
                #fig = plt.figure(figsize=(16, 8))
                ax = fig.add_subplot(122, projection='3d')
                ax.grid(False)
                ax.view_init(viewang[0], viewang[1])
                for nh in range(len(xs_t)):
                    nh_bool = xs_t[nh]>-100000.
                    ax.plot(xs_t[nh][nh_bool], ys_t[nh][nh_bool], zs_t[nh][nh_bool], linewidth=1.)
                    if nh==0:
                        ax.scatter(xs_t[nh, -1], ys_t[nh, -1], zs_t[nh, -1], s=10., alpha=0.)
                    else:
                        ax.scatter(xs_t[nh, -1], ys_t[nh, -1], zs_t[nh, -1], s=10.)
                u = np.linspace(0, 2*pi, 100)
                v = np.linspace(0, pi, 100)
                x = np.outer(np.cos(u), np.sin(v))
                y = np.outer(np.sin(u), np.sin(v))
                z = np.outer(np.ones(np.size(u)), np.cos(v))
                ax.plot_surface(rvirs_i[c_h,-1] * x, rvirs_i[c_h,-1] * y, rvirs_i[c_h,-1] * z, color='k', alpha=0.15)
                ax.plot_surface(x_host[-1]+rvirs_i[inf_id[i],-1] * x, 
                        y_host[-1]+rvirs_i[inf_id[i],-1] * y, 
                        z_host[-1]+rvirs_i[inf_id[i],-1] * z, color='b', alpha=0.2)
                ax.set_xlim(-2950, 2950)
                ax.set_ylim(-2950, 2950)
                ax.set_zlim(-2950, 2950)
                ax.set_xlabel(r'$x$ / kpc/$h$')
                ax.set_ylabel(r'$y$ / kpc/$h$')
                ax.set_zlabel(r'$z$ / kpc/$h$')
                ax.xaxis.labelpad=40
                ax.yaxis.labelpad=10
                ax.zaxis.labelpad=20
                ax.tick_params(axis='x', pad=15)
                ax.tick_params(axis='z', pad=10)
                plt.tight_layout()
                plt.savefig(('data_out/1group_path_clust_%03d_halo_'%c)+str(halo_id)+'.pdf')
                plt.savefig(('data_out/1group_path_clust_%03d_halo_'%c)+str(halo_id)+'.png', dpi=300)
                plt.show()


"""