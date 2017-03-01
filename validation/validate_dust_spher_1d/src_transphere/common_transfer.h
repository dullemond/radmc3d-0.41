      common/irays/rays_ncnr,rays_ncex,rays_ncst,
     %             rays_nr,rays_ir,rays_ns,rays_istar,
     %             rays_imu,rays_itang,rays_nrefs,rays_node
      integer rays_ncnr,rays_ncex,rays_nr,rays_ncst,rays_nrefs
      integer rays_istar(FRSIZE_RAYS)
      integer rays_ir(FRSIZE_S,FRSIZE_RAYS)
      integer rays_ns(FRSIZE_RAYS)
      integer rays_itang(FRSIZE_RAYS)
      integer rays_imu(FRSIZE_S,FRSIZE_RAYS)
      integer rays_node(FRSIZE_S,FRSIZE_RAYS)
      common/rays/rays_b,rays_s,rays_dr
      doubleprecision rays_b(FRSIZE_RAYS)
      doubleprecision rays_s(FRSIZE_S,FRSIZE_RAYS)
      doubleprecision rays_dr(FRSIZE_S,FRSIZE_RAYS)
c
      common/oneray/oneray_dr,oneray_s
      doubleprecision oneray_dr(FRSIZE_S)
      doubleprecision oneray_s(FRSIZE_S)
      common/ioneray/oneray_ir,oneray_node
      integer oneray_ir(FRSIZE_S)
      integer oneray_node(FRSIZE_S)
c
      common/inodes/node_nmu,node_iray
      integer node_nmu(FRSIZE_R)
      integer node_iray(-FRSIZE_RAYS:FRSIZE_RAYS,FRSIZE_R)
      common/nodes/node_mu,node_dmu
      doubleprecision node_mu(-FRSIZE_RAYS:FRSIZE_RAYS,FRSIZE_R)
      doubleprecision node_dmu(-FRSIZE_RAYS:FRSIZE_RAYS,FRSIZE_R)
c
      common/iiter/iter_max,iter_method,iter_qdr,iter_fluxcorr
      integer iter_max,iter_method,iter_qdr,iter_fluxcorr
      common/iter/iter_convcrit
      doubleprecision iter_convcrit
c
      common/radfield/rays_intensity,node_intensity
      doubleprecision rays_intensity(
     %               2*FRSIZE_R+3,FRSIZE_RAYS)
      doubleprecision node_intensity(
     %            -FRSIZE_RAYS:FRSIZE_RAYS,FRSIZE_R)
c
c     The planck-mean dust opacity
c
      common/dopac/dopac_temp,dopac_kappa
      doubleprecision dopac_temp(NR_TEMPERATURE)
      doubleprecision dopac_kappa(NR_TEMPERATURE)
      common/idopac/dopac_nt
      integer dopac_nt
c
c     Commons for the spectrum and images
c
      common/specim/spim_b,spim_bi,spim_intens,spim_spectrum
      doubleprecision spim_b(FRSIZE_RAYS),spim_bi(0:FRSIZE_RAYS)
      doubleprecision spim_intens(FRSIZE_RAYS,FRSIZE_FREQ)
      doubleprecision spim_spectrum(FRSIZE_FREQ)
      common/ispecim/spim_nb
      integer spim_nb
c
