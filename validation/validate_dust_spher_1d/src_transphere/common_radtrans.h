      common/spherrt/srt_rho,srt_temp,srt_qfpi
      doubleprecision srt_rho(FRSIZE_R)
      doubleprecision srt_temp(FRSIZE_R)
      doubleprecision srt_qfpi(FRSIZE_R)
c
      common/varedd/vet_jj,vet_hh,vet_jnu,vet_hnu,vet_fj,
     %      vet_fnu,vet_kapt,vet_kapj,vet_kaph,
     %      vet_jjme,vet_hhme,vet_corfact
      doubleprecision vet_hhme(FRSIZE_R)
      doubleprecision vet_jjme(FRSIZE_R)
      doubleprecision vet_hh(FRSIZE_R)
      doubleprecision vet_jj(FRSIZE_R)
      doubleprecision vet_fj(FRSIZE_R)
      doubleprecision vet_kapt(FRSIZE_R)
      doubleprecision vet_kapj(FRSIZE_R)
      doubleprecision vet_kaph(FRSIZE_R)
      doubleprecision vet_jnu(FRSIZE_FREQ,FRSIZE_R)
      doubleprecision vet_hnu(FRSIZE_FREQ,FRSIZE_R)
      doubleprecision vet_fnu(FRSIZE_FREQ,FRSIZE_R)
      doubleprecision vet_corfact(FRSIZE_R)
c
      common/iglobpar/global_algorithm,global_idump,
     %            global_init,global_dump_intens,
     %            global_ncst,global_ncex,global_ncnr,
     %            global_itypemw
      integer global_algorithm,global_idump,
     %    global_init,global_dump_intens,
     %    global_ncst,global_ncex,global_ncnr,
     %    global_itypemw
c
      common/opactable/optab_kappa
      doubleprecision optab_kappa(FRSIZE_FREQ)
c
      common/vet/vet_convcrit
      doubleprecision vet_convcrit
      common/ivet/vet_itermax
      integer vet_itermax
c
c     Commons for the central star
c
      common/star/star_t,             ! Temperature of central star
     %            star_r,             ! Radius of central star
     %            star_m,             ! Mass of star
     %            star_lum            ! Spectrum of the star in lum
      doubleprecision star_r,star_t,star_m
      doubleprecision star_lum(FRSIZE_FREQ)
c
