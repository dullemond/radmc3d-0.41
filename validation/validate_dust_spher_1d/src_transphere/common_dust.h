c     --------------------------------------------------------------
c                THE GLOBAL COMMONS FOR THE DUST TEMPERATURES
c                              AND OPACITIES
c     --------------------------------------------------------------
      common/idustopac/dust_done_read,dust_opacity_tempdep
      integer dust_done_read,dust_opacity_tempdep
      common/dustopac/dust_kappawgt_abs,dust_kappawgt_scat,
     %                dust_freq_wgt,dust_tmin,dust_tmax,
     %                dust_dtmin,dust_dtmax,dust_temprange_high,
     %                dust_temprange_low
      doubleprecision dust_kappawgt_abs(1:FRSIZE_FREQ,1:DUST_TRANGE_MAX,
     %                    1:DUST_SIZE_MAX,1:DUST_SPECIES_MAX)
      doubleprecision dust_kappawgt_scat(1:FRSIZE_FREQ,1:DUST_TRANGE_MAX,
     %                    1:DUST_SIZE_MAX,1:DUST_SPECIES_MAX)
      doubleprecision dust_freq_wgt(1:FRSIZE_FREQ)
      doubleprecision dust_temprange_low(1:DUST_TRANGE_MAX,
     %                    1:DUST_SPECIES_MAX)
      doubleprecision dust_temprange_high(1:DUST_TRANGE_MAX,
     %                    1:DUST_SPECIES_MAX)
      doubleprecision dust_tmin(1:DUST_SPECIES_MAX)
      doubleprecision dust_tmax(1:DUST_SPECIES_MAX)
      doubleprecision dust_dtmin(1:DUST_SPECIES_MAX)
      doubleprecision dust_dtmax(1:DUST_SPECIES_MAX)
c
      common/idust/dust_nr_species,dust_nr_size,
     %             dust_frwgt_read,dust_warn_few_freqs,
     %             dust_nr_temp,dust_warn_zero_temp
      integer dust_nr_species,dust_warn_zero_temp
      integer dust_frwgt_read,dust_warn_few_freqs
      integer dust_nr_size(1:DUST_SPECIES_MAX)
      integer dust_nr_temp(1:DUST_SPECIES_MAX)
c
      common/idustsetup/dust_setup_nrspecies,dust_setup_nrsizes
      integer dust_setup_nrspecies
      integer dust_setup_nrsizes(1:DUST_SPECIES_MAX)
c
