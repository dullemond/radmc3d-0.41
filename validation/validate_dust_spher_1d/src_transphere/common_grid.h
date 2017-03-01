c
c     Commons for the grid
c
      common/grid/grid_r,             ! Radial grid
     %            freq_nu             ! Frequency grid
      doubleprecision grid_r(FRSIZE_R)
      doubleprecision freq_nu(FRSIZE_FREQ)
      common/igrid/grid_nr,           ! Number of radial grid points
     %             freq_nr            ! Number of frequency points
      integer grid_nr,freq_nr
c
