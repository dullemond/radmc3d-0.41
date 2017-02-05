@readradmc
@natconst

xsize=17.78
ysize=16.70
thick=2.0      ; Postscript line thickness
size=1.7       ; Postscript char size

a=read_data(/dtemp)

!p.thick=thick
!x.thick=thick
!y.thick=thick
!z.thick=thick
!p.charthick=thick
!p.charsize=size
set_plot,'ps'
device,file='temperature.ps',xsize=xsize,ysize=ysize
surface,a.temp,(a.grid.r-a.grid.r[0])/1d5,(!dpi/2-a.grid.theta)*180./!dpi,$
        xtitle='Height [km]',ytitle='Latt (toward star)',ax=50,/ys,$
        ztitle='T [K]'
device,/close
set_plot,'x'
!p.thick=1
!x.thick=1
!y.thick=1
!z.thick=1
!p.charthick=1
!p.charsize=1


end

