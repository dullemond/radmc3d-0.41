@readradmc
@natconst

xsize=17.78
ysize=12.70
thick=3.0      ; Postscript line thickness
size=1.4       ; Postscript char size

a=read_data(/dtemp)
cd,current=current
cd,'../run_spher1d_2_trsph/'
openr,1,'envstruct.dat'
nr=0
readf,1,nr
data=dblarr(3,nr)
readf,1,data
close,1
cd,current

!p.thick=thick
!x.thick=thick
!y.thick=thick
!z.thick=thick
!p.charthick=thick
!p.charsize=size
set_plot,'ps'
device,file='compare_temp.ps'
plot,a.grid.r/au,a.temp,/xl
oplot,transpose(data[0,*])/au,transpose(data[2,*]),psym=6
device,/close
set_plot,'x'
!p.thick=1
!x.thick=1
!y.thick=1
!z.thick=1
!p.charthick=1
!p.charsize=1


end

