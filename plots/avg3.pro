
; File: avg3.pro
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; This is a slave code for avg_multi.pro.  See that file to understand
; what this does. 

openr, lun, 'cs2_51_znbyfe.dat', /get_lun
cs_dat = fltarr(2,5)
readf, lun, cs_dat
rs = cs_dat[0,*]
cbyfe_avg = cs_dat[1,*]
plotsym, 0, 0.7, /FILL
; plot, rs, cbyfe_avg, psym=-8, yrange=[-1.0,1.0], xtitle='z',
; ytitle='[Zn/Fe]'
plot, rs, cbyfe_avg, psym=-8, yrange=[-1.0,1.0], ytitle='[Zn/Fe]', xrange=[2,6.5], xstyle=1
close, lun

openr, lun, 'cs2_49_znbyfe.dat', /get_lun
cs_dat = fltarr(2,5)
readf, lun, cs_dat
rs = cs_dat[0,*]
cbyfe_avg = cs_dat[1,*]
plotsym, 0, 0.7, /FILL
oplot, rs, cbyfe_avg, psym=-8, color=2
close, lun

openr, lun, 'cs2_50_znbyfe.dat', /get_lun
cs_dat = fltarr(2,5)
readf, lun, cs_dat
rs = cs_dat[0,*]
cbyfe_avg = cs_dat[1,*]
plotsym, 0, 0.7, /FILL
oplot, rs, cbyfe_avg, psym=-8, color=3
close, lun

;; device, /close_file
;; set_plot, 'X'
