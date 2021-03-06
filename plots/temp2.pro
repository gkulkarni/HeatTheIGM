
; File: temp.pro 
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; Plot temperature evolution. 

;; set_plot, 'ps'
;; device, filename='t0.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

set_plot, 'ps'
device, filename='t0.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

;; window, xsize=1000, ysize=1000
;; Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4
TvLCT, 127, 0, 153, 5
TvLCT, 255, 255, 0, 6
TvLCT, 1, 3, 64, 7
!P.charsize = 1.5
!P.thick = 1.0
!P.charthick = 1

restore, 'reionfiletemplate.sav'
; reiondata = read_ascii('set89/reion.out', template=reionfiletemplate)
reiondata = read_ascii('set5/reion.out', template=reionfiletemplate)
redshift = reiondata.z
q = reiondata.q
tempc = reiondata.field06
temph = reiondata.field05
tempav = reiondata.field07
temphvolav = reiondata.field14
tva = temphvolav 

for i = 0, 99 do begin 
   
   tva[i] = q[i]*temphvolav[i] + (1.0-q[i])*tempc[i] 

endfor

plot, redshift, tva, /xlog, /ylog, xrange=[1,100], xtitle='!6z', ytitle='T!D0!N (K)' 
oplot, redshift, tempc, linestyle=5, color=3
oplot, redshift, temphvolav, linestyle=5, color=2

vline, 8.5, linestyle=2
xyouts, 7.5, 1.0e3, 'z!Dreion!N', orientation=90.0, charsize=2.0, alignment=0.5

readcol, '../data/t0_schaye.dat', x, y, dxms, dxps, dyms, dyps 
rs_schaye = x 
t0_schaye = y*1.0e3  
plotsym, 0, 0.5, /FILL
oplot, rs_schaye, t0_schaye, psym=8, color=5
dyms = -dyms*1.0e3 
dyps = dyps*1.0e3
plot_err, rs_schaye, t0_schaye, dyms, dy2=dyps, dx1=-dxms, dx2=dxps, color=5

readcol, '../data/t0_lidz.dat', x, y, dxml, dxpl, dyml, dypl 
rs_lidz = x 
t0_lidz = y*1.0e3  
oplot, rs_lidz, t0_lidz, psym=8, color=2
dyml = -dyml*1.0e3 
dypl = dypl*1.0e3
plot_err, rs_lidz, t0_lidz, dyml, dy2=dypl, dx1=-dxml, dx2=dxpl, color=2

readcol, '../data/t0_becker.dat', zlow, zhigh, z, t0, t0_err 
dzm = z-zlow 
dzp = zhigh-z 
dy=t0_err*0.5 
oplot, z, t0, psym=8, color=3
plot_err, z, t0, dy, dx1=dzm, dx2=dzp, color=3

legend, ['Schaye et al. 00', 'Lidz et al. 10', 'Becker et al. 11',  'HI regions', $
         'HII regions', 'average'], linestyle=[0,0,0,5,5,0], psym=[8,8,8,0,0,0], $
        color=[5,2,3,3,2,-1], /bottom, charsize=1

;; plot, redshift, tva, xrange=[1.5,5.5], xstyle=1
;; oplot, rs_schaye, t0_schaye, psym=8, color=4
;; plot_err, rs_schaye, t0_schaye, dyms, dy2=dyps, dx1=-dxms, dx2=dxps, color=4 
;; oplot, rs_lidz, t0_lidz, psym=8, color=2
;; plot_err, rs_lidz, t0_lidz, dyml, dy2=dypl, dx1=-dxml, dx2=dxpl, color=2
;; oplot, z, t0, psym=8, color=3
;; plot_err, z, t0, dy, dx1=dzm, dx2=dzp, color=3

device, /close_file
set_plot, 'X'

END
