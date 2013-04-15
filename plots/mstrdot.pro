
window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4 
TvLCT, 1, 3, 64, 5 
!P.charsize = 2
!P.multi=1

restore, 'halo_template.sav'
mstrdot_data = read_ascii('set17/halos_aux.out', template=stars_template)
mstrdot = mstrdot_data.field001[col,*]*1.0e10 ; Msun / yr (Check!)
z = mstrdot_data.field001[0,*] 
plot, z, mstrdot, /xlog, /ylog, xrange=[0.1,100] 


