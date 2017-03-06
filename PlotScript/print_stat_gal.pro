PRO print_stat_gal

nx=['_120','_225','_320']
namez = ['0.5','0.9','1.5']

for iz=0,2 do begin
   ntot = 0.
   nmean= 0.
   ngal = dblarr(5)
   for i=0,4 do begin
      m5=mrdfits("/sps/lsst/data/rcecile/Planck_BAO_grids/grids_"+nx+"_z"+namez[iz]+"_5cubes.fits",i,h)
      ngal  = sxpar(h,'NGRID')

      print,'iz = ', namez[iz], 'grid n',i
      print,'  minmax z = ',minmax(m5)
      print,ngal
   endfor
endfor
END
