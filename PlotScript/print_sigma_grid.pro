PRO print_sigma_grid

loadct,12
lcol  = [95,35,120]

dir='/sps/lsst/data/rcecile/Planck_BAO_grids/'

restore,'temp_ascii.sav'        ; contient dir

nameerr=['Photo-z','Photo-z BDT 90%','Photo-z BDT 80%']
err=['','_BDT90','_BDT80']
gerr=['Photo-z','Photo-z BDT 90%','Photo-z BDT 80%']
namez = ['0.5','0.9','1.5']
nx= ['_120','_225', '_320']

xh = dindgen((0.06-0.01)/0.0005+1)*0.0005+0.01
for id=0,n_elements(namez)-1 do begin
   grid=dir+'grids'+ strtrim(nx[id],2)+'_z'+namez[id]+'_5cubes.fits'
   c=mrdfits(grid,6,h)
   window,id
   plot,minmax(xh),[.1,100],/xs,/ys,th=3,yra=[0,30],/nodata
   for i=0,n_elements(err) -1 do begin
      
      f='/sps/lsst/data/rcecile/Planck_BAO_grids/Sigma_photoZ'+err[i]+'.txt'
      s = read_ascii(f, template =  TEMP_SEL_FUNC)
      
      sig_grid = interpol(s.(1),s.(0),c) ; fill with IQR
      sig_grid = sig_grid / 1.34         ; translate in sigma

      h = histogram(sig_grid,min=0.01,max=0.06,bin=0.0005)
      oplot,xh,h/total(h)*100.,th=3,col=lcol[i]
      print,namez[id], err[i], median(c),median(sig_grid),mean(c),mean(sig_grid),stddev(sig_grid)
   endfor

endfor

END

