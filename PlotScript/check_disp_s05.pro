PRO check_disp_s05


namez = '0.5'
nx= '_120'

loadct,12
nsim=10
restore,'temp_ascii.sav'   
dir="/sps/lsst/data/rcecile/Planck_BAO_PS/"
dirn="/sps/lsst/data/rcecile/Planck_noBAO_PS/"

!p.thick=3


plot,[0.01,0.1],[1,1],/xs,/ys,xra=[0.015,0.19],/nodata,yra=[140.,1050],xtit='wavenumber k [Mpc^-1]',ytit='(Undamped PS) x k',xma=[4.5,1],yma=[3.5,.5],ytickn=replicate(" " ,10),/yl,/xl
mycol=0
for igd =0,4 do begin
                                 ; < spectre simule no osc >
   for j=0,nsim-1 do begin
      fnobao  = dirn + 'PS_SN_'+strtrim(j,2)+nx+'_z'+namez+'_G'+strtrim(igd,2)+'_wngal.txt' 
      p = read_ascii(fnobao, template =  TEMP_POW_SPEC_TXT)  
      if (mycol eq 0 ) then plot,p.(0),p.(1)*p.(0),/xs,/ys,xra=[0.015,0.19],/xl,/yl,yra=[200,1.5e3]    else oplot,p.(0),p.(1)*p.(0),col=mycol
      mycol = mycol+5

   endfor
endfor
 write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/disp_z05.jpg' ,tvrd(true=3),true=3


END
