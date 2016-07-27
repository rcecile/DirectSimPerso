
PRO save_ps_mean

namez = ['0.5','0.65','1.0','1.3','1.5']
nx= ['_300', '_432', '_675', '_864', '_960']

nameerr = ['', '_err0.03', '_errPpodds']
lcol  = [80, 100, 150, 210, 240]
lerr=['spectroZ','Gaussian 0.03','photoZ podds']
lpsym=[0,-5,-4]
llin = [0,3,2]

nz = n_elements(namez)
nerr = n_elements(nameerr)
nsim=10

restore,'temp_ascii.sav'   
dir="/sps/lsst/data/rcecile/PaperBAO_cell8_PS/"
dirn="/sps/lsst/data/rcecile/PaperNoBAO_cell8_PS/"
myformat = '(F9.7," ",G14.6," ",G14.6," ",G14.6," ",G14.6," ",G14.6," ",G14.6," ",G14.6," ",I6," ",I6," ",G14.6)'
loadct,39
!p.thick=3
!p.charsize=2
plot,[0.01,0.1],[1,1],/xs,/ys,xra=[0.005,0.15],/nodata,yra=[1000,.5e5],/yl,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum'

err=dblarr(nz,nerr,nsim)
for iz = 0,nz-1 do begin
   
   ; spectre théorique no osc, no error
   ftheo  = dir + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
   ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
   xt=ptheo.(0)
   
   for ir = 0, nerr-1 do begin
      suff = nameerr[ir]

      ; spectre observé
      fobs  = dir + 'PS_G2sfmean'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
  ;    print,fobs
      pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_TXT)     
      xs = pobs.(0)
      ok = intarr(n_elements(xs))
      for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))

      ; < spectre simulé no osc, with error >
      pnobao=dblarr(n_elements(xs))
      for j=0,nsim-1 do begin
         fnobao  = dirn + 'PS_G2sfmean_nobao'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
         p = read_ascii(fnobao, template =  TEMP_POW_SPEC_TXT)     
         pnobao = pnobao + p.(1)
      endfor
      pnobao = pnobao / double(nsim)
      okerr = where(p.(0) gt 0.05 and p.(0) lt 0.15)
      for j=0,nsim-1 do begin
         fnobao  = dirn + 'PS_G2sfmean_nobao'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
         p = read_ascii(fnobao, template =  TEMP_POW_SPEC_TXT)    
         rap = p.(1)/pnobao
         err[iz,ir,j] = stddev(rap[okerr])
      endfor
  ;    print,iz,ir,' ',mean(err[iz,ir,*])


;      p2fit = (pobs.(1)-pobs.(6)) / (pnobao-pobs.(6)) *(ptheo.(2)[ok]-pobs.(6))
      p2fit = (pobs.(1)-pobs.(6)) / (pnobao-pobs.(6)) *(ptheo.(2)[ok])

      oplot,xs,p2fit,col=lcol[iz],thick=5-3*ir;psym=lpsym[ir]
      oplot,xt,ptheo.(2),li=2

      fname = dir + 'PS_G2sfmean2fitErr'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
      openw,lun1, fname, /get_lun
 
      newp1 = p2fit + pobs.(6)

      newp7 = (newp1-pobs.(6))/(pobs.(1)-pobs.(6))*pobs.(7) ; erreur "undampee", comme le signal
      newp7 = newp7*sqrt(1.+1./10.)                     ; car 10 simus sans BAOs pour la base
      newp7 = newp7 * mean(err[iz,ir,*])/mean(err[iz,0,*]); augmentation de la dispersion en cas d'erreur non spectro : facteur introduit par err sur z
      print,'  ' ,suff,(pobs.(0))[40],(pobs.(7))[40],newp7[40],(newp1[40]-(pobs.(6))[40])/((pobs.(1)[40]-pobs.(6))[40])
      for i=0,n_elements(p2fit)-1 do  printf, lun1, (pobs.(0))[i], newp1[i],(pobs.(2))[i], (pobs.(3))[i], (pobs.(4))[i], $
                                                     (pobs.(5))[i], (pobs.(6))[i], newp7[i], (pobs.(8))[i], (pobs.(9))[i], (pobs.(10))[i],format=myformat
      close, lun1
      free_lun, lun1
  ;    print,'  ',fname
   endfor
endfor
legend,'z = '+ namez,col=lcol,box=1,line=0,/fill,/left,/bottom,charsize=1.5
legend,lerr,li=0,box=1,/fill,/center,/bottom,charsize=1.5,thick=[5,3,1]

for iz = 0,nz-1 do begin
   for ir = 0, nerr-1 do begin
      print,namez[iz],nameerr[ir],mean(err[iz,ir,*])/mean(err[iz,0,*])*100.,' %'
   endfor
endfor

;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/plot_ps2fit_'+strtrim(n_elements(namez),2)+'.jpg' ,tvrd(true=3),true=3

END
