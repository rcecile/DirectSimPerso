PRO plot_ps_4suff_z,iw,iz,icase
; plot_ps_4suff_z,10,0,'' (BAO,1er z)
; plot_ps_4suff_z,21,1,'5' (noBAO, 2ieme z)

!p.charsize=2
h= 0.679
h3 = h*h*h

; si nz =200
; kz =   0. 0.00392699   0.00785398    0.0117810
; si nz =300
; kz =   0.  0.00261799   0.00523599   0.00785398

loadct,39

  restore,'temp_ascii.sav' ; contient dir
  if (strlen(icase) eq 0) then begin
     dir="/sps/lsst/data/rcecile/TJP_BAO_PS/"
     ssuff = ''
  endif
  if (strlen(icase) gt 0) then begin
     dir="/sps/lsst/data/rcecile/TJP_noBAO_PS/"
     ssuff = '_'+icase
  endif

  nz = ['0.7','1.4']
  nx= ['_450', '_875']

  suff=['_err0.03','_errP','_errPpodds']
  n = n_elements(suff)

  namez = nz[iz]

  t  = dir + 'simu_ps'+nx[iz]+'_z'+namez+'_ntpk.txt' 
  print,'theory ',t
  p = read_ascii(t, template =  TEMP_POW_SPEC_TH)
  pt = dblarr(2,n_elements(p.(0)))
  pt[0,*]=p.(0)
  if (strlen(icase) gt 0) then  pt[1,*]=p.(2) else  pt[1,*]=p.(1)

  t  = dir + 'PS_G2'+ssuff+nx[iz]+'_z'+namez+'_wngal.txt'
  print,t
  check = FILE_TEST(t)
  if (check eq 1 ) then begin
     p = read_ascii(t, template =  TEMP_POW_SPEC_TXT)
     ps = dblarr(2,n_elements(p.(0)))
     ps[0,*]=p.(0)
     ps[1,*]=p.(1)-p.(6)
  endif

  for i=0,n-1 do begin
     t  = dir + 'PS_G2'+ssuff+nx[iz]+'_z'+namez+suff[i]+'_wngal.txt'
     print,t
     check = FILE_TEST(t)
     if (check eq 0 ) then stop
     print,'    read'
     p = read_ascii(t, template =  TEMP_POW_SPEC_TXT)
     if (i eq 0) then  ps3 = dblarr(2*n,n_elements(p.(0)))
     ps3[2*i,*]=p.(0)
     ps3[2*i+1,*]=p.(1)-p.(6)
  endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  !p.multi=0
  
  ltext ='z = '+ namez 
  lcol  = 60*findgen(5)
  lline= 0
  lpsym=[0,4,0,0,0]
  lerr=['theoretical shape','spectroZ',' Gaussian 0.03','photoZ','photoZ with podds cut']

!p.thick=2
 
window,iw

lcol  = [0,80,150,195,210];60*findgen(5)

plot,ps[0,*],ps[1,*],/xs,/ys,xra=[0.0025,0.2],/yl,yra=[20,10e4],xtit='wavenumber [Mpc^-1]',ytit='Power spectrum [(Mpc)^-3]',/nodata
oplot,ps[0,*] , ps[1,*],col=lcol[1],psym=lpsym[1],symsize=2
for i=0,n-1 do oplot,ps3[2*i,*],ps3[2*i+1,*],col=lcol[i+2],th=4
oplot,pt[0,*], pt[1,*],col=lcol[0]

legend,lerr,col=lcol,psym=lpsym,box=1,/fill,/left,/bottom,charsize=1.5,th=4
legend,ltext,box=1,/fill,/right,/top,charsize=1.5


end
