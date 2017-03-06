PRO test_plot_ps_osc,icase

nsuff=['','_err0.03','_errP','_errPBDT9','_errPBDT8']
proj = '_cube'
;proj = '_shell'

!p.charsize=2
!p.thick=3
loadct,12
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_noBAO_PS/"

namez = ['0.5','0.9','1.5']
namezmean = ['0.525','0.96','1.54']
nx=['_160','_300','_225']

nline=[2,2,0,0];,0,0]
nz = n_elements(namez)
nk = n_elements(namek)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0

lcol  = [0,95,35,  135, 120]
lcol  = [0,95,210,35,135,  120]


!p.multi=0
plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.0015,0.19],/yl,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum [(Mpc)^-3]',/nodata,yra=[2,5e4],xma=[9,1],yma=[3,1]

for iz=0,nz-1 do begin


   for isuff=0,n_elements(nsuff)-1 do begin

      suff = nsuff[isuff]

      nsim=0
      t  = dir + 'PS'+proj+'_'+strtrim(icase,2)+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
      print,t
      check = FILE_TEST(t)
      if (check eq 0 ) then continue
      nsim = nsim +  1
      p = read_ascii(t, template =  TEMP_POW_SPEC_TXT) 
      xs = p.(0)
      pref = (p.(1)-p.(6))
      print,suff,xs[10],pref[10]
      oplot,xs,pref,col=lcol[isuff+1];,li=nline[iz];,th=lth[isuff];,psym=lpsym[isuff+1],symsize=2
     
   endfor
   
   
   t  = dir + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   print,(ps.(1))[10:12]
   xt=ps.(0)
   pt=ps.(1)
   oplot,xt,pt,li=2
  

endfor   

mytext=['spectroZ','Gauss 0.03','photoZ','photoZ with BDT 90% cut','photoZ with BDT 80% cut']
mytext = ['theoretical spectrum with BAO',mytext]
legend,mytext,line=[2,0,0,0,0,0],col=lcol,box=1,/fill,/left,/bottom,charsize=1.5
;legend,'z = '+ namez,line=indgen(nz)*2,box=1,/fill,/right,/top,charsize=1.5,th=4



end
