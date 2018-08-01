
PRO plot_lf,ifile,itypeplot
  common sch_var, Mstar, phi_star, alpha

file=[ 'lfparamRamosAll.txt', 'lfparamZuccaAll.txt', 'lfparamDahlenAll_ini.txt','lfparamZuccaAllFixed.txt','lfparamDahlenAll.txt']

use_all =0

dir='/sps/lsst/users/rcecile/BAOProgs/DirectSimPerso/PipeScripts/'

;#   LF_Params = [ phistar mstar alpha ]1
;#  zmin zmax LF_Params_Ell LF_Params_Spiral LF_Params_StarBurst LF_Params_All 

restore,'temp_ascii_new.sav'
!p.charsize=1.5
!p.multi=0
!p.thick=3
lf = read_ascii(dir+file[ifile], template =  TEMP_LF)
dMstar = 5.*alog(0.70)-5.*alog(0.679)
cphi_star = (0.679/0.7)^3
PRINT,' MODIF h70 -> hPlanck : dM, cphi = ',dMstar,cphi_star
niz = n_elements(lf.(0))

dz = 0.05
zmin=0.1
zmax=2.5
z = findgen( (zmax-zmin)/dz+1)*dz+zmin
nz = n_elements(z)

;Mmax = dindgen(1)-13.
;Mmax = dindgen(1)-19.
Mmax = dindgen(1)-20.5
Mmin = dblarr(n_elements(Mmax))-24.

lfz = dblarr(nz)
lfallz = dblarr(nz)
lfzE = lfz
lfzS = lfz
lfzL = lfz
window,10+ifile*2+10
for im=0,n_elements(Mmin)-1 do begin
   for iz=0,nz-1 do begin
      
      ibinz = 0
      while (z[iz] ge lf.(1)[ibinz] and ibinz lt niz-1) do  ibinz = ibinz+1
      
      phi_star =1.d * lf.(11)[ibinz] * cphi_star
      Mstar =1.d * lf.(12)[ibinz] + dMstar
      alpha = 1.d * lf.(13)[ibinz]
      si_all = QROMB('fsch_par', Mmin[im], Mmax[im], /double,jmax=1000,K=5)
      lfallz[iz] = si_all
      
      phi_star =1.d * lf.(2)[ibinz] * cphi_star
      Mstar =1.d * lf.(3)[ibinz] + dMstar
      alpha = 1.d * lf.(4)[ibinz]
      si_E = QROMB('fsch_par', Mmin[im], Mmax[im], /double,jmax=1000,K=5)
      
      phi_star =1.d * lf.(5)[ibinz] * cphi_star
      Mstar =1.d * lf.(6)[ibinz] + dMstar
      alpha = 1.d * lf.(7)[ibinz]
      si_L = QROMB('fsch_par', Mmin[im], Mmax[im], /double,jmax=1000,K=5)
      
      phi_star =1.d * lf.(8)[ibinz]  * cphi_star
      Mstar =1.d * lf.(9)[ibinz] + dMstar
      alpha = 1.d * lf.(10)[ibinz]
      si_S = QROMB('fsch_par', Mmin[im], Mmax[im], /double,jmax=1000,K=5)
      
      lfz[iz] = (si_E + si_L + si_S)
      lfzE[iz] = si_E
      lfzL[iz] = si_L
      lfzS[iz] = si_S
   endfor

   mycol=[225,125,50]

   if (itypeplot eq 0) then begin ; Fraction
      if (im eq 0) then plot,z,lfz,/xs,/ys,yra=[0.01,1.2],xra=[0.1,1.125],/nodata,xtit='redshift',ytit='LF integrale',tit=file[ifile]+' Mag min='+strtrim(Mmin,2)+', Mag max= '+strtrim(Mmax,2),/yl
      oplot,z,lfzE/lfz,col=mycol[0]
      oplot,z,lfzL/lfz,col=mycol[1]
      oplot,z,lfzS/lfz,col=mycol[2]
      print,max(lfz), Mmin[im], Mmax[im]
   endif

   if (itypeplot eq 1) then begin ; LF
      if (im eq 0) then plot,z,lfallz,/xs,/ys,xra=[0.1,2.5],/yl,yra=[min(lfallz)/1000.,max(lfallz)*10.],xtit='redshift',ytit='LF integrale fractions',tit=file[ifile]
      oplot,z,lfzE,col=mycol[0]
      oplot,z,lfzL,col=mycol[1]
      oplot,z,lfzS,col=mycol[2]
      print,max(lfz), Mmin[im], Mmax[im]
   endif
print,lfzE/lfz
print,lfzL/lfz
endfor
legend,['Early','Late','SB'],col=mycol,li=0,/fill,/bottom,/right

;write_jpeg,'/sps/lsst/users/rcecile/BAO_InOut/frac_'+file[ifile]+strtrim(Mmax,2)+'.jpg' ,tvrd(true=3),true=3

END
