PRO save_ps_mean,saveplot

dir="/sps/lsst/data/rcecile/Planck_BAO2_PS/"
root = "_lfZuccaAllFalse"

print,'root = ',root
namez = ['0.5','0.9','1.5']
nx= ['_120','_225', '_320']
nx_cell= [120.,225., 320.]
nz_cell= [125.,125.,175.]
cell=8.
namezmean=['0.51', '0.93', '1.57']
lpsym=[0,4,0]
llin = [2,0,0]

nameerr = ['', '_err0.03', '_errP','_errPBDT9','_errPBDT8']
lcol  = [95,210,35,135,  120]
lerr=['spectroZ','Gaussian 0.03','photoZ','photoZ BDT 90%','photoZ BDT 80%']

loadct,12

nz = n_elements(namez)
nerr = n_elements(nameerr)
nsim=10

kmax = intarr(n_elements(namez), n_elements(nameerr))
kmax[*,0] = [20,20,12]
kmax[*,1] = [20,20,12]
kmax[*,2] = [20,20,12]
kmax[*,3] = [18,18,10]
kmax[*,4] = [16,16,10]

restore,'/sps/lsst/users/rcecile/BAOProgs/DirectSimPerso/PlotScript/temp_ascii_new.sav'   
myformat = '(F9.7," ",G14.6," ",G14.6," ",G14.6)'
!p.thick=3
!p.charsize=1.5

if (saveplot ) then begin
; current plotting device.
      mydevice = !D.NAME
      SET_PLOT, 'PS'
      DEVICE, FILENAME='/sps/lsst/users/rcecile/Fig/ps_ready2fit'+root+'.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=6.6,FONT_SIZE=5
   endif else window,10

;plot,[0.01,0.1],[1,1],/xs,/ys,xra=[0.02,0.19],/nodata,yra=[.2,150],xtit='wavenumber k [Mpc^-1]',ytit='Power spectrum x k^2 [(Mpc)^-5]',xma=[4.5,1],yma=[3.5,2.5],ytickn=replicate(" " ,10),/yl,/xl
plot,[0.01,0.1],[1,1],/xs,/ys,xra=[0.02,0.19],/nodata,yra=[30,80000],xtit='wavenumber k [Mpc^-1]',ytit='Power spectrum [(Mpc)^-5]',xma=[4.5,1],yma=[3.5,2.5],ytickn=replicate(" " ,10),/yl,/xl

err=dblarr(nz,nerr,nsim)
minSN=[0.3,0.15,0.15,0.1,0.1]
for iz = 0,nz-1 do begin
   t  = dir + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   xt=ps.(0)
   pt=ps.(1)
   ; oplot,xt,pt*xt*xt,li=2
   oplot,xt,pt,li=2

   for ir = 0, nerr-1 do begin
      suff = nameerr[ir]
      for igd =0,4 do begin
         
                                ; < spectre observe with osc >
         fobs  = dir + 'PS'+root+nx[iz]+'_z'+namez[iz]+suff+'_G'+strtrim(igd,2)+'_wngal.txt' 
         print,fobs
         
         pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_4col)     
         if (igd eq 0) then begin
            xs = pobs.(0)
            p_osc =pobs.(1)
            p_sn = pobs.(2)
         endif else begin
            p_osc = p_osc + pobs.(1)
            p_sn = p_sn + pobs.(2)
         endelse
      endfor

      p_osc = p_osc / 5.      
      p_sn = p_sn / 5. 
  
      SN = mean(p_sn[where(pobs.(0) gt minSN[ir])])
      p2fit = p_osc - SN

  ;    oplot,xs,(p2fit)*xs*xs,col=lcol[ir],li=llin[iz],psym=lpsym[iz]
      oplot,xs,(p2fit),col=lcol[ir],li=llin[iz],psym=lpsym[iz]

      delta_k = ( (pobs.(0))[15]-(pobs.(0))[10] ) / 5. ; pour moyenner eventuelle erreur d'arrondi
      Vsurvey = nx_cell[iz]*nx_cell[iz]*nz_cell[iz]*cell*cell*cell * 5.
      CteSigmaPk = 2*!pi/sqrt(delta_k * Vsurvey)
      sig_osc = CteSigmaPk /pobs.(0)* (p2fit + SN)
      
                                ; errplot,xs,(p2fit-sig_osc)*xs*xs,(p2fit+sig_osc)*xs*xs,col=lcol[ir]
   endfor
endfor

legend,[lerr],li=0,col=[lcol],th=3,box=1,/fill,/right,/top,charsize=1.25
legend,'z='+namez,li=[2,0,0],psym=[-3,4,-3],box=1,col=lcol[4],/fill,/left,/bottom,charsize=1.5,th=3

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif
                                ;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/ps_cata.jpg' ,tvrd(true=3),true=3

END
