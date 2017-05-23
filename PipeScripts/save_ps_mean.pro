PRO save_ps_mean2,saveplot,writefile 

dir="/sps/lsst/data/rcecile/Planck_BAO2_PS/"

namez = ['0.5','0.9','1.5']
;namezmean=['0.525', '0.96', '1.54']
;nx= ['_160','_300', '_225']
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

restore,'temp_ascii.sav'   
myformat = '(F9.7," ",G14.6," ",G14.6," ",G14.6," ",G14.6)'
!p.thick=3
!p.charsize=1.5

if (saveplot ) then begin
; current plotting device.
      mydevice = !D.NAME
      SET_PLOT, 'PS'
      DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/ps_ready2fit_Nz.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=6.6,FONT_SIZE=5
   endif


plot,[0.01,0.1],[1,1],/xs,/ys,xra=[0.01,0.19],/nodata,yra=[.2,94],xtit='wavenumber k [Mpc^-1]',ytit='Power spectrum x k^2 [(Mpc)^-5]',xma=[4.5,1],yma=[3.5,.5],ytickn=replicate(" " ,10),/xl,/yl
plot,[0.01,0.1],[1,1],/xs,/ys,xra=[0.02,0.19],/nodata,yra=[1,94],xtit='wavenumber k [Mpc^-1]',ytit='Power spectrum x k^2 [(Mpc)^-5]',xma=[4.5,1],yma=[3.5,.5],ytickn=replicate(" " ,10),/xl,/yl

err=dblarr(nz,nerr,nsim)
for iz = 0,nz-1 do begin
   
   for ir = 0, nerr-1 do begin
      suff = nameerr[ir]
      for igd =0,4 do begin

        ; < spectre observe with osc >
         fobs  = dir + 'PS_nZ'+nx[iz]+'_z'+namez[iz]+suff+'_G'+strtrim(igd,2)+'_wngal.txt' 
         pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_TXT)     

         if (igd eq 0) then begin
            xs = pobs.(0)
            SN = mean((pobs.(3))[where(pobs.(0) gt 0.2 and pobs.(0) lt 0.4)])
            p_osc =pobs.(1)
         endif else p_osc = p_osc + pobs.(1)
      endfor
      p_osc = p_osc / 5.      
      p2fit = p_osc - SN

      oplot,xs,(p2fit)*xs*xs,col=lcol[ir],li=llin[iz],psym=lpsym[iz]
 
      delta_k = (pobs.(0))[15]-(pobs.(0))[10] / 5. ; pour moyenner eventuelle erreur d'arrondi
      Vsurvey = nx_cell[iz]*nx_cell[iz]*nz_cell[iz]*cell*cell*cell * 5.
      CteSigmaPk = 2*!pi/sqrt(delta_k * Vsurvey)
      sig_osc = CteSigmaPk /pobs.(0)* (p2fit + SN)
 
     ; errplot,xs,(p2fit-sig_osc)*xs*xs,(p2fit+sig_osc)*xs*xs,col=lcol[ir]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      if (writefile) then begin
         fname = dir + 'PS_nZ'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
         print,'OUTPUT in ', fname
         openw,lun1, fname, /get_lun
       
         for i=0,n_elements(p2fit)-1 do  printf, lun1, (pobs.(0))[i], $
                                                 p2fit[i] + SN, $ ;spectre
                                                 (pobs.(2))[i], SN, sig_osc[i],format=myformat
         close, lun1 
         free_lun, lun1
         print,'  ',fname
         endif
    endfor
endfor
legend,[lerr],li=0,col=[lcol],th=3,box=1,/fill,/left,/bottom,charsize=1.25
legend,'z='+namez,li=[2,0,0],psym=[-3,4,-3],box=1,col=lcol[4],/fill,/left,/top,charsize=1.5,th=3

   if (saveplot) then begin
      DEVICE, /CLOSE
      SET_PLOT, mydevice
   endif

END
 ;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/ps_k2.jpg' ,tvrd(true=3),true=3
