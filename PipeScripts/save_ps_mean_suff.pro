PRO save_ps_mean_suff,refsuff,psuff

; save_ps_mean_suff,"_noBias_noCata",""
; save_ps_mean_suff,"_noBias",""
; save_ps_mean_suff,"_noCata",""
; save_ps_mean_suff,"",""

; save_ps_mean_suff,"_SFspec",""

namez = ['0.5','0.9','1.5']
namezmean=['0.51', '0.93', '1.57']
nx= ['_120','_225', '_320']
nameerr = ['_errP','_errPBDT9','_errPBDT8']
;nameerr = ['', '_err0.03', '_errPBDT9','_errPBDT8']

loadct,12

lcol  = [35,135,  120,0]
;lcol  = [95,210,135,  120,0]
lerr=['photoZ','photoZ BDT 90%','photoZ BDT 80%','fiducial, no oscillation']
;lerr=['spectroZ','Gaussian 0.03','photoZ BDT 90%','photoZ BDT 80%','theoretical, no oscillation']
lpsym=[0,-5,-4]
llin = [0,3,2]

nz = n_elements(namez)
nerr = n_elements(nameerr)
nsim=10

restore,'temp_ascii.sav'   
dir="/sps/lsst/data/rcecile/Planck_BAO_PS/"
dirn="/sps/lsst/data/rcecile/Planck_noBAO_PS/"
myformat = '(F9.7," ",G14.6," ",G14.6," ",G14.6," ",G14.6," ",G14.6," ",G14.6," ",G14.6," ",I6," ",I6," ",G14.6)'
!p.thick=3
!p.charsize=2

plot,[0.01,0.1],[1,1],/xs,/ys,xra=[0.005,0.15],/nodata,yra=[1100.,.5e5],/yl,xtit='wavenumber [Mpc^-1]',ytit='Undamped Power spectrum',xma=[4,0.5],yma=[3.,.5],ytickn=replicat

err=dblarr(nz,nerr,nsim)
for iz = 0,nz-1 do begin
   
   ; spectre th??orique no osc, no error
   ftheo  = dirn + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
   xt=ptheo.(0)
   
   for ir = 0, nerr-1 do begin
      suff = nameerr[ir]
      for igd =0,4 do begin
                                ; spectre observ??
         fobs  = dir + 'PS'+psuff+nx[iz]+'_z'+namez[iz]+suff+'_G'+strtrim(igd,2)+'_wngal.txt' 
                                ;    print,fobs
         pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_TXT)   
         if (igd eq 0) then p2fit =pobs.(1) * 0.
         xs = pobs.(0)
         ok = intarr(n_elements(xs))
         for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))
         
      ; < spectre simul?? no osc, with error >
         pnobao=dblarr(n_elements(xs))
         pshot=dblarr(n_elements(xs),nsim)
         mynsim = 0
         for j=0,nsim-1 do begin
            fnobao  = dirn + 'PS'+refsuff+'_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_G'+strtrim(igd,2)+'_wngal.txt' 
            print,fnobao
            check = FILE_TEST(fnobao)
            if (check eq 0) then continue
            mynsim = mynsim+1
            p = read_ascii(fnobao, template =  TEMP_POW_SPEC_TXT)     
            pnobao = pnobao + p.(1) - p.(6)
            pshot[*,j] = p.(6)
         endfor
         pnobao = pnobao / double(mynsim)
         print,suff, ' N sim OK =',mynsim
         std_shot = dblarr(n_elements(xs))
         for i=0,n_elements(xs)-1 do std_shot[i] = sqrt(mean(pshot[i,where(pshot[i,*] ne 0)]))
         okerr = where(p.(0) gt 0.02 and p.(0) lt 0.1)
         for j=0,nsim-1 do begin
            fnobao  = dirn + 'PS'+refsuff+'_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_G'+strtrim(igd,2)+'_wngal.txt' 
            check = FILE_TEST(fnobao)
            if (check eq 0) then continue
            p = read_ascii(fnobao, template =  TEMP_POW_SPEC_TXT)    
            rap = p.(1)/pnobao
            err[iz,ir,j] = stddev(rap[okerr])
                                ;   print,j
         endfor
                                ;    print,iz,ir,' ',mean(err[iz,ir,*])
         
         
        p2fit =  p2fit + (pobs.(1)-pobs.(6)) / (pnobao) *(ptheo.(2)[ok])
      endfor
      p2fit = p2fit / 5.

      oplot,xs,p2fit,col=lcol[ir]
      
      fname = dir + 'PS'+refsuff+psuff+'_2fit'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
                                ;    print,fname
      openw,lun1, fname, /get_lun
      
      newp1 = p2fit + pobs.(6)
      p7shot = sqrt(pobs.(7)*pobs.(7) + std_shot*std_shot)
      print,'SHOT  ',iz, ir, (pobs.(7))[25],p7shot[25]
      newp7 = (newp1-pobs.(6))/(pobs.(1)-pobs.(6))*p7shot ; erreur "undampee", comme le signal
      newp7 = newp7*sqrt(1.+1./10.)                       ; car 10 simus sans BAOs pour la base
;      newp7 = newp7 * mean(err[iz,ir,where(err[iz,ir,*] ne 0)])/mean(err[iz,0,where(err[iz,ir,*] ne 0)]); augmentation de la dispersion en cas d'erreur non spectro : facteur introduit par err sur z
                                ;    print,'ICI ',  mean(err[iz,ir,*]), ' et ', mean(err[iz,0,okerr])
                                ;    print,'  ' ,suff,(pobs.(0))[40],(pobs.(7))[40],newp7[40],(newp1[40]-(pobs.(6))[40])/((pobs.(1)[40]-pobs.(6))[40])
      
;read,xx
      for i=0,n_elements(p2fit)-1 do  printf, lun1, (pobs.(0))[i], $
                                              newp1[i], $ ;spectre
                                              (pobs.(2))[i], (pobs.(3))[i], (pobs.(4))[i], (pobs.(5))[i], (pobs.(6))[i], $
                                              newp7[i], $ ; erreur
                                              (pobs.(8))[i], (pobs.(9))[i], (pobs.(10))[i],format=myformat
      close, lun1 
      free_lun, lun1
                                ;    print,'  ',fname
   endfor
   oplot,xt,ptheo.(2),li=2
endfor
;legend,'z = '+ namez+nsuff,col=lcol,box=1,line=0,/fill,/left,/bottom,charsize=1.5
;legend,'z = '+ namez,col=lcol,box=1,line=0,/fill,/left,/bottom,charsize=1.5
legend,lerr,li=[replicate(0,n_elements(nameerr)),2],col=lcol,box=1,/fill,/right,/top,charsize=1.25

 

END
