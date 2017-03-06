PRO plot_ps_noosc2,saveplot
nsuff=['','_errPBDT8']
mycol=[80,205,240]

!p.charsize=2

loadct,39
restore,'temp_ascii.sav'        ; contient dir
dirn="/sps/lsst/data/rcecile/Planck_noBAO_PS/"

namez = ['0.5','0.9','1.5']
namezmean = ['0.51','0.93','1.57']
nx=['_120','_225','_320']
nz = n_elements(namez)
nk = n_elements(namek)
mytext=['spectroZ','Gauss 0.03','photoZ with podds cut']
mytext = ['theoretical spectrum wo BAO',mytext]
kBAO =2.*!pi/150.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps_noosc_disp.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif 

!p.multi=[0,1,2]
minmy=[0.,4]
maxmy=[0.5,0.]
myxt=['','wavenumber [Mpc^-1]']
txt=['Spectroscopic redshift','Photometric redshift' ]
txt2=['','after 80% BDT cut' ]
for isuff=0,1 do begin
   if (isuff eq 0) then plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.01,0.3],xtit=myxt[isuff],ytit=' PS / <PS>',/nodata,$
        yra=[0.1,1.9],xma=[9,1],yma=[minmy[isuff],maxmy[isuff]],/xl,xtickname=replicate(" ",10) else $

      plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.01,0.3],xtit=myxt[isuff],ytit=' PS / <PS>',/nodata,$
      yra=[0.1,1.9],xma=[9,1],yma=[minmy[isuff],maxmy[isuff]],/xl

   suff = nsuff[isuff]
   for iz=2,0,-1 do begin
      
      for ig=0,4 do begin
         nref=0
         for j=0,9 do begin
            t  = dirn + 'PS_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_G'+strtrim(ig,2)+'_wngal.txt' 
                                ; print,t
            check = FILE_TEST(t)
            if (check eq 0 ) then continue
                                ; print,'                read now'
            p = read_ascii(t, template =  TEMP_POW_SPEC_TXT) 
            if (j eq 0 ) then   PSref=dblarr(n_elements(p.(0)))
            if (j eq 0 ) then   SNref=dblarr(n_elements(p.(0)))
            PSref = PSref+p.(1)
            SNref = SNref+p.(6)
            nref = nref+1
         endfor
         PSref = PSref/nref
         SNref = SNref/nref
         ;print,'REF computed from ',nref,' spectra'
         
         for j=0,9 do begin
            t  = dirn + 'PS_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_G'+strtrim(ig,2)+'_wngal.txt' 
                                ; print,t
            check = FILE_TEST(t)
            if (check eq 0 ) then continue
                                ; print,'                read now'
            p = read_ascii(t, template =  TEMP_POW_SPEC_TXT) 
            xs = p.(0)
       ;     oplot,xs,(p.(1)-p.(6))/(PSref-SNref),col=mycol[iz],th=1
            oplot,xs,(p.(6))/(SNref),col=mycol[iz],th=1
         endfor
      endfor
   endfor
   xyouts,.016,0.325,txt[isuff],charsize=2
   xyouts,.016,0.15,txt2[isuff],charsize=2
   if (isuff eq 0) then legend,'z = '+ namez,col=mycol,line=0,box=1,/fill,/right,/top,charsize=1.5,th=4
   oplot,[0.01,0.4],[1,1],li=2,th=3 ;,col=210
   oplot,[kBAO,kBAO],[0,2],th=4,col=0,li=2

endfor
xyouts,0.046,1.525,'150 Mpc',ori=90 

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif
   
END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
