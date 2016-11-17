PRO self,doplot
; doplot = 1 : plot on screen, = 2 : save .eps file

restore,'temp_ascii.sav'
!p.charsize=2.
!p.multi=0
!p.thick=3
loadct,12
lcol  = [105,35,210,  135, 120]
lcol2=[105, 210,135,120]

;dirg='/sps/lsst/data/rcecile/TJP_BAO_grids/'
;dirng='/sps/lsst/data/rcecile/TJP_noBAO_grids/'
dirg='/sps/lsst/data/rcecile/Planck_BAO_grids/'
dirng='/sps/lsst/data/rcecile/Planck_noBAO_grids/'

; may want to add some selection function

;st = [dirg+'SelFunc_G2_full_specz_nofz.txt']
;sg = [dirg+'SelFunc_G2_Gauss0.03_full_nofz.txt']
;sp = [dirg+'SelFunc_G2_errP_full_nofz.txt']
;sd = [dirg+'SelFunc_G2_errPpodds_full_nofz.txt']
;sb = [dirg+'SelFunc_G2_errPBDT_full_nofz.txt']

st = [dirg+'SelFunc__specz_nofz.txt']
sg = [dirg+'SelFunc__Gauss0.03_full_nofz.txt']
sp = [dirg+'SelFunc__errP_full_nofz.txt']
sd = [dirg+'SelFunc__errPpodds_full_nofz.txt']
sb = [dirg+'SelFunc__errPBDT_full_nofz.txt']

nCase=10
;for i=0,nCase-1 do st = [st, dirng+'SelFunc_G2_'+strtrim(i,2)+'_full_specz_nofz.txt']
;for i=0,nCase-1 do sg = [sg, dirng+'SelFunc_G2_'+strtrim(i,2)+'_Gauss0.03_full_nofz.txt']
;for i=0,nCase-1 do sp = [sp, dirng+'SelFunc_G2_'+strtrim(i,2)+'_errP_full_nofz.txt']
;for i=0,nCase-1 do sd = [sd, dirng+'SelFunc_G2_'+strtrim(i,2)+'_errPpodds_full_nofz.txt']
;for i=0,nCase-1 do sb = [sb, dirng+'SelFunc_G2_'+strtrim(i,2)+'_errPBDT_full_nofz.txt']

for i=0,nCase-1 do st = [st, dirng+'SelFunc__'+strtrim(i,2)+'_specz_nofz.txt']
for i=0,nCase-1 do sg = [sg, dirng+'SelFunc__'+strtrim(i,2)+'_Gauss0.03_full_nofz.txt']
for i=0,nCase-1 do sp = [sp, dirng+'SelFunc__'+strtrim(i,2)+'_errP_full_nofz.txt']
for i=0,nCase-1 do sd = [sd, dirng+'SelFunc__'+strtrim(i,2)+'_errPpodds_full_nofz.txt']
for i=0,nCase-1 do sb = [sb, dirng+'SelFunc__'+strtrim(i,2)+'_errPBDT_full_nofz.txt']


if (doplot eq 2) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
;   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/self.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=8.8,FONT_SIZE=4
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/self_G1.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=8.8,FONT_SIZE=4
endif


print,'DOPLOT ',doplot
!p.multi=[0,2,3]
ns0 = 0
for i=0,nCase do begin
   print,st[i],ns0
   check = FILE_TEST(st[i])
   if (check eq 0) then continue
   print,'           read OK '
   ns0 = ns0+1
   s = read_ascii(st[i], template =  TEMP_SEL_FUNC)
   z = s.(0)
   print,s.(1)
   if (ns0 eq 1) then s1 = s.(1)
   if (ns0 eq 1) then specz = s.(1) else specz = specz + s.(1)
   if (doplot eq 1) then $
      if (ns0 eq 1) then plot,z,s.(1)/s1,/xs,/ys,ytit='Selection fct ratio',xra=[0.1,2.65],yra=[0.9,1.1],yma=[0,4],xma=[8,0] $
      else oplot,z,s.(1)/s1,col=80*i
endfor
if (doplot eq 1) then $
   xyouts,0.5,1.05,'spectroZ'
if (ns0 ge 1) then specz = specz/ ns0

ns1 = 0
for i=0,nCase do begin
   print,sg[i],ns1
   check = FILE_TEST(sg[i])
   if (check eq 0) then continue
   print,'           read OK '
   ns1 = ns1+1
   s = read_ascii(sg[i], template =  TEMP_SEL_FUNC)
   if (ns1 eq 1) then s1 = s.(1)
   if (ns1 eq 1) then gauss = s.(1) else gauss = gauss + s.(1)
   if (doplot eq 1) then $
      if (ns1 eq 1) then plot,z,s.(1)/s1,/xs,/ys,xra=[0.1,2.65],yra=[0.9,1.1],yma=[0,4],xma=[0,4] $
      else oplot,z,s.(1)/s1,col=80*i
endfor
if (doplot eq 1) then $
   xyouts,0.5,1.05,'Gaussian error: 0.03*(1+z)'
if (ns1 ge 1) then gauss = gauss / ns1

ns2 = 0
for i=0,nCase do begin
   print,sp[i],ns2
   check = FILE_TEST(sp[i])
   if (check eq 0) then continue
   print,'           read OK '
   ns2 = ns2+1
   s = read_ascii(sp[i], template =  TEMP_SEL_FUNC)
   if (ns2 eq 1) then s1 = s.(1)
   if (ns2 eq 1) then photz  = s.(1) else photz = photz + s.(1)
   if (doplot eq 1) then $
      if (ns2 eq 1) then plot,z,s.(1)/s1,/xs,/ys,ytit='Selection fct ratio',xtit='redshift',xra=[0.1,2.65],yra=[0.9,1.1],yma=[4,0],xma=[8,0] $
      else oplot,z,s.(1)/s1,col=80*i
endfor
if (doplot eq 1) then $
   xyouts,0.5,1.05,'Photometric error'
if (ns2 ge 1) then photz = photz / ns2

ns3 = 0
for i=0,nCase do begin
   print,sd[i],ns3
   check = FILE_TEST(sd[i])
   if (check eq 0) then continue
   print,'           read OK '
   ns3 = ns3+1
   s = read_ascii(sd[i], template =  TEMP_SEL_FUNC)
   if (ns3 eq 1) then s1 = s.(1)
   if (ns3 eq 1) then podds  = s.(1) else podds = podds + s.(1)
   if (doplot eq 1) then $
      if (ns3 eq 1) then plot,z,s.(1)/s1,/xs,/ys,xtit='redshift',xra=[0.1,2.65],yra=[0.9,1.1],yma=[4,0],xma=[0,4] $
      else oplot,z,s.(1)/s1,col=80*i
endfor
if (doplot eq 1) then $
   xyouts,0.5,1.05,'Photometric error + odds cut'
if (ns3 ge 1) then podds = podds/ ns3

ns4 = 0
for i=0,nCase do begin
   print,sb[i],ns4
   check = FILE_TEST(sb[i])
   if (check eq 0) then continue
   print,'           read OK '
   ns4 = ns4+1
   s = read_ascii(sb[i], template =  TEMP_SEL_FUNC)
   if (ns4 eq 1) then s1 = s.(1)
   if (ns4 eq 1) then bdt  = s.(1) else bdt  = bdt + s.(1)
   if (doplot eq 1) then $
      if (ns4 eq 1) then plot,z,s.(1)/s1,/xs,/ys,ytit='Selection fct ratio',xtit='redshift',xra=[0.1,2.65],yra=[0.9,1.1],yma=[4,0],xma=[8,0] $
      else oplot,z,s.(1)/s1,col=80*i
endfor
if (doplot eq 1) then $
   xyouts,0.5,1.05,'Photometric error + BDT cut'
if (ns4 ge 1) then bdt = bdt/ ns4
if (doplot eq 1) then read,xx


!p.multi=[0,1,2]
plot,[0.1,2.65],[0.00002,1],/xs,/ys,ytit='Selection function',/yl,xra=[0.05,2.7],yra=[0.00002,1],yma=[0,4],/nodata,xmar=[8,1]
if (ns0 ge 1) then oplot,z,specz,col=lcol[0],li=2
;if (ns1 ge 1) then oplot,z,gauss,col=lcol[1];,li=2
if (ns2 ge 1) then oplot,z,photz,col=lcol[2]
if (ns3 ge 1) then oplot,z,podds,col=lcol[3]
if (ns4 ge 1) then oplot,z,bdt,col=lcol[4]
oplot,[2.2,2.2],[1e-5,1],th=1,col=250,li=2
oplot,[0.2,0.2],[1e-5,1],th=1,col=250,li=2

what=['spectroZ','photoZ','photoZ + odds cut','photoZ + BDT cut']
legend,what,col=lcol2,line=[2,0,0,0],box=1,psym=0,/fill,/left,/bottom,charsize=1.25

plot,[0.1,2.65],[0,2],/xs,/ys,ytit='Sel. function ratio',xtit='Redshift',xra=[0.05,2.7],yma=[4,0],yra=[0,2],/nodata,xmar=[8,1]
if (ns0 ge 1) then oplot,z,specz/specz,col=lcol[0],li=2
;if (ns1 ge 1) then oplot,z,gauss/specz,col=lcol[1];,li=2
if (ns2 ge 1) then oplot,z,photz/specz,col=lcol[2]
if (ns3 ge 1) then oplot,z,podds/specz,col=lcol[3]
if (ns4 ge 1) then oplot,z,bdt/specz,col=lcol[4]

oplot,[2.2,2.2],[0,3],th=1,col=250,li=2
oplot,[0.2,0.2],[0,3],th=1,col=250,li=2

if (doplot eq 2) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif


;stop
   
openw,lun0, '/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_speczmean_nofz.txt', /get_lun
if (ns0 ge 1) then for i=0,n_elements(z)-1 do  printf, lun0, z[i], specz[i]
close, lun0
free_lun, lun0
   
openw,lun1, '/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_gaussmean_nofz.txt', /get_lun
if (ns1 ge 1) then for i=0,n_elements(z)-1 do  printf, lun1, z[i], gauss[i]
close, lun1
free_lun, lun1
   
openw,lun2, '/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_photzmean_nofz.txt', /get_lun
if (ns2 ge 1) then for i=0,n_elements(z)-1 do  printf, lun2, z[i], photz[i]
close, lun2
free_lun, lun2
   
openw,lun3, '/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_poddsmean_nofz.txt', /get_lun
if (ns3 ge 1) then for i=0,n_elements(z)-1 do  printf, lun3, z[i], podds[i]
close, lun3
free_lun, lun3

   
openw,lun4, '/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_pbdtmean_nofz.txt', /get_lun
if (ns4 ge 1) then for i=0,n_elements(z)-1 do  printf, lun4, z[i], bdt[i]
close, lun4
free_lun, lun4

if (doplot eq 2) then stop

print,'1 if you want to check the written selection functions, q if not'
read,xx

restore,'temp_ascii.sav'
!p.multi=0
s0 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_speczmean_nofz.txt', template =  TEMP_SEL_FUNC)
s1 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_G1_specz_nofz.txt', template =  TEMP_SEL_FUNC)
s1 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_gaussmean_nofz.txt', template =  TEMP_SEL_FUNC)
s2 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_photzmean_nofz.txt', template =  TEMP_SEL_FUNC)
s3 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_poddsmean_nofz.txt', template =  TEMP_SEL_FUNC)
s4 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_G2_errPBDT_full_nofz.txt', template =  TEMP_SEL_FUNC)
plot,[0.1,2.65],[0.0002,1],/xs,/ys,ytit='Mean selection function',/yl,xra=[0.1,3.],yra=[0.00002,1],/nodata
oplot,s0.(0),s0.(1),col=lcol[0]
oplot,s1.(0),s1.(1),col=lcol[1]
oplot,s2.(0),s2.(1),col=lcol[2]
oplot,s3.(0),s3.(1),col=lcol[3]

 oplot,s4.(0),s4.(1),col=0

;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/sf_bdt.jpg' ,tvrd(true=3),true=3
s0 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_G1_specz_nofz.txt', template =  TEMP_SEL_FUNC)
s1 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc70_14__specz_nofz.txt', template =  TEMP_SEL_FUNC)
s2 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc140_28__specz_nofz.txt', template =  TEMP_SEL_FUNC)
plot,[0.1,2.65],[0.0002,1],/xs,/ys,ytit='Mean selection function',/yl,xra=[0.1,0.9],yra=[0.1,1],/nodata
oplot,s0.(0),s0.(1),col=lcol[0]
oplot,s1.(0),s1.(1),col=lcol[1],li=2
oplot,s2.(0),s2.(1),col=lcol[2]

!p.multi=[0,1,2]
s0 = read_ascii('/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_G2_full_ZONLY.txt', template =  TEMP_SEL_FUNC)
s1 = read_ascii('/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_G2_full_specz_nofz.txt', template =  TEMP_SEL_FUNC)
plot,s0.(0),s0.(1)*s1.(1),/xs,/ys,xra=[0.1,.25],psym=-1



!p.multi=0
s0 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc__ZONLY.txt', template =  TEMP_SEL_FUNC)
s1 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc__specz_nofz.txt', template =  TEMP_SEL_FUNC)
plot,s0.(0),s0.(1)*s1.(1),/xs,/ys,xra=[0.16,.4],psym=-1


s2 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc140_28__ZONLY.txt', template =  TEMP_SEL_FUNC)
s3 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc140_28__specz_nofz.txt', template =  TEMP_SEL_FUNC)
oplot,s0.(0),s2.(1)*s3.(1),col=123,psym=-1
oplot,s0.(0),s2.(1)*s1.(1),col=234,li=2

plot,s0.(0),s0.(1)*s1.(1)-s0.(1)*s3.(1),/xs,/ys,xra=[0.16,.36],psym=-1


!p.multi=0
s1 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc__ZONLY.txt', template =  TEMP_SEL_FUNC)
plot,s1.(0),s1.(1),/xs,/ys,xra=[0.16,.4],psym=-1

s3 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc140_28__ZONLY.txt', template =  TEMP_SEL_FUNC)
oplot,s3.(0),s3.(1),col=123,psym=-1



!p.multi=0
s1 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc__specz_nofz.txt', template =  TEMP_SEL_FUNC)
s5 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc25__specz_nofz.txt', template =  TEMP_SEL_FUNC)
s2 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_G25TA_n15__specz_nofz.txt', template =  TEMP_SEL_FUNC)
s3 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_G40TA_n15__specz_nofz.txt', template =  TEMP_SEL_FUNC)
what=['avant','hier (pas gold)','type25 gold','type40 gold']

plot,s1.(0),s1.(1),/xs,/ys,xra=[0.1,.38],yra=[0.3,1.02]
oplot,s5.(0),s5.(1),col=23
oplot,s2.(0),s2.(1),col=100
oplot,s3.(0),s3.(1),col=200,li=2

legend,what,col=[0,23,100,200],line=[0,0,0,2],box=1,psym=0,/fill,/left,/bottom,charsize=2


end
