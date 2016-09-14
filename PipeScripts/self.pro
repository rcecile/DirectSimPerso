PRO self,doplot
; doplot = 1 : plot on screen, = 2 : save .eps file

restore,'temp_ascii.sav'
!p.charsize=2.5
!p.multi=0
!p.thick=2
lcol  = [80, 150, 198, 240]

dirg='/sps/lsst/data/rcecile/TJP_BAO_grids/'
dirng='/sps/lsst/data/rcecile/TJP_noBAO_grids/'
loadct,39

; may want to add some selection function
st = [dirg+'SelFunc_G2_full_specz_nofz.txt']
sg = [dirg+'SelFunc_G2_Gauss0.03_full_nofz.txt']
sp = [dirg+'SelFunc_G2_errP_full_nofz.txt']
sd = [dirg+'SelFunc_G2_errPpodds_full_nofz.txt']

nCase=10
for i=0,nCase-1 do st = [st, dirng+'SelFunc_G2_'+strtrim(i,2)+'_full_specz_nofz.txt']
for i=0,nCase-1 do sg = [sg, dirng+'SelFunc_G2_'+strtrim(i,2)+'_Gauss0.03_full_nofz.txt']
for i=0,nCase-1 do sp = [sp, dirng+'SelFunc_G2_'+strtrim(i,2)+'_errP_full_nofz.txt']
for i=0,nCase-1 do sd = [sd, dirng+'SelFunc_G2_'+strtrim(i,2)+'_errPpodds_full_nofz.txt']


if (doplot eq 2) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/self.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=8.8,FONT_SIZE=4
endif


print,'DOPLOT ',doplot
!p.multi=[0,2,2]
ns0 = 0
for i=0,nCase do begin
   print,st[i],ns0
   check = FILE_TEST(st[i])
   if (check eq 0) then continue
   print,'           read OK '
   ns0 = ns0+1
   s = read_ascii(st[i], template =  TEMP_SEL_FUNC)
   z = s.(0)
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
if (doplot eq 1) then read,xx


!p.multi=[0,1,2]
plot,[0.1,2.65],[0.0002,1],/xs,/ys,ytit='Selection function',/yl,xra=[0.05,2.7],yra=[0.0002,1],yma=[0,4],/nodata,xmar=[8,1]
if (ns0 ge 1) then oplot,z,specz,col=lcol[0]
if (ns1 ge 1) then oplot,z,gauss,col=lcol[1];,li=2
if (ns2 ge 1) then oplot,z,photz,col=lcol[2]
if (ns3 ge 1) then oplot,z,podds,col=lcol[3]
oplot,[2.2,2.2],[1e-5,1],th=1,col=250,li=2
oplot,[0.2,0.2],[1e-5,1],th=1,col=250,li=2

what=['spectroZ','Gaussian 0.03','photoZ','photoZ + odds cut']
legend,what,col=lcol,line=0,box=1,psym=0,/fill,/left,/bottom,charsize=1.8

plot,[0.1,2.65],[0,2],/xs,/ys,ytit='Sel. function ratio',xtit='Redshift',xra=[0.05,2.7],yma=[4,0],yra=[0,2],/nodata,xmar=[8,1]
if (ns0 ge 1) then oplot,z,specz/specz,col=lcol[0]
if (ns1 ge 1) then oplot,z,gauss/specz,col=lcol[1];,li=2
if (ns2 ge 1) then oplot,z,photz/specz,col=lcol[2]
if (ns3 ge 1) then oplot,z,podds/specz,col=lcol[3]

oplot,[2.2,2.2],[0,3],th=1,col=250,li=2
oplot,[0.2,0.2],[0,3],th=1,col=250,li=2

if (doplot eq 2) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif


;stop
   
openw,lun0, '/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_speczmean_nofz.txt', /get_lun
if (ns0 ge 1) then for i=0,n_elements(z)-1 do  printf, lun0, z[i], specz[i]
close, lun0
free_lun, lun0
   
openw,lun1, '/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_gaussmean_nofz.txt', /get_lun
if (ns1 ge 1) then for i=0,n_elements(z)-1 do  printf, lun1, z[i], gauss[i]
close, lun1
free_lun, lun1
   
openw,lun2, '/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_photzmean_nofz.txt', /get_lun
if (ns2 ge 1) then for i=0,n_elements(z)-1 do  printf, lun2, z[i], photz[i]
close, lun2
free_lun, lun2
   
openw,lun3, '/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_poddsmean_nofz.txt', /get_lun
if (ns3 ge 1) then for i=0,n_elements(z)-1 do  printf, lun3, z[i], podds[i]
close, lun3
free_lun, lun3

if (doplot eq 2) then stop

print,'1 if you want to check the written selection functions, q if not'
read,xx

restore,'temp_ascii.sav'
!p.multi=0
plot,[0.1,2.65],[0.0002,1],/xs,/ys,ytit='Mean selection function',/yl,xra=[0.1,2.65],yra=[0.0002,1],/nodata
if (ns0 ge 1) then s0 = read_ascii('/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_speczmean_nofz.txt', template =  TEMP_SEL_FUNC)
if (ns1 ge 1) then s1 = read_ascii('/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_gaussmean_nofz.txt', template =  TEMP_SEL_FUNC)
if (ns2 ge 1) then s2 = read_ascii('/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_photzmean_nofz.txt', template =  TEMP_SEL_FUNC)
if (ns3 ge 1) then s3 = read_ascii('/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_poddsmean_nofz.txt', template =  TEMP_SEL_FUNC)
if (ns0 ge 1) then oplot,s0.(0),s0.(1),col=lcol[0]
if (ns1 ge 1) then oplot,s1.(0),s1.(1),col=lcol[1]
if (ns2 ge 1) then oplot,s2.(0),s2.(1),col=lcol[2]
if (ns3 ge 1) then oplot,s3.(0),s3.(1),col=lcol[3]


end
