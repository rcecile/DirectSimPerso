PRO self,doplot
; doplot = 1 : plot on screen, = 2 : save .eps file

restore,'temp_ascii.sav'
!p.charsize=2.
!p.multi=0
!p.thick=3
loadct,12
lcol  = [95,210,35,135, 120]

dirg='/sps/lsst/data/rcecile/Planck_BAO_grids/'
dirng='/sps/lsst/data/rcecile/Planck_noBAO_grids/'

st = [dirg+'SelFunc_gold__specz_nofz.txt']
sg = [dirg+'SelFunc_gold_Gauss0.03_nofz.txt']
sp = [dirg+'SelFunc_gold_errP_nofz.txt']
sb8 = [dirg+'SelFunc_gold_errPBDT8_nofz.txt']
sb9 = [dirg+'SelFunc_gold_errPBDT9_nofz.txt']

nCase=10

for i=0,nCase-1 do st = [st, dirng+'SelFunc_gold__'+strtrim(i,2)+'_specz_nofz.txt']
for i=0,nCase-1 do sg = [sg, dirng+'SelFunc_gold_'+strtrim(i,2)+'_Gauss0.03_nofz.txt']
for i=0,nCase-1 do sp = [sp, dirng+'SelFunc_gold_'+strtrim(i,2)+'_errP_nofz.txt']
for i=0,nCase-1 do sb8 = [sb8, dirng+'SelFunc_gold_'+strtrim(i,2)+'_errPBDT8_nofz.txt']
for i=0,nCase-1 do sb9 = [sb9, dirng+'SelFunc_gold_'+strtrim(i,2)+'_errPBDT9_nofz.txt']


if (doplot eq 2) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
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
   if (ns0 eq 1) then s1 = s.(1)
   if (ns0 eq 1) then specz = s.(1) else specz = specz + s.(1)
   if (doplot eq 1) then $
      if (ns0 eq 1) then plot,z,s.(1)/s1,/xs,/ys,ytit='Selection fct ratio',xra=[0.1,2.65],yra=[0.9,1.1],yma=[0,4],xma=[10,0] $
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
      if (ns2 eq 1) then plot,z,s.(1)/s1,/xs,/ys,ytit='Selection fct ratio',xtit='redshift z_s',xra=[0.1,2.65],yra=[0.9,1.1],yma=[4,0],xma=[8,0] $
      else oplot,z,s.(1)/s1,col=80*i
endfor
if (doplot eq 1) then $
   xyouts,0.5,1.05,'Photometric error'
if (ns2 ge 1) then photz = photz / ns2

ns3 = 0
for i=0,nCase do begin
   print,sb9[i],ns3
   check = FILE_TEST(sb9[i])
   if (check eq 0) then continue
   print,'           read OK '
   ns3 = ns3+1
   s = read_ascii(sb9[i], template =  TEMP_SEL_FUNC)
   if (ns3 eq 1) then s1 = s.(1)
   if (ns3 eq 1) then bdt9  = s.(1) else bdt9 = bdt9 + s.(1)
   if (doplot eq 1) then $
      if (ns3 eq 1) then plot,z,s.(1)/s1,/xs,/ys,xtit='redshift z_s',xra=[0.1,2.65],yra=[0.9,1.1],yma=[4,0],xma=[0,4] $
      else oplot,z,s.(1)/s1,col=80*i
endfor
if (doplot eq 1) then $
   xyouts,0.5,1.05,'Photometric error + BDT90 cut'
if (ns3 ge 1) then bdt9 = bdt9/ ns3

ns4 = 0
for i=0,nCase do begin
   print,sb8[i],ns4
   check = FILE_TEST(sb8[i])
   if (check eq 0) then continue
   print,'           read OK '
   ns4 = ns4+1
   s = read_ascii(sb8[i], template =  TEMP_SEL_FUNC)
   if (ns4 eq 1) then s1 = s.(1)
   if (ns4 eq 1) then bdt8  = s.(1) else bdt8  = bdt8 + s.(1)
   if (doplot eq 1) then $
      if (ns4 eq 1) then plot,z,s.(1)/s1,/xs,/ys,ytit='Selection fct ratio',xtit='redshift z_s',xra=[0.1,2.65],yra=[0.9,1.1],yma=[4,0],xma=[8,0] $
      else oplot,z,s.(1)/s1,col=80*i
endfor
if (doplot eq 1) then $
   xyouts,0.5,1.05,'Photometric error + BDT80 cut'
if (ns4 ge 1) then bdt8 = bdt8/ ns4
if (doplot eq 1) then read,xx



;!p.multi=[0,1,2]
;plot,[0.1,2.65],[0.00002,1],/xs,ys=8,ytit='Selection function',/yl,xra=[0.2,2.45],yra=[0.0001,1],yma=[0,4],/nodata,xmar=[10,10]
;if (ns0 ge 1) then oplot,z,specz,col=lcol[0]
;;if (ns1 ge 1) then oplot,z,gauss,col=lcol[1];,li=2
;if (ns2 ge 1) then oplot,z,photz,col=lcol[2]
;if (ns3 ge 1) then oplot,z,bdt9,col=lcol[3]
;if (ns4 ge 1) then oplot,z,bdt8,col=lcol[4]

;what=['spectroZ','photoZ','photoZ + 90% cut','photoZ + 80% cut']
;legend,what,col=lcol2,line=[2,0,0,0],box=1,psym=0,/fill,/left,/bottom,charsize=1.5


!p.multi=[0,1,2]
plot,[0.1,2.65],[0.00002,1],/xs,ys=8,ytit='log10(SF)',xra=[0.2,2.45],yra=[-4,0],yma=[0,.5],/nodata,xmar=[6,5.5]
if (ns0 ge 1) then oplot,z,alog10(specz),col=lcol[0]
;if (ns1 ge 1) then oplot,z,alog10(gauss),col=lcol[1],li=2
if (ns2 ge 1) then oplot,z,alog10(photz),col=lcol[2]
if (ns3 ge 1) then oplot,z,alog10(bdt9),col=lcol[3]
if (ns4 ge 1) then oplot,z,alog10(bdt8),col=lcol[4]

what=['spectroZ (or Gaussian)','photoZ','photoZ + 90% cut','photoZ + 80% cut']
legend,what,col=[lcol[0],lcol[2:4]],line=0,box=1,psym=0,/fill,/left,/bottom,charsize=1.5

; voir coef dans cat_mpc3.pro
AXIS, YAXIS=1, YRANGE = (!Y.CRANGE)+ alog10(0.10397683), ys=1,  YTITLE = 'log10(d) [Mpc^-3]'


plot,[0.1,2.65],[0,2],/xs,/ys,ytit='SF ratio',xtit='Spectroscopic redshift z_s',xra=[0.2,2.45],yma=[3,0],yra=[0.1,1.75],/nodata,xmar=[6,5.5]
if (ns0 ge 1) then oplot,z,specz/specz,col=lcol[0]
;if (ns1 ge 1) then oplot,z,gauss/specz,col=lcol[1];,li=2
if (ns2 ge 1) then oplot,z,photz/specz,col=lcol[2]
if (ns3 ge 1) then oplot,z,bdt9/specz,col=lcol[3]
if (ns4 ge 1) then oplot,z,bdt8/specz,col=lcol[4]

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
   
openw,lun3, '/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_pbdt9mean_nofz.txt', /get_lun
if (ns3 ge 1) then for i=0,n_elements(z)-1 do  printf, lun3, z[i], bdt9[i]
close, lun3
free_lun, lun3

   
openw,lun4, '/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_pbdt8mean_nofz.txt', /get_lun
if (ns4 ge 1) then for i=0,n_elements(z)-1 do  printf, lun4, z[i], bdt8[i]
close, lun4
free_lun, lun4

if (doplot eq 2) then stop

print,'1 if you want to check the written selection functions, q if not'
read,xx

restore,'temp_ascii.sav'
!p.multi=0
s0 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_speczmean_nofz.txt', template =  TEMP_SEL_FUNC)
s1 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_gaussmean_nofz.txt', template =  TEMP_SEL_FUNC)
s2 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_photzmean_nofz.txt', template =  TEMP_SEL_FUNC)
s3 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_pbdt9mean_nofz.txt', template =  TEMP_SEL_FUNC)
s4 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_pbdt8mean_nofz.txt', template =  TEMP_SEL_FUNC)
plot,[0.1,2.65],[0.0002,1],/xs,/ys,ytit='Mean selection function',/yl,xra=[0.1,3.],yra=[0.00002,1],/nodata
oplot,s0.(0),s0.(1),col=lcol[0]
oplot,s1.(0),s1.(1),col=lcol[1]
oplot,s2.(0),s2.(1),col=lcol[2]
oplot,s3.(0),s3.(1),col=lcol[3]
oplot,s4.(0),s4.(1),col=lcol[4]


plot,s0.(0),s4.(1)/s0.(1),/xs,/ys,ytit='Mean selection function',xra=[0.1,3.]
stop

end
