PRO self

restore,'temp_ascii.sav'
!p.charsize=2
!p.multi=0
!p.thick=3

dirg='/sps/lsst/data/rcecile/TJP_BAO_grids/'
dirng='/sps/lsst/data/rcecile/TJP_BAO_grids/'
loadct,39

; may want to add some selection function
st = [dirg+'SelFunc_specz_nofz.txt',    dirng+'SelFunc_nobao6_specz_nofz.txt',    dirng+'SelFunc_nobao7_specz_nofz.txt']
sg = [dirg+'SelFunc_Gauss0.03_nofz.txt',dirng+'SelFunc_nobao6_Gauss0.03_nofz.txt',dirng+'SelFunc_nobao7_Gauss0.03_nofz.txt']
sp = [dirg+'SelFunc_errP_nofz.txt',     dirng+'SelFunc_nobao6_errP_nofz.txt',     dirng+'SelFunc_nobao7_errP_nofz.txt']
sd = [dirg+'SelFunc_errPpodds_nofz.txt',dirng+'SelFunc_nobao6_errPpodds_nofz.txt',dirng+'SelFunc_nobao7_errPpodds_nofz.txt']


nsmooth=9

!p.multi=[0,2,2]
ns = 0
for i=0,2 do begin
   print,st[i],ns
   check = FILE_TEST(st[i])
   if (check eq 0) then continue
   print,'           read OK '
   ns = ns+1
   s = read_ascii(st[i], template =  TEMP_SEL_FUNC)
   z = s.(0)
   if (ns eq 1) then s1 = s.(1)
   if (ns eq 1) then specz = s.(1) else specz = specz + s.(1)
   if (ns eq 1) then plot,z,s.(1)/s1,/xs,/ys,ytit='Selection fct ratio',xra=[0.11,2.5],yra=[0.9,1.1],yma=[0,4],xma=[8,0] $
   else oplot,z,s.(1)/s1,col=80*i
endfor
xyouts,0.5,1.05,'spectroZ'
what=['case 2 / case 1','case 3 / case 1']
legend,what,col=[80,160],line=0,box=1,psym=0,/fill,/left,/bottom,charsize=1.5

specz = smooth(specz,nsmooth)/ ns

ns = 0
for i=0,2 do begin
   print,sg[i],ns
   check = FILE_TEST(sg[i])
   if (check eq 0) then continue
   print,'           read OK '
   ns = ns+1
   s = read_ascii(sg[i], template =  TEMP_SEL_FUNC)
   if (ns eq 1) then s1 = s.(1)
   if (ns eq 1) then gauss = s.(1) else gauss = gauss + s.(1)
   if (ns eq 1) then plot,z,s.(1)/s1,/xs,/ys,xra=[0.11,2.5],yra=[0.9,1.1],yma=[0,4],xma=[0,4] $
   else oplot,z,s.(1)/s1,col=80*i
endfor
xyouts,0.5,1.05,'Gaussian error: 0.03*(1+z)'
gauss = smooth(gauss,nsmooth) / ns

ns = 0
for i=0,2 do begin
   print,sp[i],ns
   check = FILE_TEST(sp[i])
   if (check eq 0) then continue
   print,'           read OK '
   ns = ns+1
   s = read_ascii(sp[i], template =  TEMP_SEL_FUNC)
   if (ns eq 1) then s1 = s.(1)
   if (ns eq 1) then photz  = s.(1) else photz = photz + s.(1)
   if (ns eq 1) then plot,z,s.(1)/s1,/xs,/ys,ytit='Selection fct ratio',xtit='redshift',xra=[0.11,2.5],yra=[0.9,1.1],yma=[4,0],xma=[8,0] $
   else oplot,z,s.(1)/s1,col=80*i
endfor
xyouts,0.5,1.05,'Photometric error'
photz = smooth(photz,nsmooth) / ns

ns = 0
for i=0,2 do begin
   print,sd[i],ns
   check = FILE_TEST(sd[i])
   if (check eq 0) then continue
   print,'           read OK '
   ns = ns+1
   s = read_ascii(sd[i], template =  TEMP_SEL_FUNC)
   if (ns eq 1) then s1 = s.(1)
   if (ns eq 1) then podds  = s.(1) else podds = podds + s.(1)
   if (ns eq 1) then plot,z,s.(1)/s1,/xs,/ys,xtit='redshift',xra=[0.11,2.5],yra=[0.9,1.1],yma=[4,0],xma=[0,4] $
   else oplot,z,s.(1)/s1,col=80*i
endfor
xyouts,0.5,1.05,'Photometric error + podds cut'
podds = smooth(podds,nsmooth)/ ns
read,xx
write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/sf_TJP_details.jpg' ,tvrd(true=3),true=3


window,0,xs=1000,ys=600
!p.multi=[0,1,2]
plot,z,specz,/xs,/ys,ytit='Mean selection fct',/yl,xra=[0.11,2.5],yra=[0.0002,1],yma=[0,4]
;oplot,z,gauss,col=80,li=2
oplot,z,photz,col=150
oplot,z,podds,col=210
;oplot,[2.2,2.2],[1e-5,1],th=1,col=250,li=2
;what=['spectroZ','Gaussian 0.03','photoZ','photoZ with podds cut']
;legend,what,col=[0,80,150,210],line=[0,1,0,0],box=1,psym=0,/fill,/left,/bottom,charsize=1.5
what=['spectroZ','photoZ','photoZ with podds cut']
legend,what,col=[0,150,210],line=0,box=1,psym=0,/fill,/left,/bottom,charsize=1.5

plot,z,specz/specz,/xs,/ys,ytit='Mean selection fct ratio',xtit='Redshift',xra=[0.11,2.5],yma=[4,0],yra=[0,2]
;oplot,z,gauss/specz,col=80,li=2
oplot,z,photz/specz,col=150
oplot,z,podds/specz,col=210

oplot,[2.2,2.2],[0,3],th=1,col=250,li=2

;stop
   
openw,lun0, '/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_speczmean_nofz.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun0, z[i], specz[i]
close, lun0
free_lun, lun0
   
openw,lun1, '/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_gaussmean_nofz.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun1, z[i], gauss[i]
close, lun1
free_lun, lun1
   
openw,lun2, '/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_photzmean_nofz.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun2, z[i], photz[i]
close, lun2
free_lun, lun2
   
openw,lun3, '/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_poddsmean_nofz.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun3, z[i], podds[i]
close, lun3
free_lun, lun3


write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/sf_TJP_mean.jpg' ,tvrd(true=3),true=3


restore,'temp_ascii.sav'
s0 = read_ascii('/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_speczmean_nofz.txt', template =  TEMP_SEL_FUNC)
s1 = read_ascii('/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_gaussmean_nofz.txt', template =  TEMP_SEL_FUNC)
s2 = read_ascii('/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_photzmean_nofz.txt', template =  TEMP_SEL_FUNC)
s3 = read_ascii('/sps/lsst/data/rcecile/TJP_BAO_grids/SelFunc_poddsmean_nofz.txt', template =  TEMP_SEL_FUNC)
plot,s0.(0),s0.(1),/xs,/ys,ytit='Mean election fct',/yl,xra=[0.11,2.5],yra=[0.0002,1]
oplot,s1.(0),s1.(1),col=50
oplot,s2.(0),s2.(1),col=100
oplot,s3.(0),s3.(1),col=150
end
