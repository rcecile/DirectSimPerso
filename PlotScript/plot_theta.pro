PRO plot_theta

loadct,39
!p.thick=3


theta= dindgen(6000)/100.
rtheta = theta *!pi/180.

ht = dindgen(6000,9)
z = dindgen(9)
for slice=5,90,10 do begin
   m=mrdfits("/sps/lsst/data/rcecile/TJP_BAO/cat_G2_Slice" + strtrim(slice,2)+".fits",1,head)
   m=mrdfits("/sps/lsst/data/rcecile/Planck_BAO/cat_zOrd_Slice" + strtrim(slice,2)+".fits",1,head)
   stheta = m.(2)
   print,slice,mean(m.(3))
   h = histogram(stheta,min=0,max=max(rtheta),bin=0.01*!pi/180.)
   for i=0,5999 do ht[i,slice/10] = total(h[0:i])
   z[slice/10] = mean(m.(3))

endfor


window,2
plot,theta,[1,.35e8],/nodata,/xs,/ys,tit='Ngal avec theta < Theta',xtit='Theta [degree]',ytit='Ngal';,xra=[0,10]
for i=0,8 do   oplot,theta,ht[*,i],col=i*50
write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/Ngal_theta.jpg' ,tvrd(true=3),true=3


window,4
plot,theta,[0,.7e8],/nodata,/xs,/ys,tit='Ngal / tan^2 avec theta < Theta',xtit='Theta';,xra=[0,10]
for i=0,8 do   oplot,theta,ht[*,i]/tan(rtheta)/tan(rtheta),col=i*50
;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/geometrie_theta.jpg' ,tvrd(true=3),true=3


window,5
plot,theta,[2e3,.7e8],/nodata,/xs,/ys,tit='Ngal / (2x(1-cos(Theta))) avec theta < Theta',xtit='Theta [degree]',/yl;,xra=[0,10]
for i=0,7 do   oplot,theta,ht[*,i]/(2.*(1-cos(rtheta))),col=i*50
legend,'z = '+ strtrim(z,2),col=indgen(9)*50,li=0,box=1,/fill,/right,/bottom,charsize=1.5

write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/geometrie_theta.jpg' ,tvrd(true=3),true=3

window,3
plot,theta,2.*(1-cos(rtheta)),yra=[1e-3,1],tit='noir = surface coquille, vert = surface disque, meme redshift',xtit='Theta'
oplot,theta,tan(rtheta)*tan(rtheta),col=123
;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/Ngal_geometrie.jpg' ,tvrd(true=3),true=3

stop


END
