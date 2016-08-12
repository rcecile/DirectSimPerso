PRO error_distrib_plot

loadct,39
!p.charsize=3
!p.thick=3
lcol  = [150, 240, 210]

restore,'/sps/lsst/data/rcecile/TJP_BAO/error_distrib.sav'
suff=["_err0.03","_errP","_errPpodds"]
nsuff = n_elements(suff)
hmin=-0.09
hmax=0.09
hbin=0.001
x=findgen((hmax-hmin)/hbin+1)*hbin+hmin
nx=n_elements(x)

ymar0=[0,0,0,0,4,4]
ymar1=[1,1,0,0,0,0]
xmar0=[8,0,8,0,8,0]
xmar1=[0,4,0,4,0,4]
mytit=['Nb of galaxies','','Nb of galaxies','','Nb of galaxies','']
mxtit=['','','','','Normalized redshift','Normalized redshift']

!p.multi=[0,2,3]

for i=0,5 do begin
   plot,x,[0,0],/xs,/ys,ytit=mytit[i],xtit=mxtit[i],xra=[hmin,hmax],yra=[0,1.1],/nodata,$
        yma=[ymar0[i],ymar1[i]],xma=[xmar0[i],xmar1[i]]
   for is=0,nsuff-1 do  oplot,x,hh[*,is,i],col=lcol[is]
   oplot,[hmin,hmax],[0.5,0.5],th=1,col=200,li=2
   oplot,[0,0],[0,2],th=1,col=200,li=2
endfor

END
