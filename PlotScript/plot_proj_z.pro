PRO plot_proj_z,z,tot,nd,nErr,nz,snx,saveplot

!p.multi=0
loadct,39
myzxt='Radial distance [redshift]'
mydxt='Radial comoving distance [Mpc]'
!x.title=myzxt
!p.charsize=2
err=['spectroZ','gauss 0.03','Photo-z','Photo-z PODDS']
mycol=80*(indgen(nErr))


for id = 0,nD-1 do begin
   window,id+20,xs=1000,ys=600
   snx2 = 1.*snx[id]*snx[id]
   plot,z[id,0:nz[id]-1],tot[id,0,0:nz[id]-1]/snx2,/ys,th=3,/nodata,ytit='Ngal / cell',xs=9 ,yma=[4,4],xtit=myzxt,yra=[0.96,1.06]
   for ie=0,nErr-1 do begin
      oplot,z[id,0:nz[id]-1],tot[id,ie,0:nz[id]-1]/snx2,col=mycol[ie],th=3,psym=10
   endfor
   
   zn= findgen(6)*(max(z[id,0:nz[id]-1])-min(z[id,0:nz[id]-1]))/5+min(z[id,0:nz[id]-1])     
   zv= dloscom(zn)              
   AXIS, XAXIS=1, XSTYLE = 1,  XTITLE = mydxt, xTICKV = zn,XTicks=9,xtickn=string(fix(zv),format='(I5)')
   legend,err,col=mycol,box=1,psym=0,/fill,/left,/top,charsize=2,th=3


endfor


end
