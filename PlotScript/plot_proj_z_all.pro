PRO plot_proj_z_all,z,tot,nd,nErr,nz,snx,ngrids,saveplot

!p.multi=0
loadct,39
myzxt='Radial distance [redshift]'
mydxt='Radial comoving distance [Mpc]'
!x.title=myzxt
!p.charsize=2
err=['spectroZ','gauss 0.03','Photo-z','Photo-z PODDS']
err1=['spectroZ','gauss 0.03']
err2=['Photo-z','Photo-z PODDS']
mycol=80*(indgen(nErr))

   
   if (saveplot ) then begin
; current plotting device.
      mydevice = !D.NAME
      SET_PLOT, 'PS'
      DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/proj_z_grids.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=8.8,FONT_SIZE=5
   endif

!p.multi=[0,1,3]

for id = 0,nD-1 do begin


 ;  if (saveplot eq 0) then  window,id+20,xs=1000,ys=600
   snx2 = 1.*snx[id]*snx[id]
   plot,z[id,0:nz[id]-1],tot[0,id,0,0:nz[id]-1]/snx2,/ys,/nodata,ytit='Normalized Ngal / cell',xs=9 ,xtit=myzxt,yra=[0.97,1.07],xma=[7,3],yma=[3,3]
   for in=ngrids,0,-1 do begin
      for ie=0,nErr-1 do begin
         oplot,z[id,0:nz[id]-1],tot[in,id,ie,0:nz[id]-1]/snx2,col=mycol[ie],th=1,psym=10
      endfor
      ;if (saveplot eq 0) then read,xx
   endfor
   zn= findgen(6)*(max(z[id,0:nz[id]-1])-min(z[id,0:nz[id]-1]))/5+min(z[id,0:nz[id]-1])     
   zv= dloscom(zn)              
   AXIS, XAXIS=1, XSTYLE = 1,  XTITLE = mydxt, xTICKV = zn,XTicks=9,xtickn=string(fix(zv),format='(I5)')
;   if(id eq 0) then legend,err,col=mycol,box=1,psym=0,/fill,/right,/top,charsize=1.5,th=3 $
;    legend,err,col=mycol,box=1,psym=0,/fill,/left,/top,charsize=1.5,th=2
   
   if(id eq 0) then legend,err[0],col=mycol[0],box=1,psym=0,/fill,/left,/top,charsize=1.5,th=3
   if(id eq 0) then legend,err[1],col=mycol[1],box=1,psym=0,/fill,/right,/top,charsize=1.5,th=3
   if(id eq 1) then legend,err[2],col=mycol[2],box=1,psym=0,/fill,/left,/top,charsize=1.5,th=3
   if(id eq 1) then legend,err[3],col=mycol[3],box=1,psym=0,/fill,/right,/top,charsize=1.5,th=3


endfor
   if (saveplot) then begin
      DEVICE, /CLOSE
      SET_PLOT, mydevice
   endif

stop

end
