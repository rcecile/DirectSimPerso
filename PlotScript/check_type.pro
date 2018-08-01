PRO check_type,is

restore,'temp_ascii_new.sav'
!p.charsize=1.5
!p.multi=0
!p.thick=3

!p.multi=0

   m0=mrdfits("/sps/lsst/data/rcecile/Planck_BAO2/cat_All_Slice"+strtrim(is,2)+".fits",1)
   m1=mrdfits("/sps/lsst/data/rcecile/Planck_BAO2/cat_ELSB_Slice"+strtrim(is,2)+".fits",1)
   m3=mrdfits("/sps/lsst/data/rcecile/Planck_BAO2/cat_LFdahlen_mag24_13all_Slice"+strtrim(is,2)+".fits",1)
   m4=mrdfits("/sps/lsst/data/rcecile/Planck_BAO2/cat_LFdahlen_mag24_13_Slice"+strtrim(is,2)+".fits",1)

   binsize=0.05
   h0=histogram(m0.(4),min=0.9,max=3.6,bin=binsize)
   h1=histogram(m1.(4),min=0.9,max=3.6,bin=binsize)
   h3=histogram(m3.(4),min=0.9,max=3.6,bin=binsize)
   h4=histogram(m4.(4),min=0.9,max=3.6,bin=binsize)

   x=findgen((3.6-0.9)/binsize+1)*binsize+0.9
   plot,x,h0,psym=10,/xs,/ys,th=3,yra=[1,max([max(h0),max(h3)])*1.1],tit='Slice '+strtrim(is,2)+' - z in ['+strtrim(min(m0.(3)),2)+','+strtrim(max(m0.(3)),2)+']'
   oplot,x,h0,col=123,psym=10
   oplot,x,h1,col=23,psym=10
   oplot,x,h3,col=123,psym=10,li=2
   oplot,x,h4,col=23,psym=10,li=2
   legend,['old all=true','old all=false' ,'new all=true','new all=false'],col=[123,23,123,23],li=[0,0,2,2],/fill,/left,/top
   ntot0= total(h0)/100.
   ntot1= total(h1)/100.
   ntot4= total(h4)/100.
   ntot3= total(h3)/100.
   print,' h0 type ',total(h0[0:21]),total(h0[22:41]),total(h0[42:53]),' et en % ',total(h0[0:21])/ntot0,total(h0[22:41])/ntot0,total(h0[42:53]) /ntot0
   print,' h1 type ',total(h1[0:21]),total(h1[22:41]),total(h1[42:53]),' et en % ',total(h1[0:21])/ntot1,total(h1[22:41])/ntot1,total(h1[42:53]) /ntot1
   print,' h3 type ',total(h3[0:21]),total(h3[22:41]),total(h3[42:53]),' et en % ',total(h3[0:21])/ntot3,total(h3[22:41])/ntot3,total(h3[42:53]) /ntot3
   print,' h4 type ',total(h4[0:21]),total(h4[22:41]),total(h4[42:53]),' et en % ',total(h4[0:21])/ntot4,total(h4[22:41])/ntot4,total(h4[42:53]) /ntot4

   print, ' TOTAL ',total(h0),total(h1),total(h3),total(h4)

   write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/check_type'+strtrim(is,2)+'.jpg' ,tvrd(true=3),true=3
   stop
END
