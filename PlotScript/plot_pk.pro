PRO gfunct1, X, A, F, pder
  ax = -(X/0.1/0.679)^(1.4)
  bx = A[1]*X
  cx = X*EXP(ax)*SIN(bx)
  F = 1 + A[0]*cx
  
 IF N_PARAMS() GE 3 THEN pder= [[cx], [A[0]* X*EXP(ax)*COS(bx)*x]]

END

PRO gfunct2, X, A, F, pder
  ax = -(X/0.1)^(1.4)
  bx = A[1]*X
  cx = X*EXP(ax)*SIN(bx)
  F = 1 + A[0]*cx
  
 IF N_PARAMS() GE 3 THEN pder= [[cx], [A[0]* X*EXP(ax)*COS(bx)*x]]

END

PRO gfunct3, X, A, F, pder
  ax = -(X/0.1/0.679)^(1.4)
  bx = A[1]*X
  cx = sqrt(X)*EXP(ax)*SIN(bx)
  F = 1 + A[0]*cx
  
 IF N_PARAMS() GE 3 THEN pder= [[cx], [A[0]* X*EXP(ax)*COS(bx)*x]]

END

PRO gfunct4, X, A, F, pder
  ax = -(X/0.1)^(1.4)
  bx = A[1]*X
  cx = sqrt(X)*EXP(ax)*SIN(bx)
  F = 1 + A[0]*cx
  
 IF N_PARAMS() GE 3 THEN pder= [[cx], [A[0]* X*EXP(ax)*COS(bx)*x]]

END


PRO plot_pk

!p.charsize=2
!p.thick=3
loadct,12
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/jsouchar/Planck_BAO/"
dirn="/sps/lsst/data/jsouchar/Planck_noBAO/"

nline=[2,2,0,0];,0,0]
nz = n_elements(namez)
nk = n_elements(namek)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0

;Tables avec P(k) CLASS
   tCLASS_17  = dir + 'simu_ntpk_CLASS_17.txt' 
   print,tCLASS_17
   tCLASS_Pl15  = dir + 'simu_ntpk_CLASS_Pl15.txt' 
   print,tCLASS_Pl15
   tCLASS_17_Dmxtr  = dir + 'simu_ntpk_CLASS_17_DMxtr.txt' 
   print,tCLASS_17_DMxtr

;Tables avec P(k) Ein
   tEin_17  = dir + 'simu_ntpk_Ein_17.txt' 
   print,tEin_17
   tEin_Pl15  = dir + 'simu_ntpk_Ein_Pl15.txt' 
   print,tEin_Pl15
   tEin_17_DMxtr  = dir + 'simu_ntpk_Ein_17_DMxtr.txt' 
   print,tEin_17_DMxtr
   tEin_17_nobao  = dirn + 'simu6_ntpk.txt' 
   print,tEin_17_nobao

   pCLASS_17 = read_ascii(tCLASS_17, template =  TEMP_POW_SPEC_TH)
   pCLASS_Pl15 = read_ascii(tCLASS_Pl15, template =  TEMP_POW_SPEC_TH)
   pCLASS_17_DMxtr = read_ascii(tCLASS_17_DMxtr, template =  TEMP_POW_SPEC_TH)

   pEin_17 = read_ascii(tEin_17, template =  TEMP_POW_SPEC_TH)
   pEin_Pl15 = read_ascii(tEin_Pl15, template =  TEMP_POW_SPEC_TH)
   pEin_17_DMxtr = read_ascii(tEin_17_DMxtr, template =  TEMP_POW_SPEC_TH)
   pEin_17_nobao = read_ascii(tEin_17_nobao, template =  TEMP_POW_SPEC_TH)
   x=pCLASS_17.(0)

   
   ok=where(x ge 0.02 and x lt 0.2)
   xp=x[ok]

   pB1=(pCLASS_17.(2)-pEin_17_nobao.(2))/pEin_17_nobao.(2);/pEin_17.(1);-1.202e-3
   pB2=(pCLASS_17.(2)-pEin_17.(1))/pEin_17.(1);-1.634e-3
   pB3=(pCLASS_17.(2)-pCLASS_17.(1))/pCLASS_17.(1);+4.31e-4
   pB4=(pEin_17.(1)/pEin_17.(2)) 
   pB4p=pB4[ok]
   pB5=(pCLASS_17.(1)/pCLASS_17.(2))
   pB5p=pB5[ok]

   A4= [2.,160.]
   A5= [2.,160.]
   
   sig=dblarr(n_elements(xp))+0.01
   fit4 = curvefit( xp, pB4p, 1./sig, A4,sigma4, FUNCTION_NAME='gfunct1',/DOUBLE,chisq=c4,status=s4,iter=i4,tol=0.00001);,itmax=40,iter=i4,status=s4,tol=0.00001)
   fit5 = curvefit( xp, pB5p, 1./sig, A5,sigma5, FUNCTION_NAME='gfunct1',/DOUBLE,chisq=c5,status=s5,iter=i5,tol=0.00001);,itmax=40,iter=i5,status=s5,tol=1)
   print,'Paramètres fit Eisenstein : ',s4,c4,i4,A4
   print,'Paramètres fit CLASS : ',s5,c5,i5,A5
;   !p.multi=[0,1,2]
;   plot,xp,pB4p,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum [(Mpc)^-3]',TITLE='Plot4 : Fit Eisenstein',psym=1,/xs,/ys
;   oplot,xp,fit4,col=123
;   oPlot,[0,0.2],[1,1]
;   PLOT,xp,pB5p,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum [(Mpc)^-3]',TITLE='Plot5 : Fit CLASS',psym=4,/xs,/ys
;   oplot,xp,fit5,col=123
;   oPlot,[0,0.2],[1,1]

   !p.multi=0
   window,0
   PLOT,xp,pB5p,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum [(Mpc)^-3]',TITLE='Plot5 : Fit CLASS',psym=4,/xs,/ys,yra=[0.93,1.09]
   chi=dblarr(2)
   s=dblarr(2)
   fit5 = curvefit( xp, pB5p, 1./sig, A5,sigma5, FUNCTION_NAME='gfunct1',/DOUBLE,chisq=ch,status=s5,iter=i5,tol=0.00001);,itmax=40,iter=i5,status=s5,tol=1)
   oplot,xp,fit5,col=50
   chi[0]=ch
   s[0] = A5[1]
   fit5 = curvefit( xp, pB5p, 1./sig, A5,sigma5, FUNCTION_NAME='gfunct4',/DOUBLE,chisq=ch,status=s5,iter=i5,tol=0.00001);,itmax=40,iter=i5,status=s5,tol=1)
   oplot,xp,fit5,col=150
   chi[1]=ch
   s[1] = A5[1]
   oplot,[0,.2],[1,1],th=5,col=210,li=2
   legend,['k & h','sqrt(k)']+ ": chi2= " + strtrim(chi,2)+ ":  s= " + strtrim(s,2),col=[50,150],thick=3,lin=0,box=1,/fill,/right,/top,charsize=1.5
  print,chi

   window,1
   PLOT,xp,pB5p,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum [(Mpc)^-3]',TITLE='Plot5 : Fit CLASS',psym=4,/xs,/ys,yra=[0.93,1.09],xra=[0.03,0.05]
   chi=dblarr(2)
   s=dblarr(2)
   fit5 = curvefit( xp, pB5p, 1./sig, A5,sigma5, FUNCTION_NAME='gfunct1',/DOUBLE,chisq=ch,status=s5,iter=i5,tol=0.00001);,itmax=40,iter=i5,status=s5,tol=1)
   oplot,xp,fit5,col=50
   chi[0]=ch
   s[0] = A5[1]
   fit5 = curvefit( xp, pB5p, 1./sig, A5,sigma5, FUNCTION_NAME='gfunct4',/DOUBLE,chisq=ch,status=s5,iter=i5,tol=0.00001);,itmax=40,iter=i5,status=s5,tol=1)
   oplot,xp,fit5,col=150
   chi[1]=ch
   s[1] = A5[1]
   oplot,[0,.2],[1,1],th=5,col=210,li=2

   stop
end
