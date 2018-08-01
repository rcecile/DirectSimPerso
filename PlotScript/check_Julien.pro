PRO check_Julien,iz,isuff

!p.charsize=2.5
!p.thick=3
!p.symsize=3
h= 0.679
h3 = h*h*h

; f='/sps/lsst/data/rcecile/Planck_BAO_PS/PS_obs_Ref_120_z0.5_err0.03_wngal.txt'
; TEMP_PS_OBS = ASCII_TEMPLATE(f)
; save,TEMP_2FIT,TEMP_INFO,TEMP_INFOS,TEMP_PDF,TEMP_POW_SPEC,TEMP_POW_SPEC_FIT,TEMP_POW_SPEC_FITCHI2,TEMP_POW_SPEC_TH,TEMP_POW_SPEC_TXT,TEMP_POW_SPEC_ZXY,TEMP_SEL_FUNC,TEMP_SHOTNOISE,TEMP_PS_OBS ,file='temp_ascii.sav'

namez = ['_z0.5','_z0.9','_z1.5']
nx= ['_120','_225', '_320']

nameerr = ['', '_err0.03', '_errP','_errPBDT9','_errPBDT8']
lcol  = [95,210,35,135,  120]

loadct,12
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO_PS/"

t=dir+'PS_2fit_SN_Ref'+nx[iz]+namez[iz]+nameerr[isuff]+'_wngal.txt'
tt= dir+'PS_obs_Ref'+nx[iz]+namez[iz]+nameerr[isuff]+'_wngal.txt'         

pno = read_ascii(t, template =  TEMP_2FIT)     
pwi = read_ascii(tt, template =  TEMP_PS_OBS)     

k = pno.(0)
obs = pwi.(1)
ref = pno.(8)

plot,k,obs*k*k,/xs,/ys,xra=[0.015,.2]
oplot,k,ref*k*k,col=123

stop
END
