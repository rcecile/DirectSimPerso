PRO create_temp

restore,'temp_ascii_new.sav'

t = '/sps/lsst/data/rcecile/Planck_BAO2_PS/PS_SN_120_z0.5_G0_wngal.txt'
TEMP_POW_SPEC_11col = ASCII_TEMPLATE(t)

t = '/sps/lsst/data/rcecile/Planck_BAO2_PS/PS_fromgrids_120_z0.5_G0_wngal.txt'
TEMP_POW_SPEC_4col = ASCII_TEMPLATE(t)

t = '/sps/lsst/data/rcecile/Planck_BAO2_PS/simu_ps_120_z0.51_ntpk.txt'
TEMP_POW_SPEC_TH = ASCII_TEMPLATE(t)

t = '/sps/lsst/data/rcecile/Planck_BAO2_PS/shotnoiseErr_3z_5err.txt'
TEMP_SHOTNOISE = ASCII_TEMPLATE(t)

t = '/sps/lsst/data/rcecile/Planck_BAO2_PS/PS_nCR_120_z0.5_wngal.txt'
TEMP_POW_SPEC_save = ASCII_TEMPLATE(t)

t = '/sps/lsst/data/rcecile/Planck_BAO2_PS/fit_superSN5_k0.04_0.10_320_z1.5_baseline.txt'
TEMP_POW_SPEC_BASE = ASCII_TEMPLATE(t)

t = '/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_All__ZONLY.txt'
TEMP_SEL_FUNC = ASCII_TEMPLATE(t)

t = '/sps/lsst/data/rcecile/Planck_BAO2_PS/fit_superSN5_k0.04_0.20_120_z0.5_result.txt'
TEMP_POW_SPEC_FIT_RES = ASCII_TEMPLATE(t)

save,TEMP_POW_SPEC_FIT_RES,TEMP_POW_SPEC_11col,TEMP_POW_SPEC_5col,TEMP_POW_SPEC_4col,TEMP_POW_SPEC_TH,TEMP_SHOTNOISE,TEMP_POW_SPEC_save,TEMP_POW_SPEC_BASE,TEMP_SEL_FUNC,file='temp_ascii_new.sav'


END

