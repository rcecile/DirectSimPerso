PRO grid_theo

dir='/sps/lsst/data/rcecile/TJP_BAO_grids/'
suff='G2_'

D=['0.9','1.3','1.8','1.8']
snx=['640','900','1024','500']
nz=[125,75,65,75]
nD = n_elements(D)

for id = 0,nD-1 do begin
   grid=dir+'grid_'+suff+strtrim(snx[id],2)+'_z'+D[id]+'.fits'
   print,grid
   c=mrdfits(grid,2,h)
   read,xx
END
