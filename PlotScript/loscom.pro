FUNCTION loscom,z

  c=3e5                         ; km/s
  H0=70.0                       ; km/s/Mpc
  OM = 0.30
  OL = 0.70
  ORad = 0.0000546405+0.0000377242
  
  zp =  1.+z

  RETURN, c/H0*1./sqrt(OL + (OM +  zp * ORad)* zp^3)

END
