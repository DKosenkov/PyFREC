import math
#-------------------------------------------------------------------------------
# Conversion factors
ATOB          = 1.889725989 #Angstrom to Bohr Conversion
BTOA          = 1.0/ATOB
HartreeToKCal = 627.509
HartreeToCM1  = 219474.629232
HartreeToeV   = 27.211215
eVToHartree   = 1/HartreeToeV
#eVToCM1_old       = 8065.54429 #web
#eVToCM1       =  8065.54400479571 #our
TimeAuToS     = 2.418884326505E-17 #Atomic unit of time to seconds

qe = 1.6021766208E-19 #Charge of an electron / coulombs
c  = 299792458  #Speed of light / m/s
h  = 6.626070040E-34 #Planck constant  Jxs
eVToCM1 = qe/(100*h*c)

c_au = 137.035999139 #Speed of light 1/alpha = 1/(Fine-structure constant) in atomic units
Bohr_Radius_m = 5.2917721092E-11 # Bohr radiaus in meters

TDM_auToDebye = 2.541746 # (Transition) Dipole Moment atomic units (a.u. or ea0) to Debye
#1 Debye = 0.393430307 a.u.
#1 TDM a.u. =  2.541746 Debye

GAS_CONST_CM1_K = 0.69465783 #Gas constant in cm-1/K

#Quantum dynamics:
cm1_to_ps1_by_hbar1 = 2*1.0e2*math.pi*c*1.0e-12  #approximate value of the conversion factor: 0.19  (1 cm-1 = ps-1/h_bar)
cm1_to_ps1 = 1.0e2*c*1.0e-12  #approximate value of the conversion factor: 0.029 (1 cm-1 = 0.03 ps)

#Conversion of GROMACS GRO files in nm to ang 
nm_to_ang = 10.0

#Conversion factor   CM-1 = CM1_NM/nm; nm = CM1_NM/cm-1
CM1_NM = 1.0E7

AVOGADROS_Na = 6.0221415E23 #Avogadro's constant per mole
#-------------------------------------------------------------------------------
ELEMENTS_BY_ATOMIC_N = ['X','H','He',
'Li','Be','B','C','N','O','F','Ne',
'Na','Mg','Al','Si','P','S','Cl','Ar',
'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',
'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',
'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Uut','Fl','Uup','Lv','Uus','Uuo']
#-------------------------------------------------------------------------------
RE_FLOAT = '[-]?\d+(?:.\d+)?|\.\d+' #Regular expression to match float numbers
RE_DOUBLE = '[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eED][+\-]?\d+)?'