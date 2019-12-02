#---------------------------------------------------------------------------
# transmission.py
#
# (c) 30.1.2004 Ennes Sarradj
#
# transmission loss calculation
#---------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#  Imports:
#-------------------------------------------------------------------------------

#~ try:
    #~ from wxPython.wx import *
#~ except:
    #~ pass
from traits.api import *
#from traits.trait_handlers import TraitHandler
from numpy import *
from scipy.linalg import solve
#from LinearAlgebra import *
#from sys import maxint
from .tr_traits_helper import PosFloat
from .itertype import itertype
from .gauss import *

#-------------------------------------------------------------------------------
#  base class for all objects; does nothing (at the moment)
#-------------------------------------------------------------------------------
class TrObj(HasTraits):
    pass
    
#-------------------------------------------------------------------------------
#  base class for all materials
#-------------------------------------------------------------------------------
class Material(TrObj):
    
#    __traits__={
#        'dens': Trait(7800.0,PosFloat()) # density
#    }
    dens = Trait(7800.0,PosFloat())    
#-------------------------------------------------------------------------------
#  elastic material
#-------------------------------------------------------------------------------
class ElasticMaterial(Material):
    
#    __traits__={
#        'poisson': Trait(0.3,TraitRange(0.0,0.5)), # Poisson number     
#        'elas': Trait(2.1e11,PosFloat()) # Young's modulus
#    }
    poisson = Range(0.0,0.5,0.3)
    elas = Trait(2.1e11,PosFloat())    
    #---------------------------------------------------------------------------
    #  shear modulus
    #---------------------------------------------------------------------------
    def G(self):
        return self.elas/(2*(1+self.poisson))
    
#-------------------------------------------------------------------------------
#  base class for plate sections; useful?? (there may but one child class)
#-------------------------------------------------------------------------------
class PlateSec(TrObj):
    pass
    
#-------------------------------------------------------------------------------
#  isotropic plate section
#-------------------------------------------------------------------------------
class IsoPlateSec(PlateSec):
    
#    __traits__={
#        'thick': Trait(5e-3,PosFloat()), # plate thickness
#        'mat': Trait(ElasticMaterial(),
#                    TraitInstance(ElasticMaterial))#, # plate material
#    }
    thick = Trait(5e-3,PosFloat())
    mat = Trait(ElasticMaterial)
    Zinf = zeros((4,4,1),complex)
    k = zeros(2,float)
    om = 0.0
    #---------------------------------------------------------------------------
    #  set and define change handler 
    #---------------------------------------------------------------------------
#    def __init__(self,*kw,**args):
#        PlateSec.__init__(self,*kw,**args)
#        self.Zinf=zeros((4,4,1),complex)
#        self.k=0
#        self.om=0

    #---------------------------------------------------------------------------
    #  moment of inertia per width
    #---------------------------------------------------------------------------
    def I(self):
        return self.thick**3/12

    #---------------------------------------------------------------------------
    #  mass per area
    #---------------------------------------------------------------------------
    def m2(self):
        return self.mat.dens*self.thick

    #---------------------------------------------------------------------------
    #  bending stiffness per width
    #---------------------------------------------------------------------------
    def B1(self):
        return self.mat.elas*self.I()/(1-self.mat.poisson**2)

    #---------------------------------------------------------------------------
    #  bending wave number^4
    #---------------------------------------------------------------------------
    def kb4(self,om):
        return om**2*self.m2()/self.B1()

    #---------------------------------------------------------------------------
    #  compression wave number^2
    #---------------------------------------------------------------------------
    def kc2(self,om):
        return om**2*self.mat.dens*(1-self.mat.poisson**2)/self.mat.elas

    #---------------------------------------------------------------------------
    #  shear wave number^2
    #---------------------------------------------------------------------------
    def ks2(self,om):
        return om**2*self.mat.dens/self.mat.G()
        
    #---------------------------------------------------------------------------
    #  calculate half-infinite plate values:
    #  input: om - omega,
    #         k - vector of scalar, wave number
    #  output: (all with extra dimension for the wavenumber k)
    #   Z - wave dynamic stiffness matrix of infinite plate
    #   Vca,Vsa - matrices for extraction of compression and shear wave amplitude from
    #           velocities (the latter is not needed)
    #   Ca - vectors of incident wave velocities
    #   P - incident wave powers
    #   Ca columns, P rows: 
    #           0 - compression
    #           1 - shear
    #           2 - bending
    #---------------------------------------------------------------------------
    def calc_inf(self,om,k):
        #if self.om==om and self.k==k:
        #    return
        self.om=om
        self.k=k
        if self.Zinf.shape[2]!=len(k):
            self.Zinf=zeros((4,4,len(k)),complex)
        k2=k*k
        om2=om*om
        # bending impedance:
        kb4=self.kb4(om)
        kb2=sqrt(kb4)
        gb=sqrt(k2-complex(kb2))
        ge=sqrt(k2+kb2)
        hb=om2*self.m2()/kb4
        hb1=hb*(gb+ge)
        hb2=hb*(gb*ge+self.mat.poisson*k2)
        self.Zinf[2:4,2:4]=array((( gb*ge*hb1 , hb2 ),
                                  ( hb2       , hb1 )))
        #compression and shear impedance:
        kc2=self.kc2(om)
        ks2=self.ks2(om)
        gc=sqrt(k2-complex(kc2))
        gs=sqrt(k2-complex(ks2))
#        print gb,ge,gc,gs
        hcs=1/(k2-gc*gs)
        hcs1=om2*self.m2()*hcs
        hcs2=hcs1*k2*(2*k2-ks2-2*gc*gs)/ks2
        self.Zinf[0:2,0:2]=array((( hcs1*k2*gc , hcs2 ),
                                  ( hcs2       , hcs1*gs )))
        # matrices for extraction of compression and shear amplitudes
        self.Vca=-hcs*array(((-k2    , gs    ),
                            (-gc*k2 , gc*gs )))
        #self.Vsa=-hcs*array((( gc*gs ,-gs    ),
        #                    ( gc*k2 ,-k2    )))
        # [Ccse*vc0, Csce*vs0, Cbe*vb0]
        one=ones((len(k)),complex)
        self.Ca=zeros((4,3,len(k)),complex)
        self.Ca[0:2,0:2]=array((( one , gs/k ),
                                ( -gc  , -k    ))) 
        self.Ca[2:4,2]=array((one,gb))
        self.P=abs(om*self.m2()*array((gc/2,gs/2,gb/kb2)))
        
    #---------------------------------------------------------------------------
    #  calculate values for strip plate impedance for compression/shear:
    #  z1 - z5,P,Q 
    #---------------------------------------------------------------------------
    def strip_cs(self,om,k,b):
        k2=k*k
        S=self.mat.G()*self.thick
        gc=sqrt(k2-complex(self.kc2(om)))
        gs=sqrt(k2-complex(self.ks2(om)))
        z1=gs/k
        z2=gc
        z3=k
        z4=2*k2*gc*S
        z5=k*(gs*gs+k2)*S
        P=exp(-gc*b)
        Q=exp(-gs*b)
        return (z1,z2,z3,z4,z5,P,Q)
        
    #---------------------------------------------------------------------------
    #  calculate values for strip plate impedance for bending:
    #  z1 - z5,P,Q
    #---------------------------------------------------------------------------
    def strip_b(self,om,k,b):
        k2=k*k
        kb2=sqrt(self.kb4(om))
        gb=sqrt(k2-complex(kb2))
        ge=sqrt(k2+kb2)
        B1=self.B1()
        muk2=self.mat.poisson*k2
        z1=ones((len(k)),Float)
        z2=-gb
        z3=-ge
        z4=gb*(ge*ge-muk2)*B1
        z5=ge*(gb*gb-muk2)*B1
        P=exp(-gb*b)
        Q=exp(-ge*b)
        return (z1,z2,z3,z4,z5,P,Q)

#-------------------------------------------------------------------------------
#  base class for beam sections
#-------------------------------------------------------------------------------
class BeamSec(TrObj):

#    __traits__={
#       'mat': Trait(ElasticMaterial(),
#                TraitInstance(ElasticMaterial)) # beam material
#    }
    mat = Trait(ElasticMaterial)
    
    #---------------------------------------------------------------------------
    #  shear centre with respect to centroid
    #---------------------------------------------------------------------------
    def SC(self):
        return (0.0,0.0)
        
    #---------------------------------------------------------------------------
    #  cross sectional area
    #---------------------------------------------------------------------------
    def A(self):
        return 1.0
        
    #---------------------------------------------------------------------------
    #  mass per length
    #---------------------------------------------------------------------------
    def m1(self):
        return 1.0
        
    #---------------------------------------------------------------------------
    #  moment of inertia about y
    #---------------------------------------------------------------------------
    def Iy(self):
        return 1.0
        
    #---------------------------------------------------------------------------
    #  moment of inertia about z
    #---------------------------------------------------------------------------
    def Iz(self):
        return 1.0
    
    #---------------------------------------------------------------------------
    #  polar moment of inertia
    #---------------------------------------------------------------------------
    def Ip(self):
        return self.Iy()+self.Iz()
        
    #---------------------------------------------------------------------------
    #  longitudinal stiffness
    #---------------------------------------------------------------------------
    def EA(self):
        return self.mat.elas*self.A()

    #---------------------------------------------------------------------------
    #  shear stiffness
    #---------------------------------------------------------------------------
    def GA(self):
        return self.mat.G()*self.A()

    #---------------------------------------------------------------------------
    #  torsion stiffness
    #---------------------------------------------------------------------------
    def GJt(self):
        return 1.0

    #---------------------------------------------------------------------------
    #  Bending stiffness about y
    #---------------------------------------------------------------------------
    def By(self):
        return 1.0
        
    #---------------------------------------------------------------------------
    #  Bending stiffness about z
    #---------------------------------------------------------------------------
    def Bz(self):
        return 1.0
        
    #---------------------------------------------------------------------------
    #  shear coefficient y
    #---------------------------------------------------------------------------
    def kappay(self):
        return 1.2
        
    #---------------------------------------------------------------------------
    #  shear coefficient z
    #---------------------------------------------------------------------------
    def kappaz(self):
        return 1.2
        
    #---------------------------------------------------------------------------
    #  Line wave impedance / (I om)
    #---------------------------------------------------------------------------
    def Z(self,om,k):
        om2=om*om
        om2k=om2*k
        k2=k*k
        om2m1=om2*self.m1()
        h=k2*self.By()-om2*self.Iy()*self.mat.dens
        k2BBy=h/(1+self.kappay()*h/self.GA())
        h=k2*self.Bz()-om2*self.Iz()*self.mat.dens
        k2BBz=h/(1+self.kappaz()*h/self.GA())
        #~ print k2BBy,k2*self.By()
        #
        k2BBy=k2*self.By()
        k2BBz=k2*self.Bz()
        (ey,ez)=self.SC()
        Z=zeros((4,4,len(k)),Float)
        # L 
        Z[0,0]=k2*(k2*self.EA()-om2m1)
        # Fz
        Z[1,1]=k2*k2BBz-om2m1
        # Fy
        Z[2,2]=k2*k2BBy-om2m1
        # Tshear
        Z[3,3]=k2*self.GJt()-om2*(self.Ip()*self.mat.dens+self.m1()*(ey*ey+ez*ez))
        Z[1,3]=-ez*om2m1
        Z[2,3]=ey*om2m1
        Z[3,1]=Z[1,3]
        Z[3,2]=Z[2,3]
        return Z
    
#-------------------------------------------------------------------------------
#  rectangular beam section
#-------------------------------------------------------------------------------
class RectBeamSec(BeamSec):

#    __traits__={
#        'leny': Trait(1e-2,PosFloat()), # height of beam section (y-direction)
#        'lenz': Trait(1e-2,PosFloat()) # width of beam section (z-direction)
#    }
    leny = Trait(1e-2,PosFloat())
    lenz = Trait(1e-2,PosFloat())
    
    #---------------------------------------------------------------------------
    #  cross sectional area
    #---------------------------------------------------------------------------
    def A(self):
        return self.leny*self.lenz
        
    #---------------------------------------------------------------------------
    #  mass per length
    #---------------------------------------------------------------------------
    def m1(self):
        return self.mat.dens*self.leny*self.lenz
        
    #---------------------------------------------------------------------------
    #  moment of inertia about y
    #---------------------------------------------------------------------------
    def Iy(self):
        return (self.lenz*self.leny**3)/12
        
    #---------------------------------------------------------------------------
    #  moment of inertia about z
    #---------------------------------------------------------------------------
    def Iz(self):
        return (self.leny*self.lenz**3)/12
    
    #---------------------------------------------------------------------------
    #  torsion stiffness
    #---------------------------------------------------------------------------
    def GJt(self):
        q=self.leny/self.lenz
        if (q<1):
            q=1/q
        return self.mat.G()*((self.leny*self.lenz)**2)/(q*(3+2.12002/q**2.5522+1.9738/q**0.958992))
        
    #---------------------------------------------------------------------------
    #  Bending stiffness about y
    #---------------------------------------------------------------------------
    def By(self):
        return self.mat.elas*self.Iy()
        
    #---------------------------------------------------------------------------
    #  Bending stiffness about z
    #---------------------------------------------------------------------------
    def Bz(self):
        return self.mat.elas*self.Iz()
        
#-------------------------------------------------------------------------------
#  connection node (for both 2D and 3D)
#-------------------------------------------------------------------------------
class CNode(TrObj):
#    __traits__={
#        'x': 0.0, # location
#        'y': 0.0,
#        'z': 0.0
#    }
    x = 0.0
    y = 0.0
    z = 0.0

    
#-------------------------------------------------------------------------------
#  base class for connection members
#-------------------------------------------------------------------------------
class ConnMember(TrObj):
    
    K = Property()

    #---------------------------------------------------------------------------
    #  returns a tuple of nodes (virtual)
    #---------------------------------------------------------------------------
    def nodes(self):
        return []

    #---------------------------------------------------------------------------
    #  is the member active, i.e. carries in/outgoing waves (virtual)
    #---------------------------------------------------------------------------
    def active(self):
        return False

    def _get_K(self):
        return self.K_calc()
    #---------------------------------------------------------------------------
    #  if K was never calculated, we should call K_calc first
    #---------------------------------------------------------------------------
#    def __getattr__(self, name):
#        if name=="K":
#            self.K_calc()
#            return self.__dict__["K"]
#        else:
#            return self.__dict__[name]#HasTraits.__getattribute__(self, name)

    #---------------------------------------------------------------------------
    #  "virtual" function, will raise exception
    #---------------------------------------------------------------------------
    def K_calc(self):
        raise AttributeError

#-------------------------------------------------------------------------------
#  infinite plate connection member
#-------------------------------------------------------------------------------
class InfPlate(ConnMember):
    
#    __traits__={
#        'sec': Trait(IsoPlateSec()), #section of plate
#        'node': Trait(CNode()), #attach node
#        'y0': 0.0, # y-offset
#        'z0': 0.0, # z-offset
#        'theta': 1.0 # connection angle in degree
#    }
    sec = Trait(IsoPlateSec()) #section of plate
    node = Trait(CNode()) #attach node
    y0 = 0.0 # y-offset
    z0 = 0.0 # z-offset
    theta = 1.0 # connection angle in degree

       
    #---------------------------------------------------------------------------
    #  returns a tuple of nodes of length 1
    #---------------------------------------------------------------------------
    def nodes(self):
        return [self.node]

    #---------------------------------------------------------------------------
    #  InfPlate is active, i.e. carries in/outgoing waves
    #---------------------------------------------------------------------------
    def active(self):
        return True

    #---------------------------------------------------------------------------
    #  K matrix co-ordinate transformation, must be called if 
    #  co-ordinate attributes have changed
    #---------------------------------------------------------------------------
    def K_calc(self):
        c=cos(self.theta*pi/180)
        s=sin(self.theta*pi/180)
        y=self.y0+self.node.y
        z=self.z0+self.node.z
        #self.K=
        return array(((1.,0,0,0),(0,c,-s,0),(0,s,c,0),(0,y*s-z*c,y*c+z*s,1.)))
    
#-------------------------------------------------------------------------------
#  strip plate connection member
#-------------------------------------------------------------------------------
class StripPlate(ConnMember):
    
#    __traits__={
#        'sec': Trait(IsoPlateSec()), #section of plate
#        'node1': Trait(CNode()), #attach node 1
#        'node2': Trait(CNode()), #attach node 2
#        'y01': 0.0, # y-offset at node 1
#        'z01': 0.0, # z-offset at node 1
#        'y02': 0.0, # y-offset at node 2
#        'z02': 0.0 # z-offset at node 2
#    }   
    sec = Trait(IsoPlateSec()), #section of plate
    node1 = Trait(CNode()) #attach node 1
    node2 = Trait(CNode()) #attach node 2
    y01 = Float(0.0) # y-offset at node 1
    z01 = 0.0 # z-offset at node 1
    y02 = 0.0 # y-offset at node 2
    z02 = 0.0 # z-offset at node 2
        
    #~ #---------------------------------------------------------------------------
    #~ #  if K was never calculated, we should call K_calc first
    #~ #---------------------------------------------------------------------------
    #~ def __getattr__(self, name):
        #~ if name=="K":
            #~ self.K_calc()
            #~ return self.__dict__["K"]
        #~ else:
            #~ return HasTraits.__getattr__(self, name)
        
    #---------------------------------------------------------------------------
    #  returns a tuple of nodes of length 2
    #---------------------------------------------------------------------------
    def nodes(self):
        return [self.node1,self.node2]

    #---------------------------------------------------------------------------
    #  StripPlate is not active, i.e. carries no in/outgoing waves
    #---------------------------------------------------------------------------
    def active(self):
        return False

    #---------------------------------------------------------------------------
    #  K matrix co-ordinate transformation, must be called if 
    #  co-ordinate attributes have changed
    #---------------------------------------------------------------------------
    def K_calc(self):#,x1,x2,x3):
        y1=self.y01+self.node1.y
        z1=self.z01+self.node1.z
        y2=self.y02+self.node2.y
        z2=self.z02+self.node2.z
        b=self.width()
        #c=(self.node2.y-self.node1.y)/b
        c=(y2-y1)/b
        #s=(self.node2.z-self.node1.z)/b
        s=(z2-z1)/b
        K=zeros((8,8),Float)
        K[0:4,0:4]=array(((1.,0,0,0),(0,c,-s,0),(0,s,c,0),(0,y1*s-z1*c,y1*c+z1*s,1.)))
        K[4:8,4:8]=array(((1.,0,0,0),(0,c,-s,0),(0,s,c,0),(0,y2*s-z2*c,y2*c+z2*s,1.)))
        return K
    #---------------------------------------------------------------------------
    #  strip plate width
    #---------------------------------------------------------------------------
    def width(self):
        yl2=(self.node1.y+self.y01-self.node2.y-self.y02)**2
        zl2=(self.node1.z+self.z01-self.node2.z-self.z02)**2
        b= (yl2+zl2)**0.5
        if b==0:
            raise ValueError("strip plate must not have zero width")
        return b

    #---------------------------------------------------------------------------
    #  strip plate cross line wave impedance
    #---------------------------------------------------------------------------
    def Zstrip(self,om,k):
        ZZ=zeros((8,8,len(k)),Float)
        self.ZxfmCx1(om,k,ZZ,0)
        self.ZxfmCx1(om,k,ZZ,2)
        return ZZ
        
    #---------------------------------------------------------------------------
    #  strip plate cross line wave impedance helper
    #---------------------------------------------------------------------------
    def ZxfmCx1(self,om,k,ZZ,i):
        b=self.width()
        if i==0:
            (z1,z2,z3,z4,z5,P,Q)=self.sec.strip_cs(om,k,b)
        elif i==2:
            (z1,z2,z3,z4,z5,P,Q)=self.sec.strip_b(om,k,b)
        else:
            raise ValueError('illegal index %s (0:cs, 2:b)' % (i))
        #~ print k,z1,z2,z3,z4,z5,P,Q,b
        z6=z5/z3 #z0 is always 1
        z7=z1*z4/z2
        z01=z1 #z0 is always 1
        z03=z3 #z0 is always 1
        z12=z1*z2
        z23=z2*z3
        Z=z3*z4-z2*z5
        P2=P*P
        Q2=Q*Q
        PA=1+P2
        PB=1-P2
        QA=1+Q2
        QB=1-Q2
        za=Z*(PB*QA*z12-PA*QB*z03)
        zb=(PA*QA-4*P*Q)*z01*(z3*z4+z2*z5)-PB*QB*(z1*z1*z2*z4+z3*z5)#z0 is always 1
        zc=2*Z*(P*QB*z03-PB*Q*z12)
        zd=-2*Z*(P-Q)*(P*Q-1)*z01
        h=Z*(z01/z23)
        ze=h*(PA*QB*z12-PB*QA*z03)
        zf=2*h*(PB*Q*z03-P*QB*z12)
        h=1/(2*(PA*QA-4*P*Q)*z03*z12-PB*QB*(z12*z12+z03*z03)) # 1/PQdetCx
        za=(za*h).real
        zb=(zb*h).real
        zc=(zc*h).real
        zd=(zd*h).real
        ze=(ze*h).real
        zf=(zf*h).real
        ZZ[i  :i+2,i  :i+2]=array((( za, zb),( zb, ze)))
        ZZ[i  :i+2,i+4:i+6]=array((( zc, zd),(-zd, zf)))
        ZZ[i+4:i+6,i  :i+2]=array((( zc,-zd),( zd, zf)))
        ZZ[i+4:i+6,i+4:i+6]=array((( za,-zb),(-zb, ze)))

    #---------------------------------------------------------------------------
    #  test function: cross line wave impedance compared to product Zxfm * Cx^-1
    #---------------------------------------------------------------------------
    def test(self,om,k):
        ZZ=self.Zstrip(om,k)
        Zfm=zeros((8,8,len(k)),complex)
        Cx=zeros((8,8,len(k)),complex)
        b=self.width()
        for i in (0,2):
            z0=ones((len(k)),Float)
            if i==0:
                (z1,z2,z3,z4,z5,P,Q)=self.sec.strip_cs(om,k,b)
            elif i==2:
                (z1,z2,z3,z4,z5,P,Q)=self.sec.strip_b(om,k,b)
            else:
                raise ValueError('illegal index %s (0:cs, 2:b)' % (i))
            z6=z5/z3 #z0 is always 1
            z7=z1*z4/z2
            Cx[i  :i+2,i  :i+2]=array((( z0, z0),( z2,-z2)))
            Cx[i  :i+2,i+4:i+6]=array((( z1,-z1),( z3, z3)))
            Cx[i+4:i+6,i  :i+2]=array((( z0*P, z0/P),( z2*P,-z2/P)))
            Cx[i+4:i+6,i+4:i+6]=array((( z1*Q,-z1/Q),( z3*Q, z3/Q)))
            Zfm[i  :i+2,i  :i+2]=array((( z4,-z4),( z6, z6)))
            Zfm[i  :i+2,i+4:i+6]=array((( z5, z5),( z7,-z7)))
            Zfm[i+4:i+6,i  :i+2]=array(((-z4*P, z4/P),(-z6*P,-z6/P)))
            Zfm[i+4:i+6,i+4:i+6]=array(((-z5*Q,-z5/Q),(-z7*Q, z7/Q)))
        h1=ZZ[:,:,0]
        h2=dot(Zfm[:,:,0],inverse(Cx[:,:,0]))-h1
        for i in range(64):
            if h1.flat[i]!=0.0:
                h2.flat[i]/=h1.flat[i]
        print(array2string(h2,precision=1,max_line_width=1000,
                           suppress_small=1)) # must be an all zero matrix

#-------------------------------------------------------------------------------
#  beam connection member
#-------------------------------------------------------------------------------
class ConnBeam(ConnMember):

#    __traits__={
#        'sec': Trait(BeamSec()), #section of beam
#        'node': Trait(CNode()), #attach node
#        'y0': 0.0, # y-offset
#        'z0': 0.0, # z-offset
#        'theta': 1.0 # rotation angle in degree
#    }
    sec = Trait(BeamSec())
    node = Trait(CNode())
    y0 = 0.0
    z0 = 0.0
    theta = 1.0
       
    #---------------------------------------------------------------------------
    #  returns a tuple of nodes of length 1
    #---------------------------------------------------------------------------
    def nodes(self):
        return [self.node]

    #---------------------------------------------------------------------------
    #  K matrix co-ordinate transformation, must be called if 
    #  co-ordinate attributes have changed
    #---------------------------------------------------------------------------
    def K_calc(self):
        c=cos(self.theta*pi/180)
        s=sin(self.theta*pi/180)
        y=self.y0+self.node.y
        z=self.z0+self.node.z
        return array(((1.,0,0,0),(0,c,-s,0),(0,s,c,0),(0,y*s-z*c,y*c+z*s,1.)))

#-------------------------------------------------------------------------------
#  point line connection member
#-------------------------------------------------------------------------------
class PointLine(ConnMember):

#    __traits__={
#        'node1': Trait(CNode()), #attach node 1
#        'node2': Trait(CNode()), #attach node 2
#    }   
    node1 = Trait(CNode())
    node2 = Trait(CNode())
    
    #---------------------------------------------------------------------------
    #  K matrices co-ordinate transformation, must be called if 
    #  co-ordinate attributes have changed
    #---------------------------------------------------------------------------
    def K_calc(self):#,x1,x2,x3):
        y1=self.node1.y
        z1=self.node1.z
        y2=self.node2.y
        z2=self.node2.z
        b=((y1-y2)**2+(z1-z2)**2)**0.5
        if b<1e-15:
            c=1
            s=0
        else:
            c=(y2-y1)/b
            s=(z2-z1)/b
        return array((((1.,0,0,0),(0,c,-s,0),(0,s,c,0),(0,y1*s-z1*c,y1*c+z1*s,1.)),
        #self.Kb=array
                    ((1.,0,0,0),(0,c,-s,0),(0,s,c,0),(0,y2*s-z2*c,y2*c+z2*s,1.))))

#-------------------------------------------------------------------------------
#  base class for connections
#-------------------------------------------------------------------------------
class Connect(TrObj):
    
    member_classes=[] # type of possible members, none here in base class

    def __init__(self,*kw,**args):
        TrObj.__init__(self,*kw,**args)
        self.nodes=[] # list of all nodes in the connection
        self.members=[] # list of all members in the connection
        self.results={} # dictionary of results, key is a (freq,from,to) tuple, value
                    # is a n(to) x n(from) matrix
    
    #---------------------------------------------------------------------------
    #  add a node to __nodes if not already there and returns its number
    #---------------------------------------------------------------------------
    def add_node(self,node):
        if not node in self.nodes:
            self.nodes.append(node)
        return self.nodes.index(node)
        
    #---------------------------------------------------------------------------
    #  add a member if of legal type
    #---------------------------------------------------------------------------
    def add_member(self,member):
        if member.__class__ in self.member_classes:
            self.members.append(member)
            for node in member.nodes():
                self.add_node(node)
        else:
            raise TypeError('member type illegal for this connection %s'
                            % (member))
            
    #---------------------------------------------------------------------------
    #  returns a matrix of transmission loss factors for transmission from
    #  'frm' to 'tom' at frequency 'freq' 
    #---------------------------------------------------------------------------
    def get_trans(self,freq,frm,wt,tom):
        for m in [frm,tom]:
            if not m in self.members:
                raise ValueError('member %s not in connection' % (m))
            if not m.active():
                raise TypeError('member %s not of active type' % (m))
        if not (freq,frm,wt,tom) in self.results:
            self.calculate(freq)
        return self.results[(freq,frm,wt,tom)]
        
    #---------------------------------------------------------------------------
    #  dummy calculation function
    #---------------------------------------------------------------------------
    def calculate(self,freq):
        print("calc")
        for x1 in itertype(self.members,InfPlate):
            for wt in range(3):
                for x2 in itertype(self.members,InfPlate):
                    self.results[(freq,x1,wt,x2)]=1.0

#-------------------------------------------------------------------------------
#  base class for line connections
#-------------------------------------------------------------------------------
class LineConnect(Connect):

    #---------------------------------------------------------------------------
    #  calculation function
    #---------------------------------------------------------------------------
    def calculate(self,freq,points=50):
        om=2*pi*freq
        #
        # collect all free wavenumbers of the connection
        #
        kall=[0.]
        for x1 in itertype(self.members,InfPlate):
            for k1 in (  x1.sec.kc2(om)**0.5,
                        x1.sec.ks2(om)**0.5,
                        x1.sec.kb4(om)**0.25):
                if not k1 in kall:
                    kall.append(k1)
        #~ for x1 in itertype(self.members,StripPlate):
            #~ for k1 in (  x1.sec.kc2(om)**0.5,
                        #~ x1.sec.ks2(om)**0.5,
                        #~ x1.sec.kb4(om)**0.25):
                #~ if not k1 in kall:
                    #~ kall.append(k1)
        kall.sort()
        #
        # loop over all active systems and wavetypes
        #
        for x1 in itertype(self.members,InfPlate):
            for wt in range(3):
                #
                # determine upper integration limit
                #
                if wt==0:
                    kmax=x1.sec.kc2(om)**0.5
                elif wt==1:
                    kmax=x1.sec.ks2(om)**0.5
                elif wt==2:
                    kmax=x1.sec.kb4(om)**0.25
                #
                # integration over k
                #
                (k,weights)=BreakGauss(compress(less_equal(kall,kmax),kall),points)
                tau=self.calck(om,k,x1,wt)
                for x2 in itertype(self.members,InfPlate):
                    self.results[(freq,x1,wt,x2)]=dot(tau[x2],weights)/kmax

#-------------------------------------------------------------------------------
#  continuous line connection
#-------------------------------------------------------------------------------
class ContLineConnect(LineConnect):
    member_classes=[InfPlate,StripPlate,ConnBeam] # type of possible members

    #---------------------------------------------------------------------------
    #  calculation per k, frm, sender wave type "wt"
    #  om: omega
    #  k: N-element vector of wave numbers
    #  frm: sender member
    #  wt: sender wave type
    #  returns dictionary of N-by-3 matrices
    #---------------------------------------------------------------------------
    def calck(self,om,k,frm,wt):
        #
        # assembly of system matrix M
        #
        lk=len(k)
        dimn=len(self.nodes)*4
        M=zeros((dimn,dimn,lk),complex)
        # loop over all InfPlate members
        for x1 in itertype(self.members,InfPlate):
            i=self.nodes.index(x1.node)*4
            x1.sec.calc_inf(om,k)
            Z=x1.sec.Zinf
            K=x1.K
            KT=transpose(K)
            for n in range(lk):
                M[i:i+4,i:i+4,n]+=dot(dot(K,Z[:,:,n]),KT)
        # loop over all ConnBeam members
        for x1 in itertype(self.members,ConnBeam):
            i=self.nodes.index(x1.node)*4
            Z=x1.sec.Z(om,k)
            K=x1.K
            KT=transpose(K)
            for n in range(lk):
                M[i:i+4,i:i+4,n]+=dot(dot(K,Z[:,:,n]),KT)
        # loop over all StripPlate members
        for x1 in itertype(self.members,StripPlate):
            i1=self.nodes.index(x1.node1)*4
            i2=self.nodes.index(x1.node2)*4
            Z=x1.Zstrip(om,k)#.real
            K=x1.K
            KT=transpose(K)
            for n in range(lk):
                Zi=dot(dot(K,Z[:,:,n]),KT)
                M[i1:i1+4,i1:i1+4,n]+=Zi[0:4,0:4]
                M[i1:i1+4,i2:i2+4,n]+=Zi[0:4,4:8]
                M[i2:i2+4,i1:i1+4,n]+=Zi[4:8,0:4]
                M[i2:i2+4,i2:i2+4,n]+=Zi[4:8,4:8]
        #
        # assembly of rhs N
        #
        N=zeros((dimn,lk),complex)
        i=self.nodes.index(frm.node)*4
        for n in range(lk):
            N[i:i+4,n]=2*dot(frm.K,diagonal(frm.sec.Zinf[:,:,n])*frm.sec.Ca[:,wt,n])
        #
        # calculation of junction velocity
        #
        vg=zeros((dimn,lk),complex)
        for n in range(lk):
            vg[:,n] = solve(M[:,:,n],N[:,n])
        #~ return self.calc_tau(lk,vg,om,k,frm,wt)
        
    #~ def calc_tau(self,lk,vg,om,k,frm,wt):      
        #
        # calculation of input power
        #
        pin=frm.sec.P[wt]
        #
        # calculation of output power
        #
#        test=zeros(lk,Float)
        v=zeros((4,lk),complex)
        vca=zeros((2,lk),complex)
        vsa=zeros((2,lk),complex)
        pout=zeros((3,lk),Float)
        # loop over out members
        tau={}
        for x2 in itertype(self.members,InfPlate):
            i=self.nodes.index(x2.node)*4
            sec=x2.sec
            Z=sec.Zinf/complex(0,om)
            Vc=sec.Vca
            KT=transpose(x2.K)
            # transformation of velocity in local coordinates
            for n in range(lk):
                v[:,n]=dot(KT,vg[i:i+4,n])
            if frm==x2:
                v-=sec.Ca[:,wt]
            # calculation of output in-plane velocities
            for n in range(lk):
                vca[:,n]=dot(Vc[:,:,n],v[0:2,n])
            vsa=v[0:2]-vca
            # loop over all k
            for n in range(lk):
                pout[:,n]=array(
                    (dot(conjugate(vca[:,n]),dot(Z[0:2,0:2,n],vca[:,n])),
                     dot(conjugate(vsa[:,n]),dot(Z[0:2,0:2,n],vsa[:,n])),
                     dot(conjugate(v[2:4,n]),dot(Z[2:4,2:4,n],v[2:4,n])))).real
            tau[x2]=0.5*(pout/pin)
#            print tau
#            print "from",frm,"to",x2,"wt",wt,array2string(tau,precision=3)
#            test+=sum(tau[x2])
#        print "test",test
        return tau

#-------------------------------------------------------------------------------
#  point line connection
#-------------------------------------------------------------------------------
class PointLineConnect(LineConnect):

#    __traits__={
#        'NFT': Trait(2, TraitRange(2,maxint)), #number of sampling points for Fourier Transform
#        'dist': 1.0 #distance between point connections 
#    }   
    NFT = Range(2,2**24,2)
    dist = 1.0
    member_classes=[InfPlate,StripPlate,PointLine] # type of possible members

    #---------------------------------------------------------------------------
    #  calculation per k, frm, sender wave type "wt"
    #  om: omega
    #  k: N-element vector of wave numbers
    #  frm: sender member
    #  wt: sender wave type
    #  returns dictionary of N-by-3 matrices
    #---------------------------------------------------------------------------
    def calck(self,om,k,frm,wt):
        #
        # assembly of system matrix M
        #
        lk=len(k)
        # how many nodes w point lines are in the connection ?
        plmembers=[]
        plnodes=[]
        for x1 in itertype(self.members,PointLine):
            plmembers.append(x1)
            plnodes.append(x1.node1)
            plnodes.append(x1.node2)
        npl=len(plnodes)
        # number of nodes
        nn=len(self.nodes)
        # number of unknowns
        dimn=nn*4*self.NFT+6*npl 
        # matrix
        oneslk=ones(lk,complex)
        L=zeros((4,6,lk),complex)
        L[1,1]=oneslk
        L[2,2]=oneslk
        L[3,3]=oneslk
        identity6=identity(6,complex)
        M=zeros((dimn,dimn,lk),complex)
        kk=k
        # n loop backward 
        for l in range(self.NFT-1,-1,-1):
            k=kk+l/(2*pi*self.dist)
            L[0,0]=-1j*k
            L[2,4]=-1j*k
            L[1,5]=1j*k
            # loop over all InfPlate members
            for x1 in itertype(self.members,InfPlate):
                i=self.nodes.index(x1.node)*4+4*l*nn
                x1.sec.calc_inf(om,k)
                Z=x1.sec.Zinf
                K=x1.K
                KT=transpose(K)
                for n in range(lk):
                    M[i:i+4,i:i+4,n]+=dot(dot(K,Z[:,:,n]),KT)
            # loop over all StripPlate members
            for x1 in itertype(self.members,StripPlate):
                i1=self.nodes.index(x1.node1)*4+4*l*nn
                i2=self.nodes.index(x1.node2)*4+4*l*nn
                Z=x1.Zstrip(om,k)#.real
                K=x1.K
                KT=transpose(K)
                for n in range(lk):
                    Zi=dot(dot(K,Z[:,:,n]),KT)
                    M[i1:i1+4,i1:i1+4,n]+=Zi[0:4,0:4]
                    M[i1:i1+4,i2:i2+4,n]+=Zi[0:4,4:8]
                    M[i2:i2+4,i1:i1+4,n]+=Zi[4:8,0:4]
                    M[i2:i2+4,i2:i2+4,n]+=Zi[4:8,4:8]
            # loop over all PointLine members
            for x1 in itertype(self.members,PointLine):
                i1=self.nodes.index(x1.node1)*4+4*l*nn
                i2=self.nodes.index(x1.node2)*4+4*l*nn
                i3=plnodes.index(x1.node1)*6+4*self.NFT*nn
                i4=plnodes.index(x1.node2)*6+4*self.NFT*nn
                i5=plmembers.index(x1)*12+4*self.NFT*nn
                for n in range(lk):
                    LL=L[:,:,n]
                    LH=transpose(conjugate(LL))
                    M[i1:i1+4,i3:i3+6,n]=1j*om*dot(x1.K[0],LL)/self.dist
                    M[i2:i2+4,i4:i4+6,n]=1j*om*dot(x1.K[1],LL)/self.dist
                    M[i5:i5+6,i1:i1+4,n]=-dot(LH,transpose(x1.K[0]))
                    M[i5:i5+6,i2:i2+4,n]=-dot(LH,transpose(x1.K[1]))
                    M[i5+6:i5+12,i3:i3+6,n]=-identity6               
                    M[i5+6:i5+12,i4:i4+6,n]=identity6                 
        #
        # assembly of rhs N
        #
        N=zeros((dimn,lk),complex)
        i=self.nodes.index(frm.node)*4
        for n in range(lk):
            N[i:i+4,n]=2*dot(frm.K,diagonal(frm.sec.Zinf[:,:,n])*frm.sec.Ca[:,wt,n])
        #
        # calculation of junction velocity
        #
        vg=zeros((dimn,lk),complex)
        for n in range(lk):
            vg[:,n]=solve(M[:,:,n],N[:,n])
#                print vg[:,0]
        #
        # calculation of input power
        #
        pin=frm.sec.P[wt]
        #
        # calculation of output power
        #
#        test=zeros(lk,Float)
        v=zeros((4,lk),complex)
        vca=zeros((2,lk),complex)
        vsa=zeros((2,lk),complex)
        pout=zeros((3,lk),Float)
        # loop over out members
        tau={}
        for l in range(self.NFT):
            k=kk+l/(2*pi*self.dist)
            for x2 in itertype(self.members,InfPlate):
                i=self.nodes.index(x2.node)*4+4*l*nn
                sec=x2.sec
                sec.calc_inf(om,k)
                Z=sec.Zinf/complex(0,om)
                Vc=sec.Vca
                KT=transpose(x2.K)
                # transformation of velocity in local coordinates
                for n in range(lk):
                    v[:,n]=dot(KT,vg[i:i+4,n])
    #                    print v[:,0]
                if frm==x2 and l==0:
                    v-=sec.Ca[:,wt]
                # calculation of output in-plane velocities
                for n in range(lk):
                    vca[:,n]=dot(Vc[:,:,n],v[0:2,n])
                vsa=v[0:2]-vca
                # loop over all k
                for n in range(lk):
                    pout[:,n]=abs(array(
                        (dot(conjugate(vca[:,n]),dot(Z[0:2,0:2,n],vca[:,n])),
                         dot(conjugate(vsa[:,n]),dot(Z[0:2,0:2,n],vsa[:,n])),
                         dot(conjugate(v[2:4,n]),dot(Z[2:4,2:4,n],v[2:4,n])))).real)
                if l==0:
                    tau[x2]=0.5*(pout/pin)
                else:
                    tau[x2]+=0.5*(pout/pin)
 #               if wt==2:
#                    print choose(tau[x2][2]>1,(0,k))
    #            print "from",frm,"to",x2,"wt",wt,array2string(tau,precision=3)
    #            test+=sum(tau[x2])
    #        print "test",test
        return tau
