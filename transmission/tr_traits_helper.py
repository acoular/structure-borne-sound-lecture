#-------------------------------------------------------------------------------
# tr_traits_helper.py
# 
# (c) 12/2003 Ennes Sarradj
#
# Implementation of a trait handler for positive floats
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#  Imports:
#-------------------------------------------------------------------------------

from traits.trait_handlers import TraitHandler
import types

#-------------------------------------------------------------------------------
#  handler class for positive floats
#-------------------------------------------------------------------------------
class PosFloat(TraitHandler):
    
    def validate ( self, object, name, value ):
        try:
            cvalue=float(value)
            if ((type( cvalue ) == float) and (cvalue > 0)):
                return cvalue
        except:
            pass
        self.error( object, name, value )

    def info ( self ):
        return 'a positive real number'

