#-------------------------------------------------------------------------------
# tr_traits_helper.py
# 
# (c) 12/2003 Ennes Sarradj
#
# Implementation of traitlets helpers used by the transmission package
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#  Imports:
#-------------------------------------------------------------------------------

from traitlets import TraitError, TraitType

#-------------------------------------------------------------------------------
#  trait class for positive floats
#-------------------------------------------------------------------------------
class PosFloat(TraitType):
    
    def validate ( self, object, value ):
        try:
            cvalue=float(value)
            if ((type( cvalue ) == float) and (cvalue > 0)):
                return cvalue
        except:
            pass
        self.error( object, value )

    def info ( self ):
        return 'a positive real number'


class Range(TraitType):

    def __init__(self, low=None, high=None, value=None, **metadata):
        self.low = low
        self.high = high
        self._value_type = float if any(type(v) is float for v in (low, high, value)) else int
        super().__init__(default_value=value, **metadata)

    def validate(self, object, value):
        if self._value_type is int:
            if type(value) is not int:
                self.error(object, value)
            cvalue = value
        else:
            if type(value) is bool or not isinstance(value, (int, float)):
                self.error(object, value)
            cvalue = float(value)
        if self.low is not None and cvalue < self.low:
            self.error(object, value)
        if self.high is not None and cvalue > self.high:
            self.error(object, value)
        return cvalue

    def info(self):
        return f'{self.low} <= a {self._value_type.__name__} <= {self.high}'


class Instance(TraitType):

    def __init__(self, klass, default_value=None, **metadata):
        self.klass = klass
        super().__init__(default_value=default_value, allow_none=True, **metadata)

    def validate(self, object, value):
        if value is None or isinstance(value, self.klass):
            return value
        self.error(object, value)

    def info(self):
        return f'an instance of {self.klass.__name__} or None'


def Trait(default, *handlers):
    if handlers:
        if len(handlers) != 1 or not isinstance(handlers[0], TraitType):
            raise TraitError('unsupported Trait declaration')
        handler = handlers[0]
        handler.default_value = default
        return handler
    if isinstance(default, type):
        return Instance(default)
    return Instance(type(default), default)
