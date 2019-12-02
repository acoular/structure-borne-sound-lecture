#-------------------------------------------------------------------------------
# itertype.py
# 
# (c) 12/2003 Ennes Sarradj
#
# Implementation of an iterator iterating over that objects in a given list
# that are instances of a given class or a subclass thereof
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#  iterator over objects in list "list" that are instances of class "type"
#-------------------------------------------------------------------------------
# =============================================================================
# class itertype(object):
#     
#     def __init__(self,list,type):
#         self.list=list
#         self.type=type
#         self.index=0
#     
#     # implementation according to the Python iterator protocol of 2.3
#     def __iter__(self):
#         return self
#     
#     # return next available object or StopIteration exception
#     def next(self):
#         try:
#             while 1:
#                 x=self.list[self.index]
#                 self.index+=1
#                 if isinstance(x,self.type):
#                     break
#             return x
#         except IndexError:
#             raise StopIteration
#  
# =============================================================================
def itertype(list_,type_):
    return filter(lambda x:isinstance(x,type_),list_)