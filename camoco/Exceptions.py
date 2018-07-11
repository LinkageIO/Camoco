# Exception abstract class
class CamocoError(Exception): # pragma: no cover
    pass

class CamocoExistsError(CamocoError): # pragma: no cover
    '''
        You tried to create a camoco object which already exists
        under the same name,type combination.
    '''
    def __init__(self,expr,message='',*args):
        self.expr = expr
        self.message = (
            'You are trying to create a Camoco based object'
            'That already exists' + message.format(*args)
        )

class CamocoGeneNameError(CamocoError): # pragma: no cover
    '''
        Gene names must be beautiful snowflakes.
    '''
    def __init__(self,expr,message='',*args):
        self.expr = expr
        self.message = 'Gene names must be unique:' + message.format(args)

class CamocoGeneAbsentError(CamocoError): # pragma: no cover
    '''
        Gene is missing from dataset
    '''
    def __init__(self,expr,message='',*args):
        self.expr = expr
        self.message = (
            'This gene is not in the dataset:' + message.format(args)
        )

class CamocoAccessionNameError(CamocoError): # pragma: no cover
    '''
        Accession names must be Unique.
    '''
    def __init__(self,expr,message='',*args):
        self.expr = expr
        self.message = (
            'Accession names must be unique:' + message.format(args)
        )

class CamocoZeroWindowError(CamocoError): # pragma: no cover
    def __init__(self,expr,message,*args):
        self.expr = expr
        self.message = (
            'Operation requiring window, but window is 0:'
        )

class CamocoInteractive(CamocoError): # pragma: no cover
    def __init__(self,expr=None,message='',*args):
        self.expr = expr
        self.message = 'Camoco interactive ipython session.'
