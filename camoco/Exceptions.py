# Exception abstract class
class CamocoError(Exception):
    pass

class CamocoExistsError(CamocoError):
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

class CamocoGeneNameError(CamocoError):
    '''
        Gene names must be beautiful snowflakes.
    '''
    def __init__(self,expr,message='',*args):
        self.expr = expr
        self.message = 'Gene names must be unique:' + message.format(args)

class CamocoGeneAbsentError(CamocoError):
    '''
        Gene is missing from dataset
    '''
    def __init__(self,expr,message='',*args):
        self.expr = expr
        self.message = (
            'This gene is not in the dataset:' + message.format(args)
        )

class CamocoAccessionNameError(CamocoError):
    '''
        Accession names must be Unique.
    '''
    def __init__(self,expr,message='',*args):
        self.expr = expr
        self.message = (
            'Accession names must be unique:' + message.format(args)
        )

class CamocoZeroWindowError(CamocoError):
    def __init__(self,expr,message,*args):
        self.expr = expr
        self.message = (
            'Operation requiring window, but window is 0:' + \
            message.format(args)
        )

class CamocoInteractive(CamocoError):
    def __init__(self,expr=None,message='',*args):
        self.expr = expr
        self.message = 'Camoco interactive ipython session.'
