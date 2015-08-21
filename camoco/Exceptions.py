# Exception abstract class
class CamocoError(Exception):
    pass

class CamocoExistsError(CamocoError):
    '''
        You tried to create a camoco object which already exists
        under the same name,type combination.
    '''
    def __init__(self,expr,message=''):
        self.expr = expr
        self.message = (
            'You are trying to create a Camoco based object'
            'That already exists' + message
        )

class CamocoGeneNameError(CamocoError):
    '''
        Gene names must be beautiful snowflakes.
    '''
    def __init__(self,expr,message=''):
        self.expr = expr
        self.message = 'Gene names must be unique:' + message

class CamocoAccessionNameError(CamocoError):
    '''
        Accession names must be Unique.
    '''
    def __init__(self,expr,message=''):
        self.expr = expr
        self.message = (
            'Accession names must be unique:' + message
        )
