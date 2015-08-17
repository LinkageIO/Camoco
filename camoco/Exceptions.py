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
        self.message = message
