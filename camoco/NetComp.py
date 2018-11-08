from .COB import COB

class NetComp(Camoco):

    def __init__(self,name,networks):
        self.networks = set()

        # Add all the networks
        for n in networks:
            self.add_network(n)

    def add_network(self,net):
        '''
            Add a network (COB) to the 
            NetComp object.
        '''
        if isinstance(net,str):
            net = COB(net)
        if not isinstance(net,COB):
            raise ValueError(f'a valid network must be provided')
        self.networks.add(net)
