import camoco as co

class Analysis(object):
    '''
        Perform an analysis based on CLI arguments:
        set up, event loop, tear down
    '''
    def __init__(self):
        # Init needs to just store args and other analysis level data
        self.args = args
        self.tag = "Analysis"

    def __call__(self):
        set_up()
        event_loop()
        tear_down()

    def set_up(self):
        pass

    def event_loop(self):
        pass

    def tear_down(self):
        pass             

    # ------------------------------------------------------------------------
    # Extra methods should fit into the above methods

    def _generate_output_file(self, filename):
        if args.out != sys.stdout:
            args.out = "{}_Locality.csv".format(args.out.replace('.csv',''))        
        if os.path.dirname(args.out) != '':
            os.makedirs(os.path.dirname(args.out),exist_ok=True)
        if os.path.exists("{}_Locality.csv".format(args.out.replace('.csv',''))):
            print(
                "{}_Locality.csv exists! Skipping!".format(
                    args.out.replace('.csv','')
                )
            )
            return None

    def _build_camoco_objects(self):
        pass
