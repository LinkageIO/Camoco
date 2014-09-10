from flask import Flask, url_for, jsonify
app = Flask(__name__)

import camoco as co
import sys
import json
from math import isinf

@app.route('/')
def index():
    return app.send_static_file("index.html")

@app.route("/Camoco/available_datasets")
def all_available_datasets():
    return str(co.available_datasets())
 
@app.route("/Camoco/available_datasets/<path:type>")
def available_datasets(type=None):
    return jsonify({ "data" : list(co.available_datasets(type)[['Name','Description']].itertuples(index=False))})

        
@app.route("/Ontology/Terms/<path:term_name>")
def Ontology_Terms(term_name):
    return jsonify({
        "data": [ (term.id,term.name,len(term.gene_list)) for term in co.Ontology(term_name).iter_terms()]
    })

@app.route("/api/COB/<network_name>/<ontology>/<term>")
def COB_network(network_name,ontology,term):
    try:
        net = {}
        subnet = co.COB(network_name).subnetwork(co.Ontology(ontology)[term].gene_list)
        net['nodes'] = [ {'data':{'id':str(x)}} for x in set(subnet.source).union(subnet.target) ]
        net['edges'] = [
            {
                'data':{
                    'source': source,
                    'target' : target,
                    'score' : score,
                    'distance' : fix_inf(distance)
                }
            } for source,target,score,distance in subnet.itertuples(index=False)
        ]
        return jsonify(net)
        
        return(str(network_name)+str(ontology_term))
    except Exception as e:
        return 500 

@app.errorhandler(500)
def server_error(error):
    return "Nopeers:{}".format(error),500

def fix_inf(val):
    if isinf(val):
        return -1
    else:
        return val

if __name__ == "__main__":
    app.debug=True
    app.run('0.0.0.0')
