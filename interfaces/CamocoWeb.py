from flask import Flask, url_for, jsonify, request
app = Flask(__name__)

import camoco as co
import sys
import json
from math import isinf
import numpy as np

networks = {x:co.COB(x) for x in ['ZmRoot','ZmSAM','ZmPAN']}
ZM = co.RefGen('Zm5bFGS')

@app.route('/')
def index():
    return app.send_static_file("index.html")

@app.route("/Camoco/available_datasets")
def all_available_datasets():
    return str(co.available_datasets())
 
@app.route("/Camoco/available_datasets/<path:type>")
def available_datasets(type=None,*args):
    return jsonify({ "data" : list(co.available_datasets(type)[['Name','Description']].itertuples(index=False))})

@app.route("/api/Expression/<network_name>/<ontology_name>/<term_name>",methods=['GET','POST'])
def expression(network_name,ontology_name,term_name):
    expr = co.COB(network_name).heatmap(co.COB(network_name).gene_expression_vals(co.Ontology(ontology_name)[term_name].gene_list,zscore=True),cluster_x=True,cluster_y=False,png_encode=True)
    return """<html><body>
        <img src="data:image/png;base64,{}"/>
        </body></html>""".format(expr)

    return jsonify({
        'rownames' : list(expr.index.values),
        'colnames' : list(expr.columns.values),
        'data' :
            [{'x': rowid,
              'y': colid,
              'value': fix_val(value)} for rowid,row in enumerate(expr.itertuples(index=False)) for colid,value in enumerate(row) ]
    })
        
@app.route("/Ontology/Terms/<path:term_name>")
def Ontology_Terms(term_name):
    return jsonify({
        'data': [ (term.id,term.name,
                    len(term.snp_list),
                    len(term.gene_list),
                    #len([gene for gene in term.gene_list if gene in networks['ROOT'].refgen])
                ) for term in co.Ontology(term_name).iter_terms()]
    })

@app.route('/api/Annotations/<network_name>/<ontology_name>/<term_name>',methods=['GET','POST'])
def Annotations(network_name,ontology_name,term_name):
    # Retrieve SNPs from 
    cob = co.COB(network_name)
    term = co.Ontology(ontology_name)[term_name]
    if len(term.gene_list) == 0:
        genes = term.flanking_genes(cob.refgen,gene_limit=4,window_size=100000)
    else:
        genes = term.gene_list
    try:
        gene_annots = co.Annotation('ZMFunc')[genes]
    except ValueError as e:
        return jsonify({})
    for net in networks.values():
        gene_annots.insert(5,'In {}'.format(net.name), ['true' if gene in net else 'false' for gene in genes])
    gene_annots.insert(5,'Term SNPs',["\n".join([snp.summary() for snp in sorted(term.flanking_snps(gene))]) for gene in genes])
    gene_annots.insert(5,'Global Degree',[str(cob.global_degree(gene)) for gene in genes])
    return jsonify({
        'data' : list(gene_annots.itertuples(index=True)),
        'header' : ['Gene'] + list(gene_annots.columns.values)
    })

@app.route("/api/COB/<network_name>/<ontology>/<term>")
def COB_network(network_name,ontology,term):
    try:
        net = {}
        cob = co.COB(network_name)
        term = co.Ontology(ontology)[term]
        if len(term.gene_list) == 0:
            genes = term.flanking_genes(cob.refgen,gene_limit=4,window_size=100000)
        else:
            genes = term.gene_list
        subnet = cob.subnetwork(genes,min_distance=100000)
        net['nodes'] = [ {'data':{'id':str(x.id), 'index':i, 'gdegree':cob.global_degree(x) }} for i,x in enumerate(genes) ]
        net['edges'] = [
            {
                'data':{
                    'source': source,
                    'target' : target,
                    'score' : score,
                    'distance' : fix_val(distance)
                }
            } for source,target,score,distance in subnet.itertuples(index=False)
        ]
        return jsonify(net)
    except Exception as e:
        return 500 

def fix_val(val):
    if isinf(val):
        return -1
    if np.isnan(val):
        # because Fuck JSON
        return "null"
    else:
        return val

if __name__ == "__main__":
    app.debug = True
    app.run('0.0.0.0',8080)
