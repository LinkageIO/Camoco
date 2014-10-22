from flask import Flask, url_for, jsonify, request
app = Flask(__name__)

import camoco as co
import sys
import json
from math import isinf

networks = {x:co.COB(x) for x in ['ROOT','SAM','PAN']}
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
        
@app.route("/Ontology/Terms/<path:term_name>")
def Ontology_Terms(term_name):
    return jsonify({
        "data": [ (term.id,term.name,
                    len(term.snp_list),
                    len(term.gene_list),
                    len(set(term.flanking_genes(networks['ROOT'].refgen)).union([gene for gene in term.gene_list if gene in networks['ROOT'].refgen]))
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
        genes = term.gene_list
    gene_annots = co.Annotation('Func')[genes]
    for net in networks.values():
        gene_annots.insert(5,'In {}'.format(net.name), ['true' if gene in net else 'false' for gene in genes])
    gene_annots.insert(5,'Term SNPs',["\n".join([snp.summary() for snp in sorted(term.flanking_snps(gene))]) for gene in genes])
    gene_annots.insert(5,'Global Degree',[str(cob.global_degree(gene)) for gene in genes])
    return jsonify({
        'data' : list(gene_annots.itertuples(index=False))
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
                    'distance' : fix_inf(distance)
                }
            } for source,target,score,distance in subnet.itertuples(index=False)
        ]
        # 
        return jsonify(net)
    except Exception as e:
        return 500 

def fix_inf(val):
    if isinf(val):
        return -1
    else:
        return val

if __name__ == "__main__":
    app.run('0.0.0.0',8080)
