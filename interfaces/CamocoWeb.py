from flask import Flask, url_for, jsonify, request, send_from_directory
app = Flask(__name__)

import camoco as co
import sys
import json
from math import isinf
import numpy as np
 
networks = {x:co.COB(x) for x in ['ZmRoot']}
ZM = co.RefGen('Zm5bFGS')

@app.route('/')
def index():
    return app.send_static_file("index.html")

@app.route('/static/<path:path>')
def send_js(path):
    return send_from_directory('static',path)

@app.route("/available_datasets")
def all_available_datasets():
    return str(co.available_datasets())
 
@app.route("/available_datasets/<path:type>")
def available_datasets(type=None,*args):
    return jsonify({ 
        "data" : list(
            co.available_datasets(type)[
                ['Name','Description']
            ].itertuples(index=False))
        }
    )

@app.route("/Expression/<network_name>/<ontology_name>/<term_name>",
        methods=['GET','POST'])
def expression(network_name,ontology_name,term_name):
    pass
    expr = networks[network_name].heatmap(
        network[network_name].gene_expression_vals(
            co.Ontology(ontology_name)[term_name].gene_list,
            zscore=True
        ),
        cluster_x=True,
        cluster_y=False,
        png_encode=True
    )
    return """<html><body>
        <img src="data:image/png;base64,{}"/>
        </body></html>""".format(expr)

    return jsonify({
        'rownames' : list(expr.index.values),
        'colnames' : list(expr.columns.values),
        'data' :
            [{'x': rowid,
              'y': colid,
              'value': fix_val(value)} for rowid,row in \
               enumerate(expr.itertuples(index=False)) for colid,value in \
               enumerate(row) 
            ]
    })
        
@app.route("/Ontology/Terms/<path:term_name>")
def Ontology_Terms(term_name):
    return jsonify({
        'data': [ 
            (term.name,
             term.desc,
             len(term.locus_list),
             len(ZM.candidate_genes(term.effective_snps(window_size=50000)))
                ) 
            for term in co.Ontology(term_name).iter_terms()]
    })

@app.route('/Annotations/<network_name>/<ontology_name>/<term_name>',
        methods=['GET','POST'])
def Annotations(network_name,ontology_name,term_name):
    # Retrieve SNPs from 
    cob = networks[network_name]
    term = co.Ontology(ontology_name)[term_name]
    genes = cob.refgen.from_ids(request.args.get('genes').split(','))
    try:
        gene_annots = co.Annotation('ZMFunc')[genes]
    except ValueError as e:
        return jsonify({})
    for net in networks.values():
        gene_annots.insert(
            5,'In {}'.format(net.name),
            ['true' if gene in net else 'false' for gene in genes]
        )
    gene_annots.insert(
        5,'Term SNPs',
        ["\n".join([snp.summary() for snp in \
            sorted(term.flanking_snps(gene))]) for gene in genes]
    )
    gene_annots.insert(
        5,'Global Degree',
        [str(cob.global_degree(gene)) for gene in genes]
    )
    return jsonify({
        'data' : list(gene_annots.itertuples(index=True)),
        'header' : ['Gene'] + list(gene_annots.columns.values)
    })

@app.route("/COB/<network_name>/<ontology>/<term>")
def COB_network(network_name,ontology,term):
    net = {}
    cob = networks[network_name]
    term = co.Ontology(ontology)[term]
    nodes = []
    seen = set()
    effective_snps = term.effective_snps(window_size=100000)
    candidate_genes = cob.refgen.candidate_genes(
        effective_snps,gene_limit=10,chain=False
    )
    locality = cob.locality(
        set([item for sublist in candidate_genes \
            for subsublist in sublist \
            for item in subsublist]
        )
    )
    for snp,genes in zip(effective_snps,candidate_genes):
        # Put the SNP in there 
        #nodes.append({'data':{'id':str(snp.summary())}})
        # append 
        for gene in [item for sublist in genes for item in sublist]:
            if gene.id not in seen:
                try:
                    local_degree = locality.ix[gene.id]['local']
                    global_degree = locality.ix[gene.id]['global']
                except KeyError as e:
                    local_degree = 0
                gene_locality = local_degree / global_degree
                if np.isnan(gene_locality):
                    gene_locality = 0
                nodes.append(
                        {'data':{
                            'id':str(gene.id), 
                            'locality':int(gene_locality*100),
                            'degree':int(local_degree), 
                            'snp':str(snp.summary()),
                            'gdegree':int(global_degree) 
                        }}  
                )
            seen.add(gene.id)
    # Now do the edges
    subnet = cob.subnetwork(
        cob.refgen.candidate_genes(
            term.effective_snps(window_size=100000),
            gene_limit=10
        )
    )
    subnet.reset_index(inplace=True)
    net['nodes'] = nodes
    net['edges'] = [
        {
            'data':{
                'source': source,
                'target' : target,
                'score' : score,
                'distance' : fix_val(distance)
            }
        } for source,target,score,significant,distance in subnet.itertuples(index=False)
    ]
    return jsonify(net)

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
    app.run('0.0.0.0',50000)
