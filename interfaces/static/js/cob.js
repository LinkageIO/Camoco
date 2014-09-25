    function Graph(params){
        defaults = {
            'div':$('<div>'),
        }    
        this.params = $.extend(true,defaults,params)
        this.selected = []

        this.cy = cytoscape(options = {
            container : this.params.div[0],    
            // General Options 
            hideEdgesOnViewport: true,
            hideLabelsOnViewport: true,
            // Style
            style: cytoscape.stylesheet()
            .selector('node')
            .css({
                'background-color': '#144566',
                'border-width': 1,
                'border-color': '#000',
                'height' : 'data(degree)',
                'width' : 'data(degree)',
                'content':'data(id)',
                'text-halign':'right',
                'min-zoomed-font-size':1
            })
            .selector(':selected')
            .css({
                'border-width': 3,
                'border-color': '#09BA00'
            })
            .selector('edge')
            .css({
                'opacity': '0.25',
                'width': 'mapData(score, 3, 7, 1, 20)',
                'curve-style': 'haystack' // fast edges!
            }),
            layout: {
                name: 'concentric',
                concentric: function(){ return this.data('score'); },
                levelWidth: function( nodes ){ return 10; },
                padding: 10
            },
            ready : function(){
                console.log('Cytoscape Web, ready to rock!')
            }
        })
        this.show = function(score){
            this.cy.edges('edge[score < '+score+']').hide()
            this.cy.edges('edge[score >= '+score+']').show()
        }
    }

    function Tab(params){
        this.name = params.name
        this.div = params.div || $('<div>')
    }
    
    function Bar(params){
        // expand with defaults
        params = $.extend(true,{
            'title': 'Hello',
            'body' : 'World'
        },params)
        this.div = params.div || $('<div>');
        this.div.append($('<a>',{class:'collapse'}))
                .append($('<span>',{class:'title'}).html(params.title))
                .append($('<span>',{class:'body'}).html(params.body))
    }

    function Menu(params){
        // Initialize Tables
        this.tabs = []
        this.handle = []
        this.tables = []
        defaults = {
            'div' : $('<div>',{'class':'menu'})
        }
        this.params = $.extend(true,defaults,params)
        // Add the ontology and term table
        // Initialize Actions
        /*
        */

    } // End Menu
    Menu.prototype.add_table = function(name, headers, params){
        this.params.div
            .append(
                $('<table>',{class: name+' display',cellspacing:'0',width:'100%'})
                .append($('<thead>')
                        .append($('<tr>')
                        )
                )
            )
        for (var i=0; i<headers.length; i++){
            this.params.div.find('table.'+name+' thead tr')
            .append($('<th>'+headers[i]+'</th>'))
        }
        this[name] = $('#cob .'+name).DataTable(
            $.extend(true,{
            "processing" : true,
            "autoWidth": false, 
            "sScrollY":  '300px',
            "bPaginate": false,
            "bJQueryUI": true,
            "bScrollCollapse": true,
            "bAutoWidth": false,
            "sScrollXInner": "100%",
            "sScrollX": true
            },params)
        )
    }

    function Footer(params){
    }

    function Header(params){
    }

    function COB(params){
        // Extend the default parameters
        defaults = {
            'div': $('<div>'),
        }    
        this.params = $.extend(true,defaults,params)

        this.params.div
            .append($('<div>',{class:'graph'}))
            .append($('<div>',{class:'menu'}))
            .append($('<div>',{class:'loci'}))
            .append($('<div>',{class:'footer'}))
            .append($('<div>',{class:'header'}))
        this.graph = new Graph({
            'div' :$("#cob .graph")
        });
        this.menu = new Menu({
            'div' :$('#cob .menu')
        });
        this.menu.add_table('OntologyTable',['Ontology','Description'],{
            "ajax":"Camoco/available_datasets/Ontology",
        })
        this.menu.add_table('TermTable',['Name','Desc','Num SNPs','Num Genes'])
        // loci table
        this.loci = new Menu({
            'div' : $('#cob .loci')
        })
        this.loci.highlighted_rows = []
        this.loci.add_table('LociTable',
            [   'Locus',
                'Transcript',
                'Chr',
                'Start',
                'End',
                'Strand',
                'MaizeSeq.org Annot',
                'BFGR Annot',
                'PFAM Domain',
                'Rice Orth',
                'Rice Annot',
                'Sorghum Orth',
                'Sorghum Annot',
                'Panicum Orth',
                'Panicum Annot',
                'Grassius TF',
                'MapMan',
                'classical',
                'EFP Link',
                'Functional Annot'
            ],
            {}
        )
        this.footer = new Footer({});
        this.header = new Header({});

        // Register top level events
        $('#cob .OntologyTable tbody').on('click','tr', function() {
            var name = $('td',this).eq(0).text();
            cob.menu.TermTable.clear().ajax.url("Ontology/Terms/"+name).load().draw()
            cob.menu.LoadedOntology = name
            $('#cob .OntologyTable .selected').toggleClass('selected')
            $(this).toggleClass('selected')
            cob.menu.OntologyTable.loadedontolgy = name
        });
        $('#cob .TermTable tbody').on('click','tr', function(){
            document.LoadedTerm = $('td',this).eq(0).text();
            $('#cob .TermTable .selected').toggleClass('selected')
            $(this).toggleClass('selected')
            $.getJSON('api/COB/ROOT/'+cob.menu.LoadedOntology+'/'+document.LoadedTerm)
                .done(function(data){
                    a = data
                    console.log('loading data')
                    cob.graph.cy.load(data,
                        function(){
                            console.log('Loading Data')
                        },
                        function(){
                            console.log('Fitting ');
                            cob.load_annotations()
                        })
                    .layout(arbor_options)
                })
                .fail(function(data){
                    console.log("Nopers")
                })
        })
        this.graph.cy.on('click','node',{},function(evt){
            var node = evt.cyTarget
            console.log("CLICKED "+node.id())
            //unhighlight old rows
            cob.loci.LociTable.rows(cob.loci.highlighted_rows)
                .nodes()
                .to$()
                .toggleClass('selected')
            cob.loci.highlighted_rows = cob.loci.LociTable.rows().flatten()
                .filter(function(rowIdx){
                    return cob.loci.LociTable.cell(rowIdx,0).data() == node.id() ? true : false;
                })
            $('.loci .dataTables_scrollBody').scrollTo(
                cob.loci.LociTable.rows(cob.loci.highlighted_rows)
                    .nodes()
                    .to$()
                    .toggleClass('selected')
            )
        })
        var timeout;
        this.graph.cy.on('select',{},function(evt){
            cob.graph.selected = []
            clearTimeout(timeout)
            timeout = setTimeout(function(){
                cob.graph.cy.elements(':selected')
                .filter(function(){return this.isNode()})
                .each(function(){
                    cob.graph.selected.push(this.id()) 
                })
                cob.loci.LociTable.search(cob.graph.selected.join("|"),true).draw()
            },100)
        })

        this.load_annotations = function(){
            // get a list of loaded nodes
            var nodes = this.graph.cy.nodes()
            var node_ids = []
            for(var i=0; i < nodes.length; i++){
                node_ids.push(nodes[i].id())
            }
            this.loci.LociTable.clear()
                .ajax.url("api/Annotations?genes="+node_ids.join(','))
                .load().columns.adjust().draw()
        }
    }
