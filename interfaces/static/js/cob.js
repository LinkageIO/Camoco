    function SavedLayout(params){
        
    }

    function Graph(params){
        defaults = {
            'div':$('<div>'),
        }    
        this.params = $.extend(true,defaults,params)
        this.selected = []
        // set default filters
        this.edge_filter = 3 // edge score
        this.node_filter = 1 // degree
        
        this.params.div
        .append($('<div>',{class:'cy'}))
        .append($('<ul>',{class:'graph_controls'})
            .append($('<li>')
                .append($('<button>Fit</button>',{})
                    .on({'click':function(){cob.graph.cy.fit()}})
                )
            )
            .append($('<li>')
                .append('<span>').html('Edge Filter')
                    .append($('<input>',{'value':3})
                        .on({'change':function(){
                            cob.graph.edge_filter = this.value
                            cob.graph.filter()
                        }})
                    )
            )
            .append($('<li>')
                .append('<span>').html('Degree Filter')
                    .append($('<input>',{'value':0})
                        .on({'change':function(){
                            cob.graph.node_filter = this.value
                            cob.graph.filter()
                        }})
                    )
            )

        )

        this.cy = cytoscape(options = {
            container : this.params.div.find('.cy')[0],    
            // General Options 
            hideEdgesOnViewport: true,
            hideLabelsOnViewport: true,
            minZoom: 1e-50,
            maxZoom: 1e50,
            zoomingEnabled: true,
            userZoomingEnabled: true,
            panningEnabled: true,
            userPanningEnabled: true,
            boxSelectionEnabled: true,
            // Style
            style: cytoscape.stylesheet()
            .selector('node')
            .css({
                'background-color': '#144566',
                'border-width': 1,
                'border-color': '#000',
                'height' : 'mapData(gdegree,0,100,50,10)',
                'width' : 'mapData(gdegree,0,100,50,10)',
                'content':'data(id)',
                'text-halign':'right',
                'font-size' : '12pt',
                'min-zoomed-font-size':1
            })
            .selector(':selected')
            .css({
                'border-width': 10,
                'border-color': '#09BA00'
            })
            .selector('.neighbors')
            .css({
                'border-width': 10,
                'border-color': '#BA0900'
            })
            .selector('.highlighted')
            .css({
                'border-width': 10,
                'border-color': '#0900BA'
            })
            .selector('edge')
            .css({
                'opacity': '0.50',
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
        this.filter = function(){
            // filter values need to be stored in the graph object because
            // interface programming sucks.
            try{
                this.rem_edges.restore() 
                this.rem_nodes.restore()
                this.rem_con_edges.restore()
            }
            catch(e){
                if(e instanceof TypeError){}
            }
            // Get the edges under the filter
            this.rem_edges = this.cy.collection('edge[score <  '+this.edge_filter+']').remove()
            var nodes = this.cy.collection('node[[degree < '+this.node_filter+']]')
            this.rem_con_edges = nodes.connectedEdges().remove()
            this.rem_nodes = nodes.remove()
        }
    }

    function Bar(params){
        defaults = {
            'name' : 'yuckyuck',
            'div'  :  $('<div>',{class:'bar'})
        }
        this.params = $.extend(true,defaults,params)
    }

    function Tab(params){
        defaults = {
            'name' : 'nope',
            'div' : $('<div>',{class:'tab'})
        }
        this.params = $.extend(true,defaults,params)
    
        this.update_table = function(params){
            // delete oldtable
            this[params.name].destroy()
            $('#cob .'+params.name).empty().append(
                $('<thead>').append($('<tr>'))
            )
            a=1
            for (var i=0; i<params.header.length; i++){
                $('#cob table.'+params.name+' thead tr')
                .append($('<th>'+params.header[i]+'</th>'))
            }
            this[params.name] = $('#cob .'+params.name).DataTable(
                $.extend(true,{
                "processing" : true,
                "autoWidth": true, 
                "bPaginate": false,
                "bJQueryUI": true,
                "bScrollCollapse": true,
                "bAutoWidth": true,
                "sScrollXInner": "100%",
                "sScrollX": '100%',
                "sScrollY":  '100%',
                },params)
            )
        }
        this.add_table = function(params){
            this.params.div
                .append(
                    $('<table>',{class: params.name+' display',cellspacing:'0',width:'100%'})
                    .append($('<thead>')
                            .append($('<tr>')
                            )
                    )
                )
            for (var i=0; i<params.header.length; i++){
                this.params.div.find('table.'+params.name+' thead tr')
                .append($('<th>'+params.header[i]+'</th>'))
            }
            this[params.name] = $('#cob .'+params.name).DataTable(
                $.extend(true,{
                "processing" : true,
                "autoWidth": true, 
                "bPaginate": false,
                "bJQueryUI": true,
                "bScrollCollapse": true,
                "bAutoWidth": true,
                "sScrollXInner": "100%",
                "sScrollX": '100%',
                "sScrollY":  '100%',
                },params)
            )
        }

    } // End Tab
    
    function Menu(params){
        // Initialize Tables
        this.tabs = []
        this.handle = []
        //this.tables = []
        defaults = {
            'div' : $('<div>',{'class':'menu'}),
            'header' : $('<ul>',{'class':'header'}),
            'tabs' : $('<div>',{'class':'tabs'})
        }
        this.params = $.extend(true,defaults,params)
        // Add the ontology and term table
        this.params.div.append(this.params.header)
        this.params.div.append(this.params.tabs)
        this.tabwidth = this.params.div.width()
        var that = this

        this.params.header.on('click','li',function(){
            that.show_tab.call(that,$(this).index())
        });

        this.show_tab = function(tabname){
            var index = undefined
            if (!isNaN(parseFloat(tabname)) && isFinite(tabname) && tabname < this.tabs.length){
                // is tab index, return that index 
                index = tabname
                this.params.tabs.css('left','-'+(this.tabwidth*index)+'px')
            }
            for(var i=0; i< this.tabs.length;i++){
                // check each tab for tabname
                if (this.tabs[i].params.name == tabname){
                    index = i
                    this.params.tabs.css('left','-'+(this.tabwidth*index)+'px')
                }
            }
            $('#cob .menu li.selected').toggleClass('selected')
            $(this.params.header.children()[index]).toggleClass('selected')
            return undefined
        }
   
        this.add_tab = function(tab){
            /*  
                This function adds a tab to the menu div.
            */
            // make the tab as wide as the tab section
            tab.params.div.css('width',this.tabwidth+'px')
            this.params.tabs.css('width',(this.tabs.length+1)*this.tabwidth+'px')
            // Append the Tabs name to the header section
            this.params.header.append($('<li>').html(tab.params.name))
            this.params.tabs.append(tab.params.div)
            this.tabs.push(tab)
        }

        this.get_tab = function(tabname){
            ///
            if (!isNaN(parseFloat(tabname)) && isFinite(tabname) && tabname < this.tabs.length)
                // is tab index, return that index 
                return this.tabs[tabname]
            for(var i=0; i< this.tabs.length;i++){
                // check each tab for tabname
                if (this.tabs[i].params.name == tabname)
                    return this.tabs[i]
            }
            return undefined
        }

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
            .append($('<div>',{class:'footer'}))
            .append($('<div>',{class:'header'}))
        this.graph = new Graph({
            'div' :$("#cob .graph")
        });
        this.menu = new Menu({
            'div' :$('#cob .menu')
        });


        this.menu.add_tab(new Tab({'name':'Dataset'}))
        this.menu.add_tab(new Tab({'name':'Network'}))
        this.menu.add_tab(new Tab({'name':'Genes'}))
        // Choose the first tab
        this.menu.show_tab(0)

        this.menu.get_tab("Dataset").add_table({
            "name" : 'OntologyTable',
            "header" : ['Ontology','Description'],
            "ajax":"Camoco/available_datasets/Ontology",
            'sScrollY': this.menu.params.div.innerHeight()/4
        })
        this.menu.get_tab('Dataset').add_table({
            "name":'TermTable',
            'header':['Name','Desc','Num SNPs','Num Genes'],//,'Root Genes'],
            'sScrollY': this.menu.params.div.innerHeight()/4
        })
        this.menu.get_tab('Network').add_table({
            "name":'NetworkTable',
            "header":['Network','Description'],
            "ajax":"Camoco/available_datasets/Expr",
            'sScrollY':this.menu.params.div.innerHeight()/4
        })

        this.menu.get_tab('Genes').highlighted_rows = []
        this.menu.get_tab('Genes').add_table({
            'name': 'LociTable',
            'header': ['Locus',
                    'Chr',
                    'Start',
                    'End',
                    'Strand',
                    'Global Degree',
                    'Term SNPs'
                ],
            'sScrollY': this.menu.params.div.innerHeight()/2,
        })
        // Fix the gene column in the table
        this.menu.get_tab('Genes').FixedColumn = new $.fn.dataTable.FixedColumns(this.menu.get_tab('Genes').LociTable)

        
        this.footer = new Footer({});
        this.header = new Header({});

        // Register top level events
        $('#cob .OntologyTable tbody').on('click','tr', function() {
            var name = $('td',this).eq(0).text();
            cob.menu.get_tab('Dataset').TermTable.clear().ajax.url("Ontology/Terms/"+name).load().draw()
            cob.menu.LoadedOntology = name
            $('#cob .OntologyTable .selected').toggleClass('selected')
            $(this).toggleClass('selected')
        });
        $('#cob .TermTable tbody').on('click','tr', function(){
            // Load the available networks for the Term
            cob.menu.LoadedTerm = $('td',this).eq(0).text();
            $('#cob .TermTable .selected').toggleClass('selected')
            $(this).toggleClass('selected')
            cob.menu.show_tab('Network')
        })
        $('#cob .NetworkTable tbody').on('click','tr',function(){
            $('#cob .NetworkTable .selected').toggleClass('selected')
            $(this).toggleClass('selected')
            cob.menu.LoadedNetwork = $('td',this).eq(0).text();
            $.getJSON('api/COB/'+cob.menu.LoadedNetwork+'/'+cob.menu.LoadedOntology+'/'+cob.menu.LoadedTerm)
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
                    .layout(cola_options)
                })
                .fail(function(data){
                    console.log("Nopers")
                })
        })
        $('#cob .LociTable tbody').on('click','tr',function(){
            $(this).toggleClass('selected')
            gene = $('td',this).eq(0).text();
            cob.graph.cy.center(
                cob.graph.cy.nodes('node[id="'+gene+'"]').select()
            )
        })
        this.graph.cy.on('cxttap','node',{},function(evt){
            var node = evt.cyTarget
            node.toggleClass('highlighted')
        })
        this.graph.cy.on('click',function(evt){
            if(evt.cyTarget == cob.graph.cy){
                //remove all non-sticky decorators
                $('#cob .LociTable .selected').toggleClass('selected')
                cob.graph.cy.$('.neighbors').removeClass('neighbors')
                //cob.loci.LociTable.search(cob.graph.selected.join("|"),true).draw()
            }
        })
        this.graph.cy.on('click','node',{},function(evt){
            var node = evt.cyTarget
            console.log("CLICKED "+node.id())
            // highlight neighbors
            //unhighlight old rows
            a =1 
            cob.menu.get_tab('Genes').LociTable.rows(cob.menu.get_tab('Genes').highlighted_rows)
                .nodes()
                .to$()
                .toggleClass('selected')
            // highlight new rows
            cob.menu.get_tab('Genes').highlighted_rows = cob.menu.get_tab('Genes')
                .LociTable.rows().flatten()
                .filter(function(rowIdx){
                    return cob.menu.get_tab('Genes').LociTable.cell(rowIdx,0).data() == node.id() ? true : false;
                })
            // scroll to row
            $('#cob .dataTables_scrollBody .LociTable').scrollTo(
                cob.menu.get_tab('Genes').LociTable.rows(cob.menu.get_tab('Genes').highlighted_rows)
                    .nodes()
                    .to$()
                    .toggleClass('selected')
            )
            cob.graph.cy.$('.neighbors').removeClass('neighbors')
            node.neighborhood().addClass('neighbors')
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
                //cob.loci.LociTable.search(cob.graph.selected.join("|"),true).draw()
            },100)
        })
        $.contextMenu({
            selector : '.LociTable th',
            callback : function(key,options){
                a=1
            },
            items : {
                "relabel" : {name:'ReLabel Genes'} 
            }
        })

        this.load_annotations = function(){
            // get a list of loaded nodes
            var nodes = this.graph.cy.nodes()
            var node_ids = []
            for(var i=0; i < nodes.length; i++){
                node_ids.push(nodes[i].id())
            }
            try{
                // Fetch ajax
                $.getJSON("api/Annotations/"
                    +cob.menu.LoadedNetwork
                    +"/"+cob.menu.LoadedOntology
                    +"/"+cob.menu.LoadedTerm+"?genes="
                    +node_ids.join(','))
                .success(
                    function(data){
                        data.name = 'LociTable'
                        data['sScrollY'] = cob.menu.params.div.innerHeight()-200
                        cob.menu.get_tab('Genes').update_table(data)
                    }
                )
            }
            catch(err){
                
            }
        }
    }
