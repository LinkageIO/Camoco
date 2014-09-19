    function Graph(params){
        defaults = {
            'div':$('<div>'),
        }    
        this.params = $.extend(true,defaults,params)

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
        defaults = {
            'div' : $('<div>',{'class':'menu'})
        }
        this.params = $.extend(true,defaults,params)
        // Add the ontology and term table
        this.params.div
            .append(
                $('<table>',{class:'OntologyTable display',cellspacing:'0',width:'100%'})
                .html('<thead><th>Ontology</th><th>Description</th></thead>')
            )
            .append(
                $('<table>',{class:'TermTable display',cellspacing:'0',width:'100%'})
                .append($('<thead>')
                    .append('<tr>')
                        .append($('<tr>')
                            .append($('<th>Name</th>'))
                            .append($('<th>Desc</th>'))
                            .append($('<th>Num SNPs</th>'))
                            .append($('<th>Num Genes</th>'))
                        )
                )

            )
            .append(
                $('<table>',{class:'LocusTable display',cellspacing:'0',width:'100%'})
                .append($('<thead>')
                    .append('<tr>')
                        .append($('<tr>')
                            .append($('<th>Locus</th>'))
                            .append($('<th>Gene</th>'))
                            .append($('<th>Expressed?</th>'))
                            .append($('<th>Local Degree</th>'))
                            .append($('<th>Global Degree</th>'))
                        )
                )
            )
        // Initialize Actions
        this.ontologytable = $('#cob .menu .OntologyTable').DataTable({
            "ajax":"Camoco/available_datasets/Ontology",
            "processing" : true,
            "paginate" : false,
            'scrollCollapse': true,
            'scrollY' : '300px'
        });
        this.termtable = $('#cob .TermTable').DataTable({
            "processing" : true,
            "paginate" : false,
            'scrollCollapse': true,
            'scrollY' : '300px'
        })
        this.locustable = $('#cob .LocusTable').DataTable({
            "processing" : true,
            "paginate" : false,
            'scrollCollapse': true,
            'scrollY' : '300px'
        })
        $('#cob .OntologyTable tbody').on('click','tr', function() {
            var name = $('td',this).eq(0).text();
            cob.menu.termtable.clear().ajax.url("Ontology/Terms/"+name).load().draw()
            cob.menu.LoadedOntology = name
            $('#cob .OntologyTable .selected').toggleClass('selected')
            $(this).toggleClass('selected')
            cob.menu.ontologytable.loadedontolgy = name
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
                        function(){console.log('Loading Data')},
                        function(){
                            console.log('Fitting ');
                        })
                    .layout(arbor_options)
                })
                .fail(function(data){
                    console.log("Nopers")
                })
        })


    } // End Menu

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
        this.footer = new Footer({});
        this.header = new Header({});

        this.load_base = function(data){
            var hello = 'world'
        }

    }
