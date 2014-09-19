

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
        this.params.div
            .append($('<div>',{class:'menu'}))
        

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
        this.menu = new Menu({});
        this.footer = new Footer({});
        this.header = new Header({});

        $('#OntologyTable').DataTable({
            "ajax":"Camoco/available_datasets/Ontology",
            "processing" : true,
            "paginate" : false,
            'scrollCollapse': true,
            'scrollY' : '100px'
        });
        $('#OntologyTable tbody').on('click','tr', function() {
            var name = $('td',this).eq(0).text();
            TermTable = $('#TermTable').DataTable()
            TermTable.clear().ajax.url("Ontology/Terms/"+name).load().draw()
            document.LoadedOntology = name
        });
        $('#TermTable').DataTable({
            "processing" : true,
            "paginate" : false,
            'scrollCollapse': true,
            'scrollY' : '100px'
        })
    
        $('#TermTable tbody').on('click','tr', function(){
            document.LoadedTerm = $('td',this).eq(0).text();
            $.getJSON('api/COB/ROOT/'+document.LoadedOntology+'/'+document.LoadedTerm)
            .done(function(data){
                a = data
                console.log('loading data')
                cob.graph.cy.load(data,function(){console.log('Loaded')},function(){console.log('Loaded2')})
                cob.graph.cy.layout(arbor_options)
            })
            .fail(function(data){
                console.log("Nopers")
            })
        })


    }
