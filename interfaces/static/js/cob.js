
var cob = (function(cob, $, undefined){

    function Graph(params){
        params.container = params.container || $('<div>')
    }

    function Menu(params){
        this.container = $(params.container) || $('<div>')
    }

    function Footer(params){
        this.container = $(params.container) || $('<div>')
    }

    function Header(params){
        this.container = $(params.container) || $('<div>')
    }

    function COB(params){
        // Extend the default parameters
        defaults = {
            'containerId':$('<div>'),
        }
        this.params = $.extend(true,params,defaults)
    
        this.graph = new Graph({});
        this.menu = new Menu({});
        this.footer = new Footer({});
        this.header = new Header({});
    }


    return COB

}(window.cob = window.cob || {}, jQuery))
