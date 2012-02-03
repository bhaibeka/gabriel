## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##
## Objects for Generation and Characterization of Gabriel Graphs
## CH - 12.29.04
##
## Version - 0.2	(1.11.2005)
##

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Class: Node
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

setClass("node",
        representation( name="character", label="character", coord="vector", degree="integer",
                        nclass="character", topology="character", n.comp="integer" ),
        prototype(name="node", label="node", coord=numeric(0), degree=0L, nclass=character(0),
                    topology=character(0), n.comp=1L))

# ----------------------------------------------------------------
# 	node				Constructor for class node
# ----------------------------------------------------------------
node <-
function( id, nclass=0, coord=numeric(0), degree=0, label=character(0), topology=character(0), n.comp=1 ) {
	
	# Make sure coordinates are given as a vector and class as a string
	x <- as.vector(coord)
	c <- as.character(nclass)
	id <- as.character(id)
	dg <- as.integer(degree)
	nc <- as.integer(n.comp)
	lb <- as.character(label)
	tp <- as.character(topology)
	
	# Check for valid choices of topology
	
	if( length(tp) != 0 ) {
        top.choices <- c("core", "border", "isolated", "mixed")
		if( !is.element(tp, top.choices) ) {
			w <- paste("Topology set to NULL.\n Valid topology values nodes are:", 
								paste(top.choices,collapse="/"), "\n")
			warning(w)
			tp <- character(0)
		}
	}
	
	val <- new( "node", name=id, coord=x, degree=dg, label=lb, nclass=c, topology=tp, n.comp=nc )
	return(val)
}

# ---------------------------------------------------------------
#	show				Show method for class node
# ---------------------------------------------------------------
setMethod("show", signature("node"),
	function(object) {
		cat("ID:", slot(object, "name"), "\nCoordinates:", slot(object, "coord"), "\nAttributes:\n")
		print( c("class"=slot(object, "nclass"), "topology"= slot(object, "topology"),
                "label"=slot(object, "label")) )
		cat("Degree:", slot(object, "degree"),"\n")
		cat("State:", ifelse(slot(object, "n.comp")==1, "single", paste("composite (", slot(object, "n.comp"),")")))
		cat("\n")
		invisible(object) 
	}
)

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Class: Edge
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

setClass("edge",
        representation(from = "node", to = "node", length = "numeric", nclass = "character",
                        topology = "character", n.comp = "integer", label = "character"),
        prototype(length=0, nclass=character(0), topology=character(0), n.comp=1L, label=character(0)))

# ----------------------------------------------------------------
# 	edge				Constructor for class edge
# ----------------------------------------------------------------

edge <-
function(node1, node2, label=NULL, n.comp=1) {
	
	if( !is.node(node1) || !is.node(node2) ) stop("node objects required.\n")
		
	# Class and topology properties based on corresponding properties of nodes
	cla <- ifelse( node1@nclass == node2@nclass, paste("class", node1@nclass, sep="-"), "tomek" )
	top <- ifelse( node1@topology == node2@topology, node1@topology, "mixed" )

	# Length
	if( length(node1@coord) != 0 && length(node2@coord) != 0 )
		dist <- vecnorm( node1@coord - node2@coord )
	else
		dist <- numeric(0)
		
	val <- new("edge", from=node1, to=node2, length=dist, nclass=cla, topology=top,
					label=as.character(label), n.comp=as.integer(n.comp))
	return(val)
}

# ---------------------------------------------------------------
#	show				Show method for class edge
# ---------------------------------------------------------------

setMethod("show", signature("edge"),
	function(object) {
		cat("Edge:", object@label, "\n")
		cat("Nodes: (", object@from@name, ",", object@to@name, ")\nAttributes:\n")
		print( c("class"=object@nclass, "topology"= object@topology, "length" = object@length) )
		cat("State:", ifelse(object@n.comp==1, "single", paste("composite (",object@n.comp,")")))
		cat("\n")
		invisible(object) 
	}
)


## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Class: 	ggraph
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


setClass("ggraph", representation( vertices = "list", edges = "matrix", edge.list = "list",
				n.v = "integer", n.e = "integer", adj.list = "list", label = "character", type = "character" ))

# ----------------------------------------------------------------
# 	ggraph				Constructor for class ggraph
# ----------------------------------------------------------------
ggraph <-
function( data=NULL, node.class=NULL, vertices=NULL, edges=NULL, edge.list=NULL, label=character(0),
            type="regular", algorithm="brute" ) {
	
	val <- NULL
	
    switch(type,
    
        "composite" = {
    		if( is.null(vertices) || is.null(edge.list) )
    			stop ("Node list and edge list expected for composite graphs" )

    		# Edge matrix
    		edges <- t( sapply(edge.list, FUN=function(x) as.numeric(c(x@from@name, x@to@name))) )

    		# Adjacency list
    		node.list <- sort(unique(edges))
    		a.list <- vector("list")
    		for( s in node.list ) {
    			a.list[[s]] <- sort(c( edges[ edges[,1]==s, 2 ], edges[ edges[,2]==s, 1 ] ))
    			# Degree of node
    			vertices[ sapply(vertices, FUN = function(x, y) x@name == y, y=s) ][[1]]@degree <-
    					length( a.list[[s]] )
    		}

    		val <- new("ggraph", vertices = vertices, edges = edges, edge.list = edge.list,
    						n.v=length(vertices), n.e=nrow(edges), adj.list=a.list, label=as.character(label),
    						type=type)
	   },

       "regular" = {
    		if( !is.null(data) ) {
                data <- as.matrix(data)
    			if( !is.matrix(data) ) stop ("Data matrix expected")
    			cat("Setting up the node objects...\n")
    			vertices <- as.nodes( cbind(data, node.class) )
    			cat("Computing the Gabriel Graph...\n")
    			edges <- switch(algorithm,
                    heuristic = gabriel.2( data ),
                    gabriel.1( data ))
    		}
    		else {
    			if( is.null(edges) ) stop ("Edges not provided")
    			if( ncol(edges) != 2 ) stop ("Edges expected in matrix of ncol=2")
    			if( is.null(vertices) ) stop ("Node list not provided")
    			if( !is(vertices,"list") || !is(vertices[[1]], "node") ) stop ("Vertices not a node list")
    		}

    		# Nodes in graph
    		node.list <- sort(unique(as.vector(edges)))

    		# Check that all node references are included in vertices
    		mis.nodes <- setdiff( node.list, sapply(vertices, FUN=function(x) x@name) )
    		if( length(mis.nodes) > 0 ) {
    			cat("Stopped because following nodes had missing references:", mis.nodes)
    			stop()
    		}

    		n.v <- length(node.list)
    		n.e <- nrow(edges)

    		# Adjacency list
    		cat("Computing adjacency list...\n")
    		a.list <- vector("list")
    		for( s in node.list )
    			a.list[[s]] <- sort(c( edges[ edges[,1]==s, 2 ], edges[ edges[,2]==s, 1 ] ))

    		# Node Topology
    		cat("Updating node topology attributes...\n")
    		class.list <- unique( sapply(vertices, FUN=function(x) x@nclass) )
    		for( s in node.list ) {

    			sc <- sapply( a.list[[s]], FUN=function(x,y) get.nodes(y, x)@nclass, y=vertices )
    			n.same <- length(sc[ sc == class.list[1] ])
    			same.label <- ifelse( n.same == 0 || n.same == length(sc), TRUE, FALSE )
    			if( same.label )
    				vertices[ sapply(vertices, FUN = function(x, y) x@name == y, y=s) ][[1]]@topology <-
    							ifelse( get.nodes(vertices,s)@class == sc[1], "core", "isolated" )
    			else
    				vertices[ sapply(vertices, FUN = function(x, y) x@name == y, y=s) ][[1]]@topology <-
    							"border"

    			# Update degree of node
    			vertices[ sapply(vertices, FUN = function(x, y) x@name == y, y=s) ][[1]]@degree <-
    					length( a.list[[s]] )
    		}

    		# Edge Objects
    		cat("Setting up the edge objects...\n")
    		e.list <- vector("list")
    		for( i in 1:n.e )
    			e.list[[i]] <- edge( get.nodes( vertices, id=edges[i,1] ), get.nodes( vertices, id=edges[i,2] ) )

    		cat("Done.\n")

    		val <- new("ggraph", vertices = vertices, edges = edges, edge.list = e.list, n.v=n.v, n.e=n.e,
    									adj.list=a.list, label=as.character(label), type=type)
  	   },
	   stop ( "Invalid graph type" ))
	
	return( val )
}

# ----------------------------------------------------------------
# 	show				Show method for class ggraph
# ----------------------------------------------------------------

setMethod("show", signature("ggraph"),
function( object ) {
	
	cat("Graph:",object@label, "\nType:", object@type, "\n")
	
	cat("Summary:\n")
	print(c(nodes=object@n.v, edges=object@n.e))
	
	cat("Degree:\n")
	deg.r <- range( sapply(object@vertices, FUN = function(x) x@degree) )
	print( c("min"=deg.r[1], "max"=deg.r[2] ) )
		
	cat("Nodes by class:\n")
	class.list <- sapply(object@vertices, FUN = function(x) x@nclass)
	stats <- sapply( split(class.list, factor(class.list)), length)
	print(stats)
		
	cat("Nodes by topology:\n")
	top.list <- sapply(object@vertices, FUN = function(x) x@topology)
	stats <- sapply( split(top.list, factor(top.list)), length)
	print(stats)
	
	cat("Edges by class:\n")
	class.list <- sapply(object@edge.list, FUN = function(x) x@nclass)
	stats <- sapply( split(class.list, factor(class.list)), length)
	print(stats)
		
	cat("Edges by topology:\n")
	top.list <- sapply(object@edge.list, FUN = function(x) x@topology)
	stats <- sapply( split(top.list, factor(top.list)), length)
	print(stats)

	invisible()
} 
)

# ----------------------------------------------------------------
# 	plot				Plot method for class ggraph
# ----------------------------------------------------------------

setMethod( "plot", signature("ggraph"),
function( x, y = NULL, ... ) {
#-----------------------------------------------------------------------------------------
# USAGE:
# 		plot(x, meth=c("xyproj","mds","gem"), ids=F, normalized=T, add=F, title=F, cex=1.0, 
#				label.offset=0.075)
#-----------------------------------------------------------------------------------------
		
	# Take care of ... arguments and defaults
	
	dot.args <- list(...)
	ids <- dot.args$ids
		if( is.null(ids) ) ids <- FALSE
	normalized <- dot.args$normalized
		if( is.null(normalized) ) normalized <- TRUE
	add <- dot.args$add
		if( is.null(add) ) add <- FALSE
	title <- dot.args$title
		if( is.null(title) ) title <- FALSE
	cex <- dot.args$cex
		if( is.null(cex) ) cex <- 1.0
	dr <- dot.args$label.offset
		if( is.null(dr) ) dr <- 0.075
	meth <- dot.args$method
		if( is.null(meth) ) meth <- "xyproj"
    legend <- dot.args$legend
        if( is.null(legend) ) legend <- TRUE
    leg.text <- dot.args$leg.text
        if( is.null(leg.text) ) leg.text <- sort(unique(sapply(x@vertices, FUN=function(x) x@nclass)))
	
	# Check whether node coordinates were supplied
	
	vert.coord <- t( sapply(x@vertices, FUN=function(x) x@coord) )
	if( diff(dim(vert.coord)) > 0 ) vert.coord <- t(vert.coord)
	mis.coord <- any(sapply(vert.coord, is.null))

    if( ncol(vert.coord) <= 2 && meth == "mds" ) meth <- "xyproj"
    if( ncol(vert.coord) >= 3) meth <- "mds"

	switch( meth,
		mds = {
			require(MASS)
			data.s <- cmdscale( dist( vert.coord ) )
			data.s[, 2] <- -data.s[, 2]
			data <- cbind( 1:x@n.v, data.s )
		},
		gem = {
			t <- 2*pi/x@n.v
			x.vec <- seq(2*pi, t, -t)-3*pi/2
			data <- matrix( c(1:x@n.v, cos(x.vec), sin(x.vec)), ncol=3 )
		},
			data <- cbind( 1:x@n.v, vert.coord )
	)
	
	if( normalized ) {
		oldpar <- par( pty = "s" )
		on.exit( oldpar )
		plot.range <- range(data[,2:3])
		cat(plot.range, "\n")
	}
	
	# Plot the nodes
	class.lab <- as.integer(sapply(x@vertices, FUN=function(x) x@nclass))
    node.pch <- 16
    node.col <- 1:3
    if (!add ) {
		if( normalized )
			plot( data[, 2], data[, 3], xlim=plot.range, ylim=plot.range, axes=FALSE, xlab="", ylab="",
                col = node.col[class.lab], pch = node.pch )
		else
			plot( data[, 2], data[, 3], axes=FALSE, xlab="", ylab="", col = node.col[class.lab], pch = node.pch )
		if( title ) title( main = x@label )
	}
	else
		points( data[, 2], data[, 3], pch = 17, cex = 0.9, col = node.col[class.lab] )
	
	# Add node ids or labels if present
	if ( ids ) {
		# Check whether node labels are provided
		node.labels <- sapply(x@vertices, FUN=function(x) x@label)
		mis.labels <- any( sapply(node.labels, length) == 0 ) 
		if( !mis.labels )
			text( data[, 2], data[, 3] + xy.offset(data[,2], data[,3], dr), labels=node.labels, cex=cex )
		else
			text( data[, 2], data[, 3] + xy.offset(data[,2], data[,3], dr), cex=cex )
	}
	
	# Add the edges
	from.node <- x@edges[,1]
	to.node <- x@edges[,2]
	segments( data[ from.node, 2], data[ from.node, 3], data[ to.node, 2], data[ to.node, 3], col=6 )

    # Add legend
    if( legend )
        legend(x = "topright", pch=node.pch, col=node.col, legend=leg.text, bty="n")
    	
	# Labels for composite graphs 
	if( x@type == "composite" ) {
		l.n <- paste( "(", sapply(x@vertices, FUN=function(x) x@n.comp), ")",sep="" )
		l.p <- sapply(x@vertices, FUN=function(x) x@topology)
		
		if( length( grep("same.class.topology", x@label) ) != 0 )
			comp.node.labels <- paste(l.n,"\n",l.p, sep="")
		else
			comp.node.labels <- l.n 
		comp.edge.labels <- paste( "(", sapply(x@edge.list, FUN=function(x) x@n.comp), ")",sep="" )		
		
		# Number of nodes in component node
		if( ids )
			text( data[,2:3] + xy.offset(data[,2],data[,3],3*dr), labels=comp.node.labels, cex=cex )
		else
			text( data[,2:3] + xy.offset(data[,2],data[,3],dr), labels=comp.node.labels, cex=cex )
			
		# Edges in component edge
		x.c <- 0.5 * ( data[from.node, 2] + data[to.node, 2] )
		y.c <- 0.5 * ( data[from.node, 3] + data[to.node, 3] )
		text( cbind(x.c, y.c) + xy.offset(x.c, y.c, 0.5*dr), labels=comp.edge.labels, col=8, cex=cex )
	}
		
	invisible()
	
} )

# ---------------------------------------------------------------------
#	subgraph				Method to extract subgraphs from ggraph objects
# ---------------------------------------------------------------------

setGeneric( "subgraph", function( x, ...)
		standardGeneric( "subgraph" ) )

setMethod( "subgraph", signature("ggraph"),
function( x, ... ) {
	
	# Take care of ... arguments and defaults
	dot.args <- list(...)
	edit.type <- dot.args$edit.type
	nclass <- dot.args$nclass
	topology <- dot.args$topology
	prune.nodes <- dot.args$prune.nodes
		if( is.null(prune.nodes) ) prune.nodes <- F
	
	edit.choices <- c("boundary", "same.class", "same.topology", "same.class.topology")
	res <- x
	
	# Edit operations on graph edges

	if( is.null(edit.type) ) {
		
		# Check if valid choices of class or topology are provided
	
		class.list <- unique( sapply(x@edge.list, FUN=function(x) x@nclass) )
		#class.list <- class.list[ class.list != "tomek" ]
		top.list <- unique( sapply(x@edge.list, FUN=function(x) x@topology) )
		if( !is.null(nclass) )
			if( !(valid.class <- is.element( nclass, class.list )) ) stop( "Invalid class" )
		if( !is.null(topology) )
			if( !(valid.top <- is.element( topology, top.list )) ) stop( "Invalid topology" )

		if( !is.null(nclass) && !is.null(topology) ) {
			new.edge.list <- get.edges( x@edge.list, nclass=nclass, topology=topology )
			new.label <- paste(x@label, nclass, topology, sep="/")
		}
		else 
			if( !is.null(nclass) ) {
			 	new.edge.list <- get.edges( x@edge.list, nclass=nclass )
				new.label <- paste(x@label, nclass, sep="/")
			}
			else {
				new.edge.list <- get.edges( x@edge.list, topology=topology )
				new.label <- paste(x@label, topology, sep="/")
			}
	}	
	else {
		new.label <- paste( x@label, edit.type, sep="/" )
		switch( edit.type,
		"boundary" = {
			new.edge.list <- get.edges( x@edge.list, nclass="tomek" )
		},
		"same.class" = {
			new.edge.list <- get.edges( x@edge.list, nclass="tomek", mode="exclude" )
		},
		"same.topology" = {
			new.edge.list <- get.edges( x@edge.list, topology="mixed", mode="exclude" )
		},
		"same.class.topology" = {
			new.edge.list <- get.edges( x@edge.list, nclass="tomek", topology="mixed", mode="exclude" )
		},

		{
			w <- paste("Valid edit.types are:", paste(edit.choices,collapse="/"), "\n")
			stop(w)
		} )
	}
	
	# New edge-related attributes	
	new.n.e <- length( new.edge.list )
	new.edge.matrix <- t(sapply( new.edge.list, FUN=function(x) as.numeric(c( x@from@name, x@to@name )) ))
	
	# Adjacency list
	node.list <- sort(unique(new.edge.matrix))
	new.a.list <- vector("list")
	for( s in node.list )
		new.a.list[[s]] <- sort(c( new.edge.matrix[ new.edge.matrix[,1]==s, 2 ], 
									  new.edge.matrix[ new.edge.matrix[,2]==s, 1 ] ))
	
	# Update edited graph attributes
	res@edges <- new.edge.matrix
	res@edge.list <- new.edge.list
	res@adj.list <- new.a.list
	res@n.e <- new.n.e
	res@label <- new.label
	
	# Prune nodes to include only those with live edges
	if( prune.nodes ) {
		res@vertices <- get.nodes( x@vertices, node.list )
		res@n.v <- length( res@vertices )
	}
	
	# Update degree of nodes
	for( s in node.list ) 
			res@vertices[ sapply(res@vertices, FUN = function(x, y) x@name == y, y=s) ][[1]]@degree <-
					length( res@adj.list[[s]] )
	
	return( res ) 	
} )

# ---------------------------------------------------------------------
#	connPlot          Plots connectivity histograms for ggraph objects
# ---------------------------------------------------------------------

setGeneric( "connHist", function( x, ...)
		standardGeneric( "connHist" ) )

setMethod( "connHist", signature("ggraph"),
function( x, ... ) {

    # Take care of ... arguments and defaults
	dot.args <- list(...)
	byclass <- dot.args$byclass
	if( is.null(byclass) ) byclass <- TRUE

	par( mfrow=c(2,2) )
	if( !byclass ) {
		e.len <- sapply(x@edge.list, FUN=function(x) x@length )
		n.deg <- sapply(x@vertices, FUN=function(x) x@degree)
		avg.e <- mean(e.len)
		sd.e <- sd(e.len)
		N.e <- length(e.len)
		avg.n <- mean(n.deg)
		sd.n <- sd(n.deg)
		N.n <- length(n.deg)

		hist( e.len, xlab="Edge Length", main=paste(x@label, "All", sep="/") )
		mtext( paste("( Mean =", format(avg.e, digits=3), ", SD =", format(sd.e,digits=3), ", N =", N.e, ")") )
		hist( n.deg, xlab="Node Degree", main=paste(x@label, "All", sep="/") )
		mtext( paste("( Mean =", format(avg.n, digits=3), ", SD =", format(sd.n,digits=3), ", N =", N.n, ")") )
	}
	else {
		class.list <- sort( unique( sapply(x@vertices, FUN=function(x) x@nclass) ) )
		for( s in class.list ) {
			st <- paste("class",s,sep="-")

			e.len <- sapply( get.edges(x@edge.list, nclass=st ), FUN=function(x) x@length )
			n.deg <- sapply( get.edges(x@vertices, nclass=s ), FUN=function(x) x@degree )
			avg.e <- mean(e.len)
			sd.e <- sd(e.len)
			N.e <- length(e.len)
			avg.n <- mean(n.deg)
			sd.n <- sd(n.deg)
			N.n <- length(n.deg)

			hist( e.len, xlab="Edge Length", main=paste(x@label, st, sep="/") )
			mtext( paste("( Mean =", format(avg.e, digits=3), ", SD =", format(sd.e,digits=3), ", N =", N.e, ")") )
			hist( n.deg, xlab="Node Degree", main=paste(x@label, st, sep="/") )
			mtext( paste("( Mean =", format(avg.n, digits=3), ", SD =", format(sd.n,digits=3), ", N =", N.n, ")") )
		}
	}
} )

# -----------------------------------------------------------------------------
#	gcomponents		Utility function to find connected components in ggraph objects
# -----------------------------------------------------------------------------

setGeneric( "gcomponents", function( x, ...)
		standardGeneric( "gcomponents" ) )

setMethod( "gcomponents", signature("ggraph"),
function( x, ... ) {
# Finds connected components of graph using non-recursive BFS

    # arguments and defaults
    dot.args <- list(...)
	verbose <- dot.args$verbose
	if( is.null(verbose) ) verbose <- FALSE

    nodes <- sapply( x@vertices, FUN=function(x) x@name )
	a.list <- x@adj.list
	component <- vector("list")

	n.components <- 0
	for( u in nodes )
		component[ u ] <- 0
	for( u in nodes )
		if( component[[ u ]] == 0 ) {
			n.components <- n.components + 1
			component[ u ] <- -2
			# Breadth-first search loop
			repeat {
				num.visited <- 0
				for( u in nodes )
					if( component[[ u ]] == -2 ) {
						num.visited <- num.visited + 1
						component[[ u ]] <- n.components
						for( v in a.list[[ as.numeric(u) ]] )
							if( component[[ as.character(v) ]] == 0 )
								component[ as.character(v) ] <- -2
					}
				if( num.visited == 0 ) break
			}
		}

	# Print components
	if(verbose)
        for( i in 1:max(unlist(component)) ) {
		  cc <- names( component[ sapply(component, FUN=function(x,y) x == y, y=i) ] )
		  cat(i, ": {", cc, "}\n")
        }

	# Output components in a list
	components.list <- vector("list")
	for( i in 1:max(unlist(component)) ) {
		cc <- names( component[ sapply(component, FUN=function(x,y) x == y, y=i) ] )
		name <- paste("Component", i, sep="-")
		components.list[[ i ]] <- cc
	}
	names( components.list ) <- paste("Component", 1:max(unlist(component)), sep="-")
	return( components.list )
} )

# -----------------------------------------------------------------------------------
# 	composite.ggraph				Builts composite graph of a given type from a ggraph object
# -----------------------------------------------------------------------------------

setGeneric( "composite.ggraph", function( x, ...)
		standardGeneric( "composite.ggraph" ) )

setMethod( "composite.ggraph", signature("ggraph"),
function( x, ... ) {

    # arguments and defaults
    dot.args <- list(...)
	type <- dot.args$type
	if( is.null(type) ) type <- "class"
	verbose <- dot.args$verbose
	if( is.null(verbose) ) verbose <- FALSE

	res <- NULL

	# Extract subgraph that corresponds to the desired type
	type.choices <- c("class", "topology", "class.topology")
	if( !is.element( type, type.choices ) ) stop ("Invalid type for composite graph")

	sub.x <- subgraph( x, edit.type = paste("same", type, sep=".") )

	# Build composite graph
	res <- comp.graph( sub.x, x@edges, verbose )

	return( res )
} )


# -----------------------------------------------------------------------------
#	comp.graph		Utility function to build composite graph from connected components
# -----------------------------------------------------------------------------

comp.graph <-
# g is the subgraph of which the components should be estimated
# edge.matrix is the edge matrix of the original graph
#                                       --------
function( g, edge.matrix, verbose=TRUE ) {

	# Find connected components
	cc <- gcomponents( g )
	n.comps <- length(cc)

	# Create node objects for composite graph
	node.list <- vector("list")

	for( i in 1:n.comps ) {
		c.coord <- NULL
		c.num <- length( cc[[i]] )
		if( c.num == 1 ) {
			node <- get.nodes( g@vertices, id = cc[[i]] )
			c.class <- node@nclass
			c.top <- node@topology
			c.coord <- node@coord
		}
		else {
			c.nodes <- get.nodes( g@vertices, id=cc[[i]] )

			# Composite class
			classes <- sapply( c.nodes, FUN=function(x) x@nclass )
			c.class <- ifelse( any( classes != classes[1] ), "mixed", classes[1] )

			# Composite topology
			topols <- sapply( c.nodes, FUN=function(x) x@topology )
			c.top <- ifelse( any( topols != topols[1] ), "mixed", topols[1] )

			#Composite coordinates
			c.xyzs <- sapply( c.nodes, FUN=function(x) x@coord )
			if( !any( sapply( c.xyzs, "is.null" ) ) )
				c.coord <- apply( c.xyzs, 1, mean )
		}
        if(verbose)
            cat( "\ncomp -",i, "\nNodes:", c.num, "\nClass:", c.class, "\nTopology:", c.top,
			  "\nCoordinates:", c.coord, "\n" )

		node.list[[i]] <- node( id=i, nclass=c.class, topology=c.top, n.comp=c.num, coord=c.coord )
	}

	# Create edge objects for composite graph
	e.list <- vector( "list" )
	edges <- edge.matrix

	if(verbose) cat( "\nEdges:\n" )
	i.e <- 0
	for( i in 1:(n.comps-1) ) {
		e1 <- edges[ sapply(1:nrow(edges),
						 FUN=function(k,x,y) any(outer(y, c(x[k,1],x[k,2]), "==")), x=edges, y=cc[[i]]), ]

		for( j in (i+1):n.comps ) {
			e2 <- e1[ sapply(1:nrow(e1),
					 FUN=function(k,x,y) any(outer(y, c(x[k,1],x[k,2]), "==")), x=e1, y=cc[[j]]), ]
			if( length(e2) > 0 ) {
				i.e <- i.e + 1
                if(verbose) cat(i, "->", j, "(", ifelse( is.vector(e2), 1, nrow(e2) ),")\n" )

				e.list[[i.e]] <- edge( node1 = node.list[[i]], node2 = node.list[[j]],
												 		n.comp = ifelse( is.vector(e2), 1, nrow(e2) ) )
			}
		}
	}

	# Create ggraph object for composite graph
	res <- ggraph( vertices = node.list, edge.list = e.list, type = "composite", label = g@label )

	return( res )
}

	
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##		Graph Algorithms and Utility Functions
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# -------------------------------------------------------------
# 	gabriel.1			Constructs Gabriel Graph from data table
# -------------------------------------------------------------

# Brute Force Algorithm for Gabriel Graph Construction
# Algorithm is O(d n^3)
# (Bhattacharya et al. 1981)

gabriel.1 <-
# Function to construct Gabriel graph from an array of points
# Brute-force algorithm
# Each row of the input matrix is an N-dimensional vector (point)
function( x ) {

	if ( !is.matrix(x) ) stop ("Matrix is expected as argument.")
	ndim <- ncol(x)
	npts <- nrow(x)

    # Function to calculate euclidian distance between vectors
    # If arguments are vectors returns ||x-y||^(1/2)

    vecnorm2<- function(x)
        sum(x^2)^0.5

	edges <- vector()

	for( i in 1:(npts-1) ) {
		for( j in (i+1):npts ) {
			d1 <- vecnorm2( x[i, ] - x[j, ] )
			sphere.empty <- TRUE
			for( k in (1:npts)[!is.element(1:npts, c(i,j))] ) {
				d2 <- vecnorm2( x[i, ] - x[k, ] )
				d3 <- vecnorm2( x[j, ] - x[k, ] )
				in.circle <- ( d1^2 > d2^2 + d3^2 )
				sphere.empty <- sphere.empty && !in.circle
				#cat(i, j, k, d1, d2, d3, in.circle, sphere.empty, "\n")
				if( in.circle ) break
			}
			if( sphere.empty ) {
				edges <- append( edges, c(i,j) )
				cat(i, "->", j, "\n")
			}
		}
	}
	return( matrix(edges, ncol=2, byrow=T) )
}

# -------------------------------------------------------------
# 	gabriel.2			Constructs Gabriel Graph from data table
# -------------------------------------------------------------

# Heuristic rejection Algorithm for Gabriel Graph Construction
# Algorithm is ~ O(d n^2)
# (Bhattacharya et al. 1981)

gabriel.2 <-
# Function to construct Gabriel graph from an array of points
# Heuristic rejection algorithm
# Each row of the input matrix is an N-dimensional vector (point)
function( x ) {
	
	if ( !is.matrix(x) ) stop ("Matrix expected as argument.")
	ndim <- ncol(x)
	npts <- nrow(x)
	edges <- vector()
	
    vecnorm2<- function(x)
        sum(x^2)^0.5

    for( i in 1:(npts-1) ) {
		# Candidate set of Gabriel neighbors, Ni
		vertices <- (i+1):npts
		# Initialize list of vertices to be excluded from Ni
		excluded <- vector()
		
		for( r in vertices ) {
			# Skip vertices in excluded list
			if ( is.element( r, excluded ) ) next
			d1 <- vecnorm2( x[i, ] - x[r, ] )
			
			for( k in (1:npts)[ !is.element( (1:npts), c(i,r) ) ] ) {
				d2 <- vecnorm2( x[i, ] - x[k, ] )
 				d3 <- vecnorm2( x[r, ] - x[k, ] )
				if( is.element( k, vertices[ !is.element( vertices, excluded ) ] ) )
					if( d2^2 > d1^2 + d3^2 ) {
						excluded <- append( excluded, k )
						#cat(k, "k-excl\n")
					}	
				if( d1^2 > d2^2 + d3^2 ) { 
					excluded <- append( excluded, r )
					#cat(r, "r-excl\n")
					break
				}
				#cat(i, r, k, "-",excluded, "N:", vertices[ !is.element( vertices, excluded ) ],"\n")
			}
		}
		adj <- vertices[ !is.element( vertices, excluded ) ]
		if( length( adj ) > 0 ) {
			cat( i, "->", adj, "\n" )
			for( s in adj )
				edges <- append( edges, c(i,s) )
		}
	}
	return( matrix(edges, ncol=2, byrow=T) )
}
	
# ---------------------------------------------------------------
#	as.nodes			Converts data matrix to a list of node objects
# ---------------------------------------------------------------
 
as.nodes <-
# Takes as input an n x (p+1) matrix of points and creates a list of node objects
# n = number of points; p = dimension of space; p+1 column contains class labels
function( x ) {
	n <- nrow(x)
	p <- ncol(x)-1
    node.list <- vector("list", n)
	for( i in 1:n )
		node.list[[i]] <- node(id = i, coord = x[i, 1:p], nclass = x[i, p+1])
	return( node.list )
}

# ---------------------------------------------------------------
#	is.node			Test relation class node
# ---------------------------------------------------------------

is.node <-
function( x )
	is(x, "node")

# ---------------------------------------------------------------
#	get.nodes			Utility function to retrieve nodes from node list
# ---------------------------------------------------------------

get.nodes <-
function( node.list, id, nclass, topology ) {
	
	if( missing(id) && missing(nclass) && missing(topology) )
			stop ("need to provide a node attribute: ID, class, or topology" )
	if( !missing(id) && (!missing(nclass) || !missing(topology)) )
			stop ("node can be looked up by ID or class/topology - not both" )
		
    res <-
    if( !missing(id) )
		if( length(id) == 1)
			node.list[ sapply(node.list, FUN = function(x, y) x@name == y, y=id) ][[1]]
		else
			node.list[ as.logical(apply( sapply(node.list,
                                FUN = function(x, y) x@name == y, y=id), 2, sum )) ]
	else 
		if( !missing(nclass) )
			if( !missing(topology) )
				node.list[ sapply(node.list,
                        FUN = function(x, y, z) (x@nclass == y && x@topology == z), y=nclass, z=topology) ]
			else
				node.list[ sapply(node.list, FUN = function(x, y) x@nclass == y, y=nclass) ]
		else
			node.list[ sapply(node.list, FUN = function(x, y) x@topology == y, y=topology) ]

	return( res )
}

# ---------------------------------------------------------------
#	get.edges			Utility function to retrieve edges from edge list
# ---------------------------------------------------------------

get.edges <-
function( edge.list, node.id, nclass, topology, mode="in" ) {
	
	if( missing(node.id) && missing(nclass) && missing(topology) )
			stop ("need to provide an edge attribute (class, or topology) or a node ID" )
		
    res <-
    if( !missing(node.id) )
		edge.list[ sapply(edge.list,
					FUN = function(x, y) x@from@name == y | x@to@name == y, y=node.id) ]
	else 
		if( !missing(nclass) )
			if( !missing(topology) )
				if( mode == "exclude" )
					edge.list[ sapply(edge.list, FUN = function(x, y, z)
								(x@nclass != y && x@topology != z), y=nclass, z=topology) ]
				else
					edge.list[ sapply(edge.list, FUN = function(x, y, z)
								(x@nclass == y && x@topology == z), y=nclass, z=topology) ]
			else
				if( mode == "exclude" )
					edge.list[ sapply(edge.list, FUN = function(x, y) x@nclass != y, y=nclass) ]
				else
					edge.list[ sapply(edge.list, FUN = function(x, y) x@nclass == y, y=nclass) ]
		else
			if( mode == "exclude" )
				edge.list[ sapply(edge.list, FUN = function(x, y) x@topology != y, y=topology) ]
			else
				edge.list[ sapply(edge.list, FUN = function(x, y) x@topology == y, y=topology) ]

	return( res )
}

# -----------------------------------------------------------------------------------
# 	xy.offset				Calculated offset for label placing on graph
# -----------------------------------------------------------------------------------

xy.offset <-
# Calculate label x-y offset for two vectors of x and y coordinates
function( x, y, dr ) {

	if( length(x) != length(y) ) stop ("vectors should be of the same length")
	
	phi <- atan2( y, x )
	offset <- dr * matrix( c(cos(phi), sin(phi)), ncol=2, byrow=FALSE )
	return( offset )
	
}

# -----------------------------------------------------------------------------------
# 	veconorm                 Vector norm from package Splus2R
# -----------------------------------------------------------------------------------

vecnorm <- function(x, p=2) {
  if(is.character(p)){
    if(charmatch(p, "maximum", nomatch = 0) == 1)
      p <- Inf
    else if(charmatch(p, "euclidean", nomatch = 0) == 1)
      p <- 2
    else
      stop("improper specification of p")
  }
  if(!is.numeric(x) && !is.complex(x))
    stop("mode of x must be either numeric or complex")
  if(!is.numeric(p))
    stop("improper specification of p")
  if(p < 1)
    stop("p must be greater than or equal to 1")

  x <- ifelse(is.numeric(x), abs(x), Mod(x))

  if(p == Inf)
    return(max(x))
  if(p == 1)
    return(sum(x))

  xmax <- max(x)
  if(!xmax)
    return(xmax)
  x <- x/xmax
  xmax * sum(x^p)^(1/p)
}
