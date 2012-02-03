## =============================================================================
##  Example File for GraphClasses
##  CH 2.3.2012
## =============================================================================

source("GraphClasses_Source.R")

## ........................................................
## nodes & edges

new("node")
new("node", name="N-1", label="new node", coord=runif(2), nclass="1")
new("edge")
new("edge", from=new("node", name="N1"), to=new("node", name="N2"))


## ........................................................
##  generate normal mixtures
##  first source the function from below (requires mvtnorm)

gm1 <- gaussmix.datagen(npoints=200, k.a=2, k.b=2, w.a=c(0.5,0.5), w.b=c(0.5,0.5),
					mean.a = cbind(c(-0.7,0.3), c(0.3,0.3)),
					mean.b = cbind(c(-0.3,0.7), c(0.4,0.7)),
					cov.a = array(c(diag(rep(0.03,2)),diag(rep(0.03,2))), dim=c(2,2,2)),
					cov.b = array(c(diag(rep(0.03,2)),diag(rep(0.03,2))), dim=c(2,2,2)))

gm2 <- gaussmix.datagen(npoints=200, k.a=2, k.b=2, w.a=c(0.8,0.2), w.b=c(0.3,0.7),
					mean.a = cbind(c(-0.7,0.3), c(0.3,0.3)),
					mean.b = cbind(c(-0.3,0.7), c(0.4,0.7)),
					cov.a = array(c(diag(rep(0.03,2)),diag(rep(0.03,2))), dim=c(2,2,2)),
					cov.b = array(c(diag(rep(0.03,2)),diag(rep(0.03,2))), dim=c(2,2,2)))

## ........................................................
##  compute Gabriel graphs

gg1 <- ggraph(gm1)

gg1@label <- "Gaussian Mix 1"
plot(gg1, title=TRUE)
plot(gg1, title=TRUE, normalized=FALSE)

gg2 <- ggraph(gm2, label="Gaussian Mix 2")
plot(gg2, title=TRUE)

## ........................................................
##  extractor functions

get.nodes(gg1@vertices, topology="isolated")
get.edges(gg1@edge.list, node.id = 1)

## ........................................................
##  filter by edge topology

gg1.s0 <- subgraph(gg1, nclass="class-0")
plot(gg1.s0, title=TRUE)

gg1.s1 <- subgraph(gg1, nclass="class-1")
plot(gg1.s1, title=TRUE)

gg1.s01 <- subgraph(gg1, nclass="tomek")
plot(gg1.s01, title=TRUE)

## ........................................................
##  filter by node topology

gg1.i <- subgraph(gg1, topology="isolated")
plot(gg1.i, title=TRUE)

gg1.b <- subgraph(gg1, topology="border")
plot(gg1.b, title=TRUE)

gg1.b0 <- subgraph(gg1, topology="border", nclass="class-0")
plot(gg1.b0, title=TRUE)

gg1.m <- subgraph(gg1, topology="mixed")
plot(gg1.m, title=TRUE)

## ........................................................
##  connectivity statistics

connHist(gg1.b, FALSE)

## ........................................................
##  connected components in graph

s <- gcomponents(subgraph(gg1, edit.type="same.class"))

## ........................................................
##  construct composite graph

gg1.ct <- composite.ggraph(gg2, type="class.topology")
plot(gg1.ct, ids=TRUE, title=TRUE, cex=.75)





## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Auxiliary Functions for Examples
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

gaussmix.datagen <-
# Generates data from 2-class population in which the class-conditional
# density functions are multivariate Normals
# Inputs - npoints: number of samples to generate
#			  d: dimension of feature space
#          p.class.b: probability of succuess (class 1 in 0/1 classification)
#          k.a, k.b: number of mixture components for class A (0) and B (1)
#			  w.a, w.b: vectors of weights for components of class A, B
#          mean.a, mean.b: d x k.a or d x k.b matrix of mean vectors
#          cov.a, cov.b: d x d x k.a (/k.b) array of covariance matrices
function( npoints, d = 2, p.class.b = 0.5, k.a, k.b, w.a, w.b, mean.a, mean.b, cov.a, cov.b ) {

	if( length(w.a) != k.a || length(w.b) != k.b )
		stop ("Incompatible dimensions for component weights: w.a, w.b")
	if( nrow(mean.a) != d || nrow(mean.b) != d || dim(cov.a)[1] != d || dim(cov.b)[1] != d )
		stop ("Incompatible dimensions in mean and covariance matrices: rows not equal to number of features")
	if( dim(cov.a)[1] != dim(cov.a)[2] || dim(cov.b)[1] != dim(cov.b)[2]  )
		stop ("Incompatible dimensions in covariance matrices: Covariance matrix not square")
	if( ncol(mean.a) != k.a || ncol(mean.b) != k.b || dim(cov.a)[3] != k.a || dim(cov.b)[3] != k.b )
		stop ("Incompatible dimensions in mean and covariance matrices: columns not equal to number of components")

    require(mvtnorm)

    # Generate class-conditional distributions

	# --- Class A (0)
	s <- sample(0:(k.a-1), npoints, T, w.a)			# multinomial to pick mixture component
	d.a <- 0
	for( i in 1:k.a ) {
		r <- rmvnorm( n=npoints, mean=mean.a[,i], sigma=cov.a[,,i])
		w <- ifelse( s == (i-1), 1, 0 )
		d.a <- d.a + w * r
	}
	# --- Class B (1)
	s <- sample(0:(k.b-1), npoints, T, w.b)			# multinomial to determine mixture component
	d.b <- 0
	for( i in 1:k.b ) {
		r <- rmvnorm( n=npoints, mean=mean.b[,i], sigma=cov.b[,,i])
		w <- ifelse( s == (i-1), 1, 0 )
		d.b <- d.b + w * r
	}

	# Bernoulli to pick component (p is the weight of class B)
	w.b <- sample( 0:1, npoints, T, c(1-p.class.b,p.class.b) )

	# Mixture joint population density
	dist.AB <- cbind( "x.val" = (1-w.b) * d.a + w.b * d.b, "component" = w.b )

	return(dist.AB)
}
