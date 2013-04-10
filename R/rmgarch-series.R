#################################################################################
##
##   R package rmgarch by Alexios Ghalanos Copyright (C) 2008-2013.
##   This file is part of the R package rmgarch.
##
##   The R package rmgarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rmgarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################
.extractmdata = function(data, warn = FALSE)
{
	xdata = try(as.xts(data), silent = TRUE)
	if(inherits(xdata, "try-error")){
		if(warn) warning("\nrmgarch-->warning: data indexing not recognized by xts...coercing to Date with origin 1970-01-01.")
		if(is.data.frame(data) | is.matrix(data)) data = as.matrix(data)
		rownames(data) = NULL
		xdata = xts(data, as.POSIXct(as.Date(1:NROW(data), origin="1970-01-01")))
	}
	obj = list()
	obj$data = coredata(xdata)
	obj$index = index(xdata)
	obj$period = median(diff(index(xdata)))
	return(obj)
}
.matrix2xts = function(data){
	data = as.matrix(data)
	return(xts(data, as.POSIXct(as.Date(1:NROW(data), origin="1950-01-01"))))
}

.genxts = function(index0, length.out = 10, period = "days"){
	Z = seq(index0, by = period, length.out=length.out)
	
}
#---------------------------------------------------------------------------------
# lag functions
.embed = function(data, k, by = 1, ascending = FALSE)
{
	# n = no of time points, k = number of columns
	# by = increment. normally =1 but if =b calc every b-th point
	# ascending If TRUE, points passed in ascending order else descending.
	# Note that embed(1:n,k) corresponds to embedX(n,k,by=1,rev=TRUE)
	# e.g. embedX(10,3)
	#if(is.null(dim(data)[1])) n = length(data) else n = dim(data)[1]
	data = matrix(data, ncol = 1)
	n = dim(data)[1]
	s = seq(1, n - k + 1, by)
	lens = length(s)
	cols = if (ascending) 1:k else k:1
	return(matrix(data[s + rep(cols, rep(lens, k)) - 1], lens))
}

.lagx = function(data, n.lag = 1, removeNA = FALSE, pad = NA)
{
	# has NAs
	data = as.matrix(data)
	n = dim(data)[1]
	d = dim(data)[2]
	if(dim(data)[2] == 1) data = matrix(data, ncol = 1)
	z = apply(data, 2, FUN = function(x) .embed(x, n.lag + 1)[, n.lag + 1])
	if(!removeNA) z = rbind(matrix(pad, ncol = d, nrow = n.lag), z)
	return(z)
}

.lagmatrix = function(data, n.lag = 1, pad = 0)
{
	n = length(as.numeric(data))
	z = matrix(NA, ncol = n.lag, nrow = n)
	for(i in 1:n.lag) z[,i] = .lagx(as.numeric(data), i, removeNA = FALSE, pad = pad)
	z = cbind(data, z)
	colnames(z) = paste("lag-", 0:n.lag, sep = "")
	return(z)
}

repmat = function(a, n, m)
{
	kronecker(matrix(1, n, m), a)
}

size = function(x, n = NULL)
{
	x = as.matrix(x)
	if(missing(n)) sol = c(n = dim(x)[1], m = dim(x)[2]) else sol = dim(x)[n]
	return(sol)
}

zeros = function(n = 1, m = 1)
{
	if(missing(m)) m = n
	sol = matrix(0, nrow = n, ncol = m)
	return(sol)
}

ones = function(n = 1, m = 1)
{
	if(missing(m)) m = n
	sol = matrix(1, nrow = n, ncol = m)
	return(sol)
}

newlagmatrix = function(x, nlags, xc)
{
	nlags = nlags+1
	xt = size(x, 1);
	newX = rbind(x, zeros(nlags, 1))
	lagmatrix = repmat(newX, nlags, 1)
	lagmatrix = matrix(lagmatrix[1:(size(lagmatrix,1)-nlags)], nrow = (xt+nlags-1), ncol = nlags)
	lagmatrix = lagmatrix[nlags:xt,]
	y = lagmatrix[,1]
	x = lagmatrix[,2:nlags]
	if(xc == 1) x = cbind(ones(size(x,1), 1), x)
	return(list(y = y, x = x))
}

.colorgradient = function(n = 50, low.col = 0.6, high.col=0.7, saturation = 0.8) {
	if (n < 2) stop("n must be greater than 2")
	n1 = n%/%2
	n2 = - n - n1
	c(hsv(low.col, saturation, seq(1, 0.5, length = n1)),
			hsv(high.col, saturation, seq(0.5, 1, length = n2)))
}

.simlayout = function(m)
{
	if(m == 1){
		nf = c(1, 1, 2, 2)
		nf = layout(matrix(nf, 2, 2, byrow = TRUE), respect = TRUE)
		middle.plot = 1
	}
	if(m == 2){
		nf = c(1, 1, 1, 1, 0, 2, 2, 0, 3, 3, 3, 3)
		nf = layout(matrix(nf, 3, 4, byrow = TRUE), respect = TRUE)
		middle.plot = 2
	}
	if(m == 3){
		nf = c(1, 0, 0, 2, 0, 3, 3, 0, 4, 0, 0, 5)
		nf = layout(matrix(nf, 3, 4, byrow = TRUE), respect = TRUE)
		middle.plot = 3
	}
	if(m == 12){
		nf = c(1, 2, 3, 4,
				5, 6, 6, 7,
				8, 6, 6, 9,
				10, 11, 12, 13)		
		nf = layout(matrix(nf, 4, 4, byrow = TRUE), respect = TRUE)
	}
}

.sdigit = function(x){
	sid = as.numeric(strsplit(format(as.numeric(x), scientific=TRUE), "e")[[1]])[2]
	10^(-sid)
}

.stars = function(testvector, levels = c(0.01, 0.05, 0.1))
{
	N = length(testvector)
	ans = vector(mode="character", length = N)
	#recursive replacement
	z = which(testvector<levels[3])
	ans[z] = c("*")
	z = which(testvector<levels[2])
	ans[z] = c("**")
	z = which(testvector<levels[1])
	ans[z] = c("***")
	ans
}