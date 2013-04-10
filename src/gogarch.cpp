/*################################################################################
##
##   R package rmgarch by Alexios Ghalanos Copyright (C) 2009-2013
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
#################################################################################*/
// Extra method for speedup
#include "gogarch.h"
using namespace Rcpp;

SEXP gogarchSigma(SEXP S, SEXP A)
{
	try{
		Rcpp::NumericMatrix RS(S);
		Rcpp::NumericMatrix RA(A);
		int m = RS.ncol(), n = RS.nrow(), mm = RA.ncol(), mn = RA.nrow(), i;
		arma::mat AS(RS.begin(), n, m, true);
		arma::mat AA(RA.begin(), mn, mm, true);
		arma::mat S(n, m);
		for(i=0;i<n;i++){
			S.row(i) = arma::trans(arma::diagvec(AA * arma::diagmat(AS.row(i))*AA.t()));
		}
		return wrap( S );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->gogarch extractor c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP gogarchCov(SEXP S, SEXP A)
{
	try{
		Rcpp::NumericMatrix RS(S);
		Rcpp::NumericMatrix RA(A);
		int m = RS.ncol(), n = RS.nrow(), mm = RA.ncol(), mn = RA.nrow(), i;
		arma::mat AS(RS.begin(), n, m, true);
		arma::mat AA(RA.begin(), mn, mm, true);
		arma::cube S(m, m, n);
		for(i=0;i<n;i++){
			S.slice(i) = AA * arma::diagmat(AS.row(i))*AA.t();
		}
		return wrap( S );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->gogarch extractor c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP gogarchCor(SEXP S, SEXP A)
{
	try{
		Rcpp::NumericMatrix RS(S);
		Rcpp::NumericMatrix RA(A);
		int m = RS.ncol(), n = RS.nrow(), mm = RA.ncol(), mn = RA.nrow(), i;
		arma::mat AS(RS.begin(), n, m, true);
		arma::mat AA(RA.begin(), mn, mm, true);
		arma::cube S(m, m, n);
		for(i=0;i<n;i++){
			arma::mat tmp1 = AA * arma::diagmat(AS.row(i))*AA.t();
			arma::mat tmp2 = arma::diagmat(1/arma::sqrt(diagvec(tmp1)));
			S.slice(i) = tmp2 * tmp1 * tmp2.t();
		}
		return wrap( S );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->gogarch extractor c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP gogarchCS(SEXP S, SEXP A)
{
	try{
		Rcpp::NumericMatrix RS(S);
		Rcpp::NumericMatrix RA(A);
		int m = RS.ncol(), n = RS.nrow(), mm = RA.ncol(), mn = RA.nrow(), i;
		arma::mat AS(RS.begin(), n, m, true);
		arma::mat AA(RA.begin(), mn, mm, true);
		arma::cube S(m, m*m, n);
		arma::mat KA = kron(AA.t(), AA.t());
		for(i=0;i<n;i++){
			arma::mat tmp = coskewind(AS.row(i));
			S.slice(i) = AA * tmp * KA;
		}
		return wrap( S );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->gogarch extractor c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP gogarchCK(SEXP K, SEXP S, SEXP A){
	try{
		Rcpp::NumericMatrix RK(K);
		Rcpp::NumericMatrix RS(S);
		Rcpp::NumericMatrix RA(A);
		int m = RS.ncol(), n = RS.nrow(), mm = RA.ncol(), mn = RA.nrow(), i, j;
		arma::mat AS(RS.begin(), n, m, true);
		arma::mat AK(RK.begin(), n, m, true);
		arma::mat AA(RA.begin(), mn, mm, true);
		arma::cube S(m, m*m*m, n);
		int m2 = m*m;
		int m3 = m*m*m;
		arma::vec cc1(m3);
		arma::vec cc2(m2);
		arma::vec cc3(m3);
		arma::vec cc4(m);
		arma::vec N(m3);
		arma::uvec idx(m);
		for(i=0;i<m;i++){
			cc1.rows((i*m2), (i+1)*m2-1) = arma::ones<arma::vec>(m2)*i;
			cc2.rows((i*m), (i+1)*m-1) = arma::ones<arma::vec>(m)*i;
			cc4(i) = i;
			idx(i) = (i)*m3 + (i)*m2 + (i)*m + i;
		}
		for(i=0;i<m3;i++){
			N(i)=i;
		}
		cc2 = arma::repmat(cc2,  m, 1);
		cc3 = arma::repmat(cc4, m2, 1);
		arma::mat Y = arma::join_rows(cc1,cc2);
		Y = arma::join_rows(Y,cc3);
		int my = Y.n_cols;
		arma::mat X(m3, my);
		for(i=0;i<m3;i++){
			X.row(i) = sort(Y.row(i));
		}
		arma::uvec X1 = (X.col(0)==X.col(1))%(X.col(1)==X.col(2));
		arma::uvec X2 = (X.col(0)<X.col(1))%(X.col(1)<X.col(2));
		arma::uvec excl1 = arma::find(X1==1);
		arma::uvec excl2 = arma::find(X2==1);
		arma::uvec excl = arma::unique(arma::join_cols(excl1, excl2));
		int nb = excl.n_rows;
		arma::uvec C(m3);
		C.ones();
		for(i=0;i<nb;i++){
			C(excl(i))=0;
		}
		arma::uvec incl = arma::find(C==1);
		N = N.elem(incl);
		arma::mat Z = X.rows(incl);
		int nx = Z.n_rows;
		arma::uvec rw(nx);
		arma::mat AAA = kron(AA.t(), kron(AA.t(), AA.t()));
		arma::mat ck(m, m3);
		for(j=0;j<n;j++){
			ck.zeros();
			arma::rowvec sx = AS.row(j);
			arma::rowvec kx = AK.row(j);
			for(i=0;i<nx;i++){
				if(Z(i,0)==Z(i,1)){
					rw(i) = Z(i,2);
				} else{
					rw(i) = Z(i,0);
				}
				arma::rowvec uq = unique(Z.row(i));
				ck(rw(i), N(i)) = sx(uq(0))*sx(uq(1));
			}
			for(i=0;i<m;i++){
				ck(idx(i)) = kx(i);
			}
			S.slice(j) = AA * ck * AAA;
		}
		return wrap( S );
		} catch( std::exception &ex ) {
			forward_exception_to_r( ex );
		} catch(...) {
			::Rf_error( "rmgarch-->gogarch cokurt c++ exception (unknown reason)" );
		}
		return R_NilValue;
}

SEXP Cov2Cor(SEXP YY, SEXP dimm){
	try{
		Rcpp::NumericMatrix Y(YY);
		Rcpp::IntegerVector dim(dimm);
		arma::cube AY(Y.begin(), dim[0], dim[1], dim[2]);
		arma::cube AX(dim[0], dim[1], dim[2]);
		int i;
		arma::mat tmp(dim[0], dim[1]);
		for(i=0;i<dim[2];i++){
			tmp = arma::diagmat(1/arma::sqrt(diagvec(AY.slice(i))));
			AX.slice(i) = tmp * AY.slice(i) * tmp.t();
		}
		return wrap( AX );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->Cov2Cor c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP ArrayDiag(SEXP YY, SEXP dimm){
	try{
		int i;
		Rcpp::NumericVector Y(YY);
		Rcpp::IntegerVector dim(dimm);
		arma::cube AY(Y.begin(), dim[0], dim[1], dim[2]);
		arma::mat AX(dim[1],dim[2]);
		for(i=0;i<dim[2];i++){
			AX.col(i) =  arma::diagvec(AY.slice(i));
		}
		AX = AX.t();
		return(wrap(AX));
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
	::Rf_error( "rmgarch-->ArrayDiag c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP Cov2Res(SEXP YY, SEXP ZZ, SEXP dimm){
	try{
		int i;
		Rcpp::NumericVector Y(YY);
		Rcpp::IntegerVector dim(dimm);
		Rcpp::NumericMatrix Z(ZZ);
		arma::cube AY(Y.begin(), dim[1], dim[1], dim[0]);
		arma::mat AZ(Z.begin(), dim[0], dim[1]);
		arma::mat AR(dim[0], dim[1]);
		for(i=0;i<dim[0];i++){
			AR.row(i) =  AZ.row(i) * diagmat(arma::sqrt(diagvec(AY.slice(i))));
		}
		return(wrap(AR));
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
	::Rf_error( "rmgarch-->Cov2Res c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP RowApplySort(SEXP YY){
	try{
		Rcpp::NumericMatrix Y(YY);
		int m = Y.ncol(), n = Y.nrow(), i;
		arma::mat AY(Y.begin(), n, m, true);
		arma::mat AX(n, m);
		for(i=0;i<n;i++){
			AX.row(i) = sort(AY.row(i));
		}
		return wrap( AX );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->RowApplySort c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP RowUnique(SEXP YY){
	try{
		Rcpp::NumericMatrix Y(YY);
		int m = Y.ncol(), n = Y.nrow(), i;
		arma::mat AY(Y.begin(), n, m, true);
		arma::mat AX(n, 2);
		for(i=0;i<n;i++){
			AX.row(i) = unique(AY.row(i));
		}
		return wrap( AX );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->RowUnique c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

arma::mat coskewind(arma::rowvec skew){
	int n = skew.n_cols, i;
	int nn = n*n;
	arma::mat S(n, nn);
	S.zeros();
	for(i=0;i<n;i++){
		int indx = (i)*nn + (i)*n + i;
		S(indx) = skew(i);
	}
	return(S);
}
