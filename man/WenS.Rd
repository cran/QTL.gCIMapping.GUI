\name{WenS}
\alias{WenS}
\title{The second step of Wen method}
\description{
An efficient multi-locus mixed model framework for the detection of small and linked QTLs in F2}
\usage{
WenS(flag,CriLOD,NUM,pheRaw,Likelihood,flagrqtl,yygg,mx,phe,chr_name,v.map,
     gen.raw,a.gen.orig,d.gen.orig,n,names.insert2,X.ad.tran.data,X.ad.t4,dir)
}
\arguments{
      \item{flag}{random or fix model.}
	\item{CriLOD}{LOD score.}
	\item{NUM}{the number of trait.}
	\item{pheRaw}{raw phenotype matrix .}
	\item{Likelihood}{likelihood function.}
	\item{flagrqtl}{do CIM or not.}
      \item{yygg}{covariate matrix.}
      \item{mx}{raw genotype matrix.}
      \item{phe}{phenotype matrix.}
	\item{chr_name}{chromosome name.}
	\item{v.map}{linkage map matrix.}
	\item{gen.raw}{raw genotype matrix.}
	\item{a.gen.orig}{additive genotype matrix.}
	\item{d.gen.orig}{dominant genotype matrix.}
      \item{n}{number of individual.}
      \item{names.insert2}{linkage map after insert.}
      \item{X.ad.tran.data}{genotype matrix after insert.}
	\item{X.ad.t4}{genotype matrix.}
      \item{dir}{file storage path.}
       }
\author{
Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G2=data(genf2)
WenS(flag=0,CriLOD=2.5,NUM=1,pheRaw=NULL,
Likelihood="REML",flagrqtl=TRUE,yygg=NULL,mx=genf2,phe=NULL,chr_name=NULL,
v.map=NULL,gen.raw=NULL,a.gen.orig=NULL,d.gen.orig=NULL,n=NULL,
names.insert2=NULL,X.ad.tran.data=NULL,X.ad.t4=NULL,dir=NULL)
}


