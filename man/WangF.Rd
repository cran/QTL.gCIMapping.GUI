\name{WangF}
\alias{WangF}
\title{To perform QTL mapping with wang method}
\description{
 Genome-wide Composite Interval Mapping
}
\usage{
WangF(pheRaw,genRaw,mapRaw1,yygg1,flagRIL,cov_en,Population,WalkSpeed,CriLOD)
}
\arguments{
	\item{pheRaw}{phenotype matrix.}
	\item{genRaw}{genotype matrix.}
	\item{mapRaw1}{linkage map matrix.}
	\item{yygg1}{the transformed covariate matrix .}
      \item{flagRIL}{if RIL or not.}
      \item{cov_en}{raw covariate matrix.}
	\item{Population}{population flag.}
	\item{WalkSpeed}{Walk speed for Genome-wide Scanning.(WalkSpeed=1).}
       \item{CriLOD}{Critical LOD scores for significant QTL (CriLOD=2.5).}
       }
\author{
Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G1=data(gen)
P1=data(phe)
M1=data(map)
WangF(phe,gen,map,yygg1=NULL,flagRIL=NULL,cov_en=NULL,Population="DH",1,2.5)
}

