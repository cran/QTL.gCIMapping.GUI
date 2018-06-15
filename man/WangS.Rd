\name{WangS}
\alias{WangS}
\title{The second step of wang method}
\description{
 Genome-wide Composite Interval Mapping
}
\usage{
WangS(flag,CriLOD,NUM,pheRaw,chrRaw_name,yygg,mx,phe,chr_name,gen,mapname)}
\arguments{
	\item{flag}{fix or random model.}
	\item{	CriLOD}{LOD score.}
	\item{NUM}{The number of trait.}
	\item{pheRaw}{Raw phenotype matrix.}
      \item{chrRaw_name}{raw chromosome name.}
      \item{yygg}{covariate matrix.}
	\item{mx}{raw genotype matrix.}
	\item{phe}{phenotype matrix.}
      \item{chr_name}{chromosome name.}
      \item{gen}{genotype matrix.}
      \item{mapname}{linkage map matrix.}
 }
\author{
Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
P1=data(phe)
WangS(flag=0,CriLOD=2.5,NUM=1,pheRaw=NULL,chrRaw_name=NULL,yygg=NULL,
mx=NULL,phe,chr_name=NULL,gen=NULL,mapname=NULL)
}


