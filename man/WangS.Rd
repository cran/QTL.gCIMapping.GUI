\name{WangS}
\alias{WangS}
\title{The second step of wang method}
\description{
 Genome-wide Composite Interval Mapping
}
\usage{
WangS(flag,CriLOD,NUM,pheRaw,chrRaw_name,yygg,mx,phe,chr_name,gen,
mapname,CLO)}
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
      \item{CLO}{Number of CPUs.}
 }
\author{
Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G1=data(gen)
P1=data(phe)
M1=data(map)
W1re<-WangF(phe,gen,map,yygg1=NULL,flagRIL=0,cov_en=NULL,Population="DH",1,2.5)
WangS(flag=1,CriLOD=2.5,NUM=1,pheRaw=phe,W1re$chrRaw_name,W1re$yygg,
W1re$mx,W1re$phe,W1re$chr_name,W1re$gen,W1re$mapname,CLO=1)
}


