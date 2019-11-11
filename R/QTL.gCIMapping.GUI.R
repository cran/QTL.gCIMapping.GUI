QTL.gCIMapping.GUI<-function(){
  
  server<-(function(input,output,session){
    
    options(shiny.maxRequestSize=-1)
    
    datasetex<-reactive({
      dataex<-data.frame(marker=c("RGA3(1)","wPt-6358","Hplc2","...", "gwm437", "gwm121","wmc157","","",""),
                         chr=c("1","1","1","...","21","21","21","trait1","trait2","Covar"),
                         pos=c("0","3.034","8.8291","...","162.5218","180.2878","197.9196","T19","T191","CovarName"),
                         DH6.10=c("B","B","A","...","A","A","A","75.33","74","A"),  
                         DH6.101=c("-","-","A","...","B","B","B","105","105.68","B"), 
                         DH6.102=c("B","-","B","...","-","-","A","96.33","97.16","B") 
      )
    })
    
    output$dataexample<-renderTable(datasetex())
    
    
    codeex<-reactive({
      coex<-data.frame(Genotype=c("AA","Aa","aa","AA+Aa(Not aa)", "Aa+aa(Not AA)", "Missing"),
                       Code=c("A","H","B","D","C","-"),
                       Meaning=c("Homozygous genotype (P1)","Heterozygous genotype (F1)","Homozygous genotype (P2)","Dominance to P1","Dominance to P2"," Missing or unclear genotype")
      )
    })
    
    output$codeexample<-renderTable(codeex())
    
    observeEvent(input$parSe, {
      updateTabsetPanel(session, "inTabset",
                        selected="PA"
                        
      )
    })
    
    
    observeEvent(input$drawPl, {
      updateTabsetPanel(session, "inTabset",
                        selected="FIG"
      )
    })
    
    
    upda1<-observeEvent(input$Resolution2,{
      
      if(input$Resolution2=="General resolution"){
        widthG<-1500
        heightG<-600
        pointG<-12
        ppiG<-72 
      }else if(input$Resolution2=="High resolution"){
        widthG<-10000
        heightG<-4000
        pointG<-12
        ppiG<-300 
      }else if(input$Resolution2=="Set by yourself"){
        widthG<-0
        heightG<-0
        pointG<-0
        ppiG<-0
      }
      updateTextInput(session, "widthGCIM", value=widthG)
      updateTextInput(session, "heightGCIM", value=heightG)
      updateTextInput(session, "pointGCIM", value=pointG)
      updateTextInput(session, "ppiGCIM", value=ppiG)
    }) 
    
    
    manual<-eventReactive(input$manl,{
      RShowDoc("Instruction",package="QTL.gCIMapping.GUI") 
    })
    output$manll<-renderUI(manual())
    
    rawGCIM<-reactive({
      req(input$fileDataset)
      genRaw1=NULL;pheRaw=NULL;mapRaw11=NULL;cov_en=NULL
      genRaw<-fread(input$fileDataset$datapath,header = FALSE,stringsAsFactors=T)
      titlenameGen<-genRaw[1,1:3]
      hapName<-c("marker","chr","pos")
      if(all(titlenameGen==hapName)==FALSE){
        showModal(modalDialog(title = "error", "please check the Linkage map's name in file!", easyClose = TRUE))
      }
      traitloc<-which(genRaw[,2]=="trait1")
      
      if(length(traitloc)==0){
        showModal(modalDialog(title = "error", "please check the phenotype in file!", easyClose = TRUE))
        
      }
      envirloc<-which(genRaw[,2]=="Covar1")
      if(length(envirloc)!=0){
        pheRaw<-t(rbind(genRaw[1,][,-(1:2)],genRaw[traitloc:(envirloc-1),][,-(1:2)]))
        cov_en<-cbind(t(genRaw[1,][,-(1:2)]),t(genRaw[envirloc:nrow(genRaw),][,-(1:2)]))
        colnames(cov_en)<-NULL;rownames(cov_en)<-NULL
      }else{
        pheRaw<-t(rbind(genRaw[1,][,-(1:2)],genRaw[traitloc:nrow(genRaw),][,-(1:2)]))  
      }
      colnames(pheRaw)<-NULL;rownames(pheRaw)<-NULL
      genRaw1<-as.matrix(genRaw[1:(traitloc-1),-c(2,3)])
      mapRaw11<-as.matrix(genRaw[1:(traitloc-1),1:3])
      output<-list(genRaw1,pheRaw,mapRaw11,cov_en) 
    })
    
    rawICIM<-reactive({
      req(input$fileDataset)
      geoo=NULL;pho=NULL;poss=NULL;parm=NULL;cov_en=NULL
      geoo<-as.matrix(read.xlsx(input$fileDataset$datapath,sheet = "Genotype",colNames = FALSE))
      pho<-as.matrix(read.xlsx(input$fileDataset$datapath,sheet = "Phenotype",colNames = FALSE))
      poss<-read.xlsx(input$fileDataset$datapath,sheet = "LinkageMap",colNames = FALSE)
      parm<-read.xlsx(input$fileDataset$datapath,sheet = "GeneralInfo",colNames = FALSE)
      if(is.null(input$fileCov)==FALSE){
        cov_en<-fread(input$fileCov$datapath,header = FALSE,stringsAsFactors=T)
        cov_en<-as.matrix(cov_en)
      }
      output<-list(geoo,pho,poss,parm,cov_en)
    })
    
    rawWIN<-reactive({
      req(input$fileDataset)
      y_jun3=NULL;y_jun00=NULL
      y_jun3<-scan(input$fileDataset$datapath,what = "",sep = "") 
      y_jun00<-scan(input$fileDataset$datapath,what = "",sep = "\n")
      output<-list(y_jun3,y_jun00)
    })
    
    output$genTable<-renderDataTable({
      if(input$dataformat=="GCIM"){
        genRaw1_1<-as.matrix(rawGCIM()[[1]])
        genR_show<-genRaw1_1[-1,]
        colnames(genR_show)<-genRaw1_1[1,]
        as.data.frame(genR_show)
        
      }else if(input$dataformat=="QTLIciMapping"){
        as.data.frame(rawICIM()[[1]])
      }
    })
    
    output$pheTable<-renderDataTable({
      if(input$dataformat=="GCIM"){
        pheRaw1_1<-as.matrix(rawGCIM()[[2]])
        pheR_show<-pheRaw1_1[-1,]
        colnames(pheR_show)<-pheRaw1_1[1,]
        as.data.frame(pheR_show)
        
      }else if(input$dataformat=="QTLIciMapping"){
        as.data.frame(rawICIM()[[2]])
      }
    })
    output$mapTable<-renderDataTable({
      if(input$dataformat=="GCIM"){
        mapRaw1_1<-as.matrix(rawGCIM()[[3]])
        mapR_show<-mapRaw1_1[-1,]
        colnames(mapR_show)<-mapRaw1_1[1,]
        as.data.frame(mapR_show)
      }else if(input$dataformat=="QTLIciMapping"){
        as.data.frame(rawICIM()[[3]])
      }
    })
    
    output$covTable<-renderDataTable({
      if(input$dataformat=="GCIM"){
        covraw<-as.matrix(rawGCIM()[[4]])
        covc<-covraw[-1,]
        colnames(covc)<-covraw[1,]
        as.data.frame(covc)
      }else if(input$dataformat=="QTLIciMapping"){
        as.data.frame(rawICIM()[[5]])
      }
    })
    
    output$TableWIN<-renderText({
      rawWIN()[[2]]
    })
    
    
    ReadData<-reactive({
      
      geoo=NULL;pho=NULL;poss=NULL;parm=NULL;y_jun3=NULL;genRaw1=NULL;pheRaw=NULL;mapRaw11=NULL;cov_en=NULL
      if(input$dataformat=="QTLIciMapping"){
        geoo<-rawICIM()[[1]]
        pho<-rawICIM()[[2]]
        poss<-rawICIM()[[3]]
        parm<-rawICIM()[[4]]
        cov_en<-rawICIM()[[5]]
      }else if(input$dataformat=="WinQTLCart"){
        y_jun3<-rawWIN()[[1]]
      }else if(input$dataformat=="GCIM"){
        genRaw1<-rawGCIM()[[1]]
        pheRaw<-rawGCIM()[[2]]
        mapRaw11<-rawGCIM()[[3]]
        cov_en<-rawGCIM()[[4]]
      }
      result<-list(geoo=geoo,pho=pho,poss=poss,parm=parm,y_jun3=y_jun3,genRaw1=genRaw1,pheRaw=pheRaw,mapRaw11=mapRaw11,cov_en=cov_en)  
    })
    
    
    plotDH<-reactive({
      
      plotgcimFunc <- function(mxmp,galaxyy1,res11,chr_name,legend_size,mainline_size,backline_size,margin_space,axis_space,logPCoff,color1,color2,lodthred)
      {
        chr_pos <- mxmp[,1:2]
        chr_num <- length(chr_name)
        chr <- matrix(0,chr_num,1)
        pos <- matrix(0,chr_num,1)
        for(i in 1:chr_num)
        {
          temp <- numeric()
          temp <- length(which(chr_pos[,1]==i))
          if(i==1)
          {
            pos[i] <- temp
            chr[i] <- chr_pos[pos[i],2]
          }else{
            pos[i] <- pos[i-1] + temp
            chr[i] <- chr_pos[pos[i],2]
          }
        }
        
        pos_acc <- matrix(0,chr_num,1)
        for(i in 1:chr_num)
        {
          if(i==1){
            pos_acc[i] <- chr[i]
          }else{
            pos_acc[i] <- pos_acc[i-1] + chr[i]
          }
        }
        
        firFil <- res11[,1:2]
        newposadd <- as.matrix(firFil[,2])
        for(i in 1:chr_num)
        {
          temp1 <- numeric()
          temp1 <- which(firFil[,1]==i)
          if(i>1)
          {
            newposadd[temp1] <- newposadd[temp1]+pos_acc[i-1]
          }
        }
        if(is.null(galaxyy1)==FALSE){
          if(is.null(dim(galaxyy1))==TRUE){
            galaxyy1<-matrix(galaxyy1,1,3)
          }
          newres_pos <- galaxyy1[,2]
          res_sumpos <- pos_acc[galaxyy1[which(galaxyy1[,1]>1),1]-1] + galaxyy1[which(galaxyy1[,1]>1),2]
          newres_pos[which(galaxyy1[,1]>1)] <- res_sumpos
          pospic<-c(newres_pos)
          lodpic<-c(galaxyy1[,3])
          resdf <- data.frame(pospic,lodpic)
        }
        
        resp<-as.matrix(res11[,3])
        pmin<-min(resp[resp!=0])
        locsub<-which(resp==0)
        if(length(locsub)!=0){
          subvalue<-10^(1.1*log10(pmin))
          res11[locsub,3]<-subvalue
        }else{
          res11<-res11
        }
        
        negloP <- -log10(as.matrix(res11[,3]))
        if(is.null(galaxyy1)==FALSE){
          par(mar=c(2*margin_space,2*margin_space,2*margin_space,2*margin_space)+margin_space,mgp=c(3*axis_space,axis_space,0))
          plot(newposadd,negloP,type="l",col=color2,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,xlim=c(0,max(newposadd)),ylim=c(0,logPCoff*max(negloP)))
          par(new=TRUE)
          plot(pospic,lodpic,type="h",col=color1,xlab="",ylab="Logarithm of odds (LOD)",cex.axis=legend_size,cex.lab=legend_size,lwd=mainline_size,xlim=c(0,max(newposadd)),ylim=c(0,max(lodpic)))
          abline(h=lodthred)
          axis(side=4,cex.axis=legend_size)
          mtext(expression('-log'[10]*'(P)'),side=4,line=3*axis_space,cex=legend_size)
          abline(v=pos_acc,lty=2,col="gray")
        }else{
          plot(newposadd,negloP,type="l",col=color1,xlab="Genome position (cM)",ylab=expression('Expected -log'[10]*'(P)'),cex.axis=legend_size,cex.lab=legend_size,lwd=mainline_size,xlim=c(0,max(newposadd)),ylim=c(0,logPCoff*max(negloP)))
        }
      }
      
    })
    
    
    plotF2<-reactive({
      
      plotgcimFuncF2 <- function(mxmp,galaxyy1,res1a,res1d,chr_name,legend_size,mainline_size,backline_size,margin_space,axis_space,logPCoff,color1,color2,color3,lodthred)
      {
        chr_pos <- mxmp[,1:2]
        chr_num <- length(chr_name)
        chr <- matrix(0,chr_num,1)
        pos <- matrix(0,chr_num,1)
        for(i in 1:chr_num)
        {
          temp <- numeric()
          temp <- length(which(chr_pos[,1]==i))
          if(i==1)
          {
            pos[i] <- temp
            chr[i] <- chr_pos[pos[i],2]
          }else{
            pos[i] <- pos[i-1] + temp
            chr[i] <- chr_pos[pos[i],2]
          }
        }
        
        pos_acc <- matrix(0,chr_num,1)
        for(i in 1:chr_num)
        {
          if(i==1){
            pos_acc[i] <- chr[i]
          }else{
            pos_acc[i] <- pos_acc[i-1] + chr[i]
          }
        }
        firFila <- res1a[,1:2]
        newposadda <- as.matrix(firFila[,2])
        for(i in 1:chr_num)
        {
          temp1a <- numeric()
          temp1a <- which(firFila[,1]==i)
          if(i>1)
          {
            newposadda[temp1a] <- newposadda[temp1a]+pos_acc[i-1]
          }
        }
        firFild <- res1d[,1:2]
        newposaddd <- as.matrix(firFild[,2])
        for(i in 1:chr_num)
        {
          temp1d <- numeric()
          temp1d <- which(firFild[,1]==i)
          if(i>1)
          {
            newposaddd[temp1d] <- newposaddd[temp1d]+pos_acc[i-1]
          }
        }
        if(is.null(galaxyy1)==FALSE){
          if(is.null(dim(galaxyy1))==TRUE){
            galaxyy1<-matrix(galaxyy1,1,3)
          }
          newres_pos <- galaxyy1[,2]
          res_sumpos <- pos_acc[galaxyy1[which(galaxyy1[,1]>1),1]-1] + galaxyy1[which(galaxyy1[,1]>1),2]
          newres_pos[which(galaxyy1[,1]>1)] <- res_sumpos
          pospic<-c(newres_pos)
          lodpic<-c(galaxyy1[,3])
          resdf <- data.frame(pospic,lodpic)
        }
        negloPa <- as.matrix(res1a[,3])
        negloPd <- as.matrix(res1d[,3])
        if(is.null(galaxyy1)==FALSE){
          par(mar=c(2*margin_space,2*margin_space,2*margin_space,2*margin_space)+margin_space,mgp=c(3*axis_space,axis_space,0))
          plot(newposadda,negloPa,type="l",col=color3,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,logPCoff*max(negloPa)))
          par(new=TRUE)
          plot(newposaddd,negloPd,type="l",col=color2,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,logPCoff*max(negloPd)))
          par(new=TRUE)
          plot(pospic,lodpic,type="h",col=color1,xlab="",ylab="Logarithm of odds (LOD)",cex.axis=legend_size,cex.lab=legend_size,lwd=mainline_size,xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,max(lodpic)))
          abline(h=lodthred)
          axis(side=4,cex.axis=legend_size)
          mtext(expression('-log'[10]*'(P)'),side=4,line=3*axis_space,cex=legend_size)
          abline(v=pos_acc,lty=2,col="gray")
        }else{
          plot(newposadda,negloPa,type="l",col=color3,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,logPCoff*max(negloPa)))
          par(new=TRUE)
          plot(newposaddd,negloPd,type="l",col=color2,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,logPCoff*max(negloPd)))
        }
      }
      
    })
    
    
    GCIM<-eventReactive(input$run,{
      
      dir2<-paste(input$SavePath,"/",sep="")
      setwd(dir2)
      
      DoData<-function(fileFormat=NULL,Population=NULL,Model=NULL,readraw=NULL){
        pheRaw=NULL;genRaw=NULL;mapRaw1=NULL;flag=NULL;yygg1=NULL;cov_en=NULL;flagRIL=NULL;
        geoo<-readraw$geoo;pho<-readraw$pho;poss<-readraw$poss;parm<-readraw$parm;y_jun3<-readraw$y_jun3
        genRaw1<-readraw$genRaw1;pheRaw<-readraw$pheRaw;mapRaw11<-readraw$mapRaw11;cov_en<-readraw$cov_en
        
        if(fileFormat=="QTLIciMapping"){
          if(parm[4,1]==1){
            pos.be<-numeric()
            for(i in 1:10){
              pos1<-poss[which(poss[,2]==i),]
              poss1<-pos1
              positi<-as.matrix(cumsum(poss1[,3]))
              chrr<-as.matrix(poss1[,1:2])
              poss2<-cbind(chrr,positi)
              pos.be<-rbind(pos.be,poss2)
            }
          }
          if(parm[4,1]==2){
            pos.be<-poss
          }
          if(parm[5,1]==2){
            posthree<-matrix(100*pos.be[,3],,1)
            postwo<-pos.be[,1:2]
            pos<-cbind(postwo,posthree)
          }
          if(parm[5,1]==1){
            pos<-as.matrix(pos.be)
          }
          pos<-as.matrix(pos)
          geo<-geoo
          
          if(Model=="Random model"){
            flag<-1
          }else if(Model=="Fixed model"){
            flag<-0
          }  
          gen_0<-geo[,-1]
          gen_0<-as.data.frame(gen_0,stringsAsFactors = F)
          gen_0<-sapply(gen_0,as.numeric)
          gen_0<-as.matrix(gen_0)
          gen_0<-matrix(as.character(gen_0),nrow(gen_0),ncol(gen_0))
          
          if(Population=="BC1"){
            gen_1<-gsub("-1","99",gen_0)
            gen_2<-gsub("1","-1",gen_1)
            gen_11<-gsub("2","1",gen_2)
            flagRIL<-0
          }else if(Population=="BC2"){
            gen_1<-gsub("-1","99",gen_0)
            gen_11<-gsub("0","-1",gen_1)
            flagRIL<-0
          }else if(Population=="DH"){
            gen_1<-gsub("-1","99",gen_0)
            gen_2<-gsub("0","-1",gen_1)
            gen_11<-gsub("2","1",gen_2)
            flagRIL<-0
          }else if(Population=="RIL"){
            gen_1<-gsub("-1","99",gen_0)
            gen_2<-gsub("0","-1",gen_1)
            gen_11<-gsub("2","1",gen_2)
            flagRIL<-1
          }else if(Population=="F2"){
            gen_1<-gsub("12","D",gen_0)
            gen_2<-gsub("10","C",gen_1)
            gen_3<-gsub("0","B",gen_2)
            gen_4<-gsub("2","A",gen_3)
            gen_5<-gsub("-1","-",gen_4)
            gen_11<-gsub("1","H",gen_5)
            flagRIL<-0
          }
          
          phett<-t(pho)
          phe_m<-as.matrix(phett[-1,])
          phe_00<-gsub(-100,NA,phe_m)
          
          seq_indiv<-seq(1,nrow(phe_00))
          seq_indiv1<-c("genotype",seq_indiv)
          seq_indiv1<-matrix(seq_indiv1,1,)
          geo1<-cbind(geo[,1],gen_11)
          genRaw<-rbind(seq_indiv1,geo1)
          
          seq_indiv2<-c("phenotype",seq_indiv)
          seq_indiv2<-matrix(seq_indiv2,,1)
          phename<-matrix(phett[1,],1,)
          phe<-rbind(phename,phe_00)
          pheRaw<-cbind(seq_indiv2,phe)
          
          colname_mapRaw1<-c("marker","chr","pos")
          colname_mapRaw1<-matrix(colname_mapRaw1,1,)
          mapRaw1<-rbind(colname_mapRaw1,pos)
          
          if(is.null(cov_en)==FALSE){
            cov_en1<-cov_en[-1,2:ncol(cov_en)]
            covnum<-t(cov_en1)
            yygg1<-numeric()
            for(i in 1:nrow(covnum)){
              otrait_ind<-unique(covnum[i,])
              cov_col<-length(otrait_ind)-1
              col_each<-numeric()
              for(j in 1:length(covnum[i,])){
                if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                  cov_0<-matrix(-1,1,cov_col)
                }else{
                  cov_0<-matrix(0,1,cov_col)
                  covnum_loc<-which(otrait_ind[]==covnum[i,j])
                  cov_0[1,covnum_loc]<-1
                }
                col_each<-rbind(col_each,cov_0)
                
              }
              yygg1<-cbind(yygg1,col_each)
            }
          }else{
            yygg1<-NULL 
          }
          
        }else if(fileFormat=="WinQTLCart"){
          if(Model=="Random model"){
            flag<-1
          }else if(Model=="Fixed model"){
            flag<-0
          }  
          start_dex<-grep("-start",y_jun3,fixed = TRUE)
          stop_dex<-grep("-stop",y_jun3,fixed = TRUE)
          
          chr_dex<-grep("-Chromosome",y_jun3,fixed = TRUE)
          
          chrdata<-c(y_jun3[chr_dex[1]:(stop_dex[1]-1)],"-Chromosome",y_jun3[stop_dex[1]])
          chrdata_dex<-grep("-Chromosome",chrdata,fixed = TRUE)
          chrdata_dexlen<-length(chrdata_dex)
          
          chr_num<-numeric()
          chrname_num<-numeric()
          chr_numfirst<-numeric()
          markername0<-numeric()
          chr_pos<-numeric()
          chrRaw_name<-numeric()
          chr_Rawnumfirst<-numeric()
          for(i in 1:(chrdata_dexlen-1)){
            chr_each<-chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]
            marker_name<-numeric()
            marker_pos<-numeric()
            
            for(j in 0:(trunc(length(chr_each)/2)-1) ){
              marker_name0<-chr_each[2*j+1]
              marker_name<-cbind(marker_name,marker_name0)
              marker_pos0<-suppressWarnings(as.numeric(chr_each[2*(j+1)]))
              marker_pos<-cbind(marker_pos,marker_pos0)
              if((y_jun3[5]=="interval")&&(y_jun3[7]=="1")&&(y_jun3[9]=="r")){
                marker_posm<-100*((-0.5)*log(1-2*marker_pos))
                
                markerlen<-length(marker_posm)
                marker_pos1<-marker_posm[1:(markerlen-1)]
                marker_pos2<-c(0,marker_pos1)
                marker_possum<-cumsum(marker_pos2)
              }
              if((y_jun3[5]=="interval")&&(y_jun3[7]=="2")&&(y_jun3[9]=="r")){
                marker_posm<-100*(0.25*log((1+2*marker_pos)/(1-2*marker_pos)))
                
                markerlen<-length(marker_posm)
                marker_pos1<-marker_posm[1:(markerlen-1)]
                marker_pos2<-c(0,marker_pos1)
                marker_possum<-cumsum(marker_pos2)
                
              }
              if((y_jun3[5]=="interval")&&(y_jun3[9]=="M")){
                marker_posm<-100*marker_pos
                
                markerlen<-length(marker_posm)
                marker_pos1<-marker_posm[1:(markerlen-1)]
                marker_pos2<-c(0,marker_pos1)
                marker_possum<-cumsum(marker_pos2)
              }
              if((y_jun3[5]=="interval")&&(y_jun3[9]=="cM")){
                marker_posm<-marker_pos
                
                markerlen<-length(marker_posm)
                marker_pos1<-marker_posm[1:(markerlen-1)]
                marker_pos2<-c(0,marker_pos1)
                marker_possum<-cumsum(marker_pos2)
              }
              if((y_jun3[5]=="position")&&(y_jun3[9]=="M")){
                marker_possum<-100*marker_pos
              }
              if((y_jun3[5]=="position")&&(y_jun3[9]=="cM")){
                marker_possum<-marker_pos
              }
            }
            markername0<-cbind(markername0,marker_name)
            markername<-matrix(markername0,,1)
            marker_possum0<-matrix(marker_possum,,1)
            chr_pos<-rbind(chr_pos,marker_possum0)
            chr_a<-suppressWarnings(as.numeric(chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]))
            chr_data<-na.omit(chr_a)
            chr_datalen<-length(chr_data)
            chr_num<-rbind(chr_num,chr_datalen)
            
            chrRawname<-chrdata[(chrdata_dex[i]+1)]
            chrname<-str_extract_all(chrRawname,"[0-9]+")
            chrRawname<-matrix(chrRawname,,1)
            chrRaw_name<-rbind(chrRaw_name,chrRawname)
            
            chrname_num<-rbind(chrname_num,chrname)
            chr_numxx<-rep(as.numeric(chrname_num[i]),chr_num[i])
            chr_numfirst<-rbind(chr_numfirst,matrix(chr_numxx,,1))
            
            chr_Rawnumxx<-rep(chrRaw_name[i],chr_num[i])
            chr_Rawnumfirst<-rbind(chr_Rawnumfirst,matrix(chr_Rawnumxx,,1))
          }
          
          chr_leng<-length(chr_pos)
          chr_numtwo<-cbind(chr_numfirst,chr_pos)
          marker_dex<-grep("markers",y_jun3,fixed = TRUE)
          marker_snp<-y_jun3[start_dex[2]:stop_dex[2]]
          
          marker_snpnum<-marker_snp
          snpa<-suppressWarnings(as.numeric(marker_snpnum))
          
          
          snpdata<-na.omit(snpa)
          indi_num<-length(snpdata)/chr_leng
          snp_data<-numeric()
          
          if_indi<-y_jun3[marker_dex[1]-1]
          if(if_indi=="individuals"){
            
            for(i in 0:(indi_num-1)){
              snp_eve<-matrix(snpdata[(chr_leng*i+1):(chr_leng*(i+1))],1,)
              snp_data<-rbind(snp_data,snp_eve)
            }
            snp_data1<-t(snp_data)
          }else{
            for(i in 0:(chr_leng-1)){
              snp_eve<-matrix(snpdata[(indi_num*i+1):(indi_num*(i+1))],1,)
              snp_data<-rbind(snp_data,snp_eve)
            }
            snp_data1<-snp_data
          }
          
          trait_total<-y_jun3[start_dex[3]:stop_dex[3]]
          
          for(i in 1:length(trait_total)){
            if(trait_total[i]=="."){
              trait_total[i]<-"0"
            }
          }
          trait_dex<-grep("traits",trait_total)
          traita<-suppressWarnings(as.numeric(trait_total))
          traitdata<-na.omit(traita)
          trait_num<-length(traitdata)/indi_num
          trait_data<-numeric()
          
          iftrait_indi<-trait_total[trait_dex[1]-1]
          if(iftrait_indi=="individuals"){
            
            for(i in 0:(indi_num-1)){
              trait_bbb<-traitdata[(trait_num*i+1):(trait_num*(i+1))]
              for(j in 1:length( trait_bbb)){
                if(trait_bbb[j]==0){
                  trait_bbb[j]<-NA
                }
              }
              trait_eve<-matrix(trait_bbb,1,)
              trait_data<-rbind(trait_data,trait_eve)
            }
            
          }else{
            
            for(i in 0:(trait_num-1)){
              trait_aaa<-traitdata[(indi_num*i+1):(indi_num*(i+1))]
              for(j in 1:length(trait_aaa)){
                if(trait_aaa[j]==0){
                  trait_aaa[j]<-NA
                }
              }
              trait_eve<-matrix(trait_aaa,,1)
              trait_data<-cbind(trait_data,trait_eve)
            }
            
          }
          if(length(start_dex)==3){yygg1<-NULL}
          if(length(start_dex)==4){
            if(y_jun3[start_dex[4]+1]=="otraits"||y_jun3[start_dex[4]+1]=="individuals"){
              cov_total<-y_jun3[start_dex[4]:stop_dex[4]]
              cov_dex<-grep("otraits",cov_total)
              cov_only<-y_jun3[(start_dex[4]+2):(stop_dex[4]-1)]
              bycross_dex<-grep("#bycross",y_jun3,fixed = TRUE)
              otrait_indi<-as.numeric(y_jun3[bycross_dex+8])
              if(y_jun3[start_dex[4]+1]=="otraits"){
                
                covnumonly<-numeric()
                for( i in 0:(otrait_indi-1)){
                  cov_each<-matrix(cov_only[(i*(indi_num+1)+1):((indi_num+1)*(i+1))],1,)
                  covnumonly<-rbind(covnumonly,cov_each)
                }
                covnum<-covnumonly[,-1]
                yygg1<-numeric()
                for(i in 1:nrow(covnum)){
                  
                  otrait_ind<-unique(covnum[i,])
                  cov_col<-length(otrait_ind)-1
                  
                  col_each<-numeric()
                  for(j in 1:length(covnum[i,])){
                    
                    if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                      cov_0<-matrix(-1,1,cov_col)
                      
                    }else{
                      cov_0<-matrix(0,1,cov_col)
                      covnum_loc<-which(otrait_ind[]==covnum[i,j])
                      cov_0[1,covnum_loc]<-1
                    }
                    col_each<-rbind(col_each,cov_0)
                    
                  }
                  yygg1<-cbind(yygg1,col_each)
                  
                }
                
              }
              if(y_jun3[start_dex[4]+1]=="individuals"){
                covdata<-cov_only[(2+otrait_indi):length(cov_only)]
                covnum<-numeric()
                otrait_name<-matrix(cov_only[2:2+(otrait_indi-1)],otrait_indi,1)
                for(m in 1:otrait_indi){
                  
                  coveach<-numeric()
                  for(n in 0:(indi_num-1)){
                    cov_each<-matrix(covdata[m+otrait_indi*n],1,1)
                    
                    coveach<-cbind(coveach,cov_each)
                  }
                  covnum<-rbind(covnum,coveach)
                }
                covnumonly<-cbind(otrait_name,covnum)
                yygg1<-numeric()
                for(i in 1:nrow(covnum)){
                  
                  otrait_ind<-unique(covnum[i,])
                  cov_col<-length(otrait_ind)-1
                  
                  col_each<-numeric()
                  for(j in 1:length(covnum[i,])){
                    
                    if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                      cov_0<-matrix(-1,1,cov_col)
                      
                    }else{
                      cov_0<-matrix(0,1,cov_col)
                      covnum_loc<-which(otrait_ind[]==covnum[i,j])
                      cov_0[1,covnum_loc]<-1
                    }
                    col_each<-rbind(col_each,cov_0)
                    
                  }
                  yygg1<-cbind(yygg1,col_each)
                  
                }
              }
              
            }
          }
          
          seq_indi<-seq(1,nrow(trait_data))
          seq_indi1<-c("genotype",seq_indi)
          seq_indi1<-matrix(seq_indi1,1,)
          snp1<-cbind(markername,snp_data1)
          seq_indi2<-c("phenotype",seq_indi)
          seq_indi2<-matrix(seq_indi2,,1)
          num_trait<-ncol(trait_data)
          seq_trait<-seq(1,num_trait)
          seq_trait<-matrix(seq_trait,1,)
          trait_data00<-rbind(seq_trait,trait_data)
          
          colnames_mapname<-c("marker","chr","pos")
          colnames_mapname<-matrix(colnames_mapname,1,)
          mapRaw1<-cbind(markername,chr_Rawnumfirst,chr_pos)
          
          mapRaw1<-rbind(colnames_mapname,mapRaw1)
          genRawqq<-rbind(seq_indi1,snp1)
          genRawq<-genRawqq[-1,]
          
          if(Population=="BC1"){
            
            marker_snp1<-gsub("1","-1",genRawq)
            marker_snp2<-gsub("2","1",marker_snp1)
            genRawh<-gsub("*","99",marker_snp2,fixed = TRUE)
            genRaw<-rbind(genRawqq[1,],genRawh)
            flagRIL<-0
          }else if(Population=="BC2"){
            
            marker_snp1<-gsub("0","-1",genRawq)
            #marker_snp2<-gsub("1","1",marker_snp1)
            genRawh<-gsub("*","99",marker_snp1,fixed = TRUE) 
            genRaw<-rbind(genRawqq[1,],genRawh)
            flagRIL<-0
          }else if(Population=="DH"){
            
            marker_snp1<-gsub("0","-1",genRawq)
            marker_snp2<-gsub("2","1",marker_snp1)
            genRawh<-gsub("*","99",marker_snp2,fixed = TRUE)  
            genRaw<-rbind(genRawqq[1,],genRawh)
            flagRIL<-0
          }else if(Population=="RIL"){
            marker_snp10<-gsub("1","99",genRawq)
            marker_snp1<-gsub("0","-1",marker_snp10)
            marker_snp2<-gsub("2","1",marker_snp1)
            genRawh<-gsub("*","99",marker_snp2,fixed = TRUE)
            genRaw<-rbind(genRawqq[1,],genRawh)
            flagRIL<-1
          }else if(Population=="F2"){
            marker_snp1<-gsub("12","D",genRawq)
            marker_snp2<-gsub("10","C",marker_snp1)
            marker_snp3<-gsub("0","B",marker_snp2)
            marker_snp4<-gsub("2","A",marker_snp3)
            marker_snp5<-gsub("-1","-",marker_snp4)
            genRawh<-gsub("1","H",marker_snp5,fixed = TRUE)  
            genRaw<-rbind(genRawqq[1,],genRawh)
            flagRIL<-0
          }
          pheRaw<-cbind(seq_indi2,trait_data00)
          
          seq_indi3<-c("covariate",seq_indi)
          seq_indi3<-matrix(seq_indi3,,1)
          if(is.null(yygg1)==FALSE){
            cov_en<-cbind(seq_indi3,t(covnumonly))
          }else{
            cov_en<-NULL
          }
        }else if(fileFormat=="GCIM"){
          genRaw1qq<-as.matrix(genRaw1)
          genRaw1<-genRaw1qq[-1,-1]
          
          pheRaw<-as.matrix(pheRaw)
          mapRaw1<-as.matrix(mapRaw11)
          if(Population=="BC1"){
            genRaw_<-gsub("-","99",genRaw1)
            genRaw_Aa<-gsub("H","-1",genRaw_)
            genRaw_AA<-gsub("A","1",genRaw_Aa)
            
            genRaw<-rbind(genRaw1qq[1,],cbind(genRaw1qq[-1,1],genRaw_AA))
            flagRIL<-0
          }else if(Population=="BC2"){
            genRaw_<-gsub("-","99",genRaw1)
            genRaw_aa<-gsub("B","-1",genRaw_)
            genRaw_Aa<-gsub("H","1",genRaw_aa)
            
            genRaw<-rbind(genRaw1qq[1,],cbind(genRaw1qq[-1,1],genRaw_Aa))
            flagRIL<-0
          }else if(Population=="DH"){
            genRaw_<-gsub("-","99",genRaw1)
            genRaw_aa<-gsub("B","-1",genRaw_)
            genRaw_AA<-gsub("A","1",genRaw_aa)
            
            genRaw<-rbind(genRaw1qq[1,],cbind(genRaw1qq[-1,1],genRaw_AA))
            flagRIL<-0 
          }else if(Population=="RIL"){
            genRaw_<-gsub("-","99",genRaw1)
            genRaw_aa<-gsub("B","-1",genRaw_)
            genRaw_AA<-gsub("A","1",genRaw_aa)
            genRaw<-rbind(genRaw1qq[1,],cbind(genRaw1qq[-1,1],genRaw_AA))
            flagRIL<-1
          }else if(Population=="F2"){
            genRaw<-genRaw1qq  
            flagRIL<-0
          }
          
          if(Model=="Random model"){
            flag<-1
          }else if(Model=="Fixed model"){
            flag<-0
          }  
          
          if(is.null(cov_en)==FALSE){
            cov_en1<-cov_en[-1,2:ncol(cov_en)]
            covnum<-t(cov_en1)
            yygg1<-numeric()
            for(i in 1:nrow(covnum)){
              otrait_ind<-unique(covnum[i,])
              cov_col<-length(otrait_ind)-1
              col_each<-numeric()
              for(j in 1:length(covnum[i,])){
                if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                  cov_0<-matrix(-1,1,cov_col)
                }else{
                  cov_0<-matrix(0,1,cov_col)
                  covnum_loc<-which(otrait_ind[]==covnum[i,j])
                  cov_0[1,covnum_loc]<-1
                }
                col_each<-rbind(col_each,cov_0)
                
              }
              yygg1<-cbind(yygg1,col_each)
            }
          }else{
            yygg1<-NULL 
          }
        }
        result<-list(pheRaw=pheRaw,genRaw=genRaw,mapRaw1=mapRaw1,flag=flag,flagRIL=flagRIL,yygg1=yygg1,cov_en=cov_en)
        return(result)
      }
    
      
      dir<-input$SavePath
      fileFormat<-input$dataformat;Model<-input$model;Population<-input$Pop
      WalkSpeed<-as.numeric(input$Walk);CriLOD<-as.numeric(input$Crilod);Likelihood<-input$likelihood;flagrqtl<-input$Rqtl
      DrawPlot<-input$drawplot;PlotFormat<-input$Plotformat;Resolution<-input$resolution
      trait<-as.numeric(input$trait);
      
      readraw<-ReadData()
      DoResult<-DoData(fileFormat,Population,Model,readraw)
      
      print("Running in progress, please be patient...")
      
      pheRaw<-DoResult$pheRaw;genRaw<-DoResult$genRaw;mapRaw1<-DoResult$mapRaw1
      flag<-DoResult$flag;flagRIL<-DoResult$flagRIL;yygg1<-DoResult$yygg1;cov_en<-DoResult$cov_en
      
      W1re<-NULL;WEN1re<-NULL
      
      if(Resolution=="General"){
        widqqvalue<-1500
        heightqqvalue<-600
        pointsizeqqvalue<-12
        resppi<-72 
      }else if(Resolution=="High"){
        widqqvalue<-10000
        heightqqvalue<-4000
        pointsizeqqvalue<-12
        resppi<-300 
      }
      
      gcimFunc<-plotDH()
      gcimFuncF2<-plotF2()
      
      
      if(length(grep(":",input$trait))!=0){
        scope<-as.numeric(unlist(strsplit(input$trait,":")))
        trait<-as.numeric(scope[1]:scope[2])
      }else if(length(grep(",",input$trait))!=0){
        scope<-as.numeric(unlist(strsplit(input$trait,",")))
        trait<-c(scope)
      }else if(length(grep(":",input$trait))==0&&length(grep(",",input$trait))==0){
        trait<-as.numeric(input$trait)
      }
      
      
      if(Population=="F2"){
        
        withProgress(message = 'Running in progress', value = 0, {
        
        WEN1re<-WenF(pheRaw,genRaw,mapRaw1,yygg1,cov_en,WalkSpeed,CriLOD,dir)
        
        for(NUM in trait){
            rewen<-NULL;mxmp=NULL;galaxyy1<-NULL;res1a=NULL;res1d=NULL;chr_name=NULL 
            TRY1<-try({
              outWEN<-WenS(flag,CriLOD,NUM,pheRaw,Likelihood,flagrqtl,WEN1re$yygg,WEN1re$mx,WEN1re$phe,WEN1re$chr_name,
                           WEN1re$v.map,WEN1re$gen.raw,WEN1re$a.gen.orig,WEN1re$d.gen.orig,WEN1re$n,WEN1re$names.insert2,WEN1re$X.ad.tran.data,WEN1re$X.ad.t4,dir)
              rewen<-outWEN$result
              mxmp<-outWEN$mxmp;res1a<-outWEN$res1a;res1d<-outWEN$res1d;chr_name<-outWEN$chr_name;galaxyy1<-outWEN$galaxyy1
              
            },silent=FALSE)    
            
            if ('try-error' %in% class(TRY1)|| !('try-error' %in% class(TRY1))){  
              TRY2<-try({ 
                write.table(rewen,paste(NUM,"_GCIM result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)
                
                if(is.null(mxmp)==FALSE){colnames(mxmp)<-c("chr","pos")}
                if(is.null(res1a)==FALSE){colnames(res1a)<-c("chr","pos","-log10(p-value)")}
                if(is.null(res1d)==FALSE){colnames(res1d)<-c("chr","pos","-log10(p-value)")}
                if(is.null(chr_name)==FALSE){chr_name<-as.matrix(chr_name);colnames(chr_name)<-c("chr")}
                if(is.null(galaxyy1)==FALSE){colnames(galaxyy1)<-c("chr","pos","LOD")}
                
                plotresult<-list(mxmp,res1a,res1d,chr_name,galaxyy1)
                
                wb <- createWorkbook("Fred")
                addWorksheet(wb, "sheet1")
                addWorksheet(wb, "sheet2")
                addWorksheet(wb, "sheet3")
                addWorksheet(wb, "sheet4")
                addWorksheet(wb, "sheet5")
                
                writeData(wb, sheet = "sheet1", plotresult[[1]])
                writeData(wb, sheet = "sheet2", plotresult[[2]])
                writeData(wb, sheet = "sheet3", plotresult[[3]])
                writeData(wb, sheet = "sheet4", plotresult[[4]])
                writeData(wb, sheet = "sheet5", plotresult[[5]])
                
                saveWorkbook(wb,paste(NUM,"_resultforplot.xlsx",sep=""), overwrite = TRUE)
                
                if(DrawPlot==TRUE){
                  if(PlotFormat=="*.png")
                  {
                    png(paste(NUM,"_resF2.png"), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
                  }else if(PlotFormat=="*.tiff"){
                    tiff(paste(NUM,"_resF2.tiff"), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
                  }else if(PlotFormat=="*.jpeg"){
                    jpeg(paste(NUM,"_resF2.jpeg"), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
                  }else if(PlotFormat=="*.pdf"){
                    pdf(paste(NUM,"_resF2.pdf"), width=16)
                  }
                  gcimFuncF2(mxmp,galaxyy1,res1a,res1d,chr_name,1.0,1.0,0.5,1.5,1.0,1.5,"red","gray50","green",CriLOD)
                  dev.off()
                }
              },silent=FALSE)  
            }
            incProgress(1/(max(trait)-min(trait)+1), detail = paste("Doing part", NUM))
            Sys.sleep(0.1)
          }
        })
        
      }else{
        
        withProgress(message = 'Running in progress', value = 0, {
        
        W1re<-WangF(pheRaw,genRaw,mapRaw1,yygg1,flagRIL,cov_en,Population,WalkSpeed,CriLOD)
        
        for(NUM in trait){
            rew<-NULL;mxmp=NULL;galaxyy1<-NULL;res11=NULL;chr_name=NULL
            TRY1<-try({
              outW<-WangS(flag,CriLOD,NUM,pheRaw,W1re$chrRaw_name,W1re$yygg,W1re$mx,W1re$phe,W1re$chr_name,W1re$gen,W1re$mapname,CLO=NULL)
              rew<-outW$result
              mxmp<-outW$mxmp;res11<-outW$res11;chr_name<-outW$chr_name;galaxyy1<-outW$galaxyy1
            },silent=FALSE)  
            
            if ('try-error' %in% class(TRY1)|| !('try-error' %in% class(TRY1))){   
              TRY2<-try({ 
                write.table(rew,paste(NUM,"_GCIM result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)
                if(is.null(mxmp)==FALSE){colnames(mxmp)<-c("chr","pos")}
                if(is.null(res11)==FALSE){colnames(res11)<-c("chr","pos","p-value")}
                if(is.null(chr_name)==FALSE){chr_name<-as.matrix(chr_name);colnames(chr_name)<-c("chr")}
                if(is.null(galaxyy1)==FALSE){colnames(galaxyy1)<-c("chr","pos","LOD")}
                plotresult<-list(mxmp,res11,chr_name,galaxyy1)
                wb <- createWorkbook("Fred")
                addWorksheet(wb, "sheet1")
                addWorksheet(wb, "sheet2")
                addWorksheet(wb, "sheet3")
                addWorksheet(wb, "sheet4")
                writeData(wb, sheet = "sheet1", plotresult[[1]])
                writeData(wb, sheet = "sheet2", plotresult[[2]])
                writeData(wb, sheet = "sheet3", plotresult[[3]])
                writeData(wb, sheet = "sheet4", plotresult[[4]])
                saveWorkbook(wb,paste(NUM,"_resultforplot.xlsx",sep=""), overwrite = TRUE)
                if(DrawPlot==TRUE){
                  if(PlotFormat=="*.png")
                  {
                    png(paste(NUM,"_res.png"), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
                  }else if(PlotFormat=="*.tiff"){
                    tiff(paste(NUM,"_res.tiff"), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
                  }else if(PlotFormat=="*.jpeg"){
                    jpeg(paste(NUM,"_res.jpeg"), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
                  }else if(PlotFormat=="*.pdf"){
                    pdf(paste(NUM,"_res.pdf"), width=16)
                  }
                  
                  gcimFunc(mxmp,galaxyy1,res11,chr_name,1.0,1.0,0.5,1.5,1.0,1.5,"red","gray50",CriLOD)
                  dev.off()
                }
                
              },silent=FALSE)  
            }  
            incProgress(1/(max(trait)-min(trait)+1), detail = paste("Doing part", NUM))
            Sys.sleep(0.1)
          }
          
        }) 
        
      }
      
    })
    
    
    
    output$result<-renderUI(GCIM())
    
    
    output$plotGCIM<-renderPlot({
      gcimFunc<-plotDH()
      gcimFuncF2<-plotF2()
      
      plotle<-as.numeric(input$leGCIM);plotma<-as.numeric(input$maGCIM);plotba<-as.numeric(input$baGCIM)
      plotmar<-as.numeric(input$marGCIM);plotaxis<-as.numeric(input$axisGCIM);plotlog<-as.numeric(input$logGCIM);
      plotco1<-input$co1GCIM;plotlod<-as.numeric(input$lodGCIM);plotco2<-input$co2GCIM;plotco3<-input$co3GCIM;plotforma<-input$plformatGCIM
      
      req(input$fileplotGCIM)
      
      if(input$plotpopGCIM=="F2"){
        mxmp<-read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet1",colNames = TRUE)
        res1a<-read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet2",colNames = TRUE)
        res1d<-read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet3",colNames = TRUE)
        chr_name<-c(read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet4",colNames = TRUE))
        galaxyy1<-read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet5",colNames = TRUE)
        
        mxmp<-sapply(mxmp,as.numeric)
        res1a<-sapply(res1a,as.numeric)
        res1d<-sapply(res1d,as.numeric)
        chr_name<-sapply(chr_name,as.numeric)
        galaxyy1<-sapply(galaxyy1,as.numeric)
        
        gcimFuncF2(mxmp,galaxyy1,res1a,res1d,chr_name,plotle,plotma,plotba,plotmar,plotaxis,plotlog,plotco1,plotco2,plotco3,plotlod)
      }else{
        mxmp<-read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet1",colNames = TRUE)
        res11<-read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet2",colNames = TRUE)
        chr_name<-c(read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet3",colNames = TRUE))
        galaxyy1<-read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet4",colNames = TRUE)
        
        mxmp<-sapply(mxmp,as.numeric)
        res11<-sapply(res11,as.numeric)
        chr_name<-sapply(chr_name,as.numeric)
        galaxyy1<-sapply(galaxyy1,as.numeric)
        
        gcimFunc(mxmp,galaxyy1,res11,chr_name,plotle,plotma,plotba,plotmar,plotaxis,plotlog,plotco1,plotco2,plotlod)
        
      }
      
      
    })
    
    
    output$downloadplotGCIM<- downloadHandler(
      
      filename = function() {
        paste("LOD", sep = ".", switch(
          input$Plotformat2, "*.png"=".png", "*.tiff"=".tiff", "*.jpeg"=".jpeg","*.pdf"=".pdf"
        ))
      },
      content = function(file) {
        
        if(input$Plotformat2=="*.png"){
          png(file,width=as.numeric(input$widthGCIM), height=as.numeric(input$heightGCIM), units= "px", pointsize =as.numeric(input$pointGCIM),res=as.numeric(input$ppiGCIM))
        }else if(input$Plotformat2=="*.tiff"){
          tiff(file,width=as.numeric(input$widthGCIM), height=as.numeric(input$heightGCIM), units= "px", pointsize =as.numeric(input$pointGCIM),res=as.numeric(input$ppiGCIM))
        }else if(input$Plotformat2=="*.jpeg"){
          jpeg(file,width=as.numeric(input$widthGCIM), height=as.numeric(input$heightGCIM), units= "px", pointsize =as.numeric(input$pointGCIM),res=as.numeric(input$ppiGCIM))
        }else if(input$Plotformat2=="*.pdf"){
          pdf(file,width=16)
        }
        
        plotle<-as.numeric(input$leGCIM);plotma<-as.numeric(input$maGCIM);plotba<-as.numeric(input$baGCIM)
        plotmar<-as.numeric(input$marGCIM);plotaxis<-as.numeric(input$axisGCIM);plotlog<-as.numeric(input$logGCIM);
        plotco1<-input$co1GCIM;plotlod<-as.numeric(input$lodGCIM);plotco2<-input$co2GCIM;plotco3<-input$co3GCIM;plotforma<-input$plformatGCIM
        
        gcimFunc<-plotDH()
        gcimFuncF2<-plotF2()
        
        if(input$plotpopGCIM=="F2"){
          mxmp<-read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet1",colNames = TRUE)
          res1a<-read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet2",colNames = TRUE)
          res1d<-read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet3",colNames = TRUE)
          chr_name<-c(read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet4",colNames = TRUE))
          galaxyy1<-read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet5",colNames = TRUE)
          
          mxmp<-sapply(mxmp,as.numeric)
          res1a<-sapply(res1a,as.numeric)
          res1d<-sapply(res1d,as.numeric)
          chr_name<-sapply(chr_name,as.numeric)
          galaxyy1<-sapply(galaxyy1,as.numeric)
          
          gcimFuncF2(mxmp,galaxyy1,res1a,res1d,chr_name,plotle,plotma,plotba,plotmar,plotaxis,plotlog,plotco1,plotco2,plotco3,plotlod)
        }else{
          mxmp<-read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet1",colNames = TRUE)
          res11<-read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet2",colNames = TRUE)
          chr_name<-c(read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet3",colNames = TRUE))
          galaxyy1<-read.xlsx(input$fileplotGCIM$datapath,sheet = "sheet4",colNames = TRUE)
          
          mxmp<-sapply(mxmp,as.numeric)
          res11<-sapply(res11,as.numeric)
          chr_name<-sapply(chr_name,as.numeric)
          galaxyy1<-sapply(galaxyy1,as.numeric)
          gcimFunc(mxmp,galaxyy1,res11,chr_name,plotle,plotma,plotba,plotmar,plotaxis,plotlog,plotco1,plotco2,plotlod)
        }
        dev.off()
      })
  })
  
  
  ui <-tagList(
    navbarPage(
      "",id = "tabs",
      tabPanel(strong("QTL.gCIMapping.GUI"),
               
               h2("QTL.gCIMapping.GUI (QTL genome-wide Composite Interval Mapping with Graphical User Interface)",align="center"),
               column(3,
                      br(),
                      h4(strong("Coding criteria")),
                      offset=3
               ),
               column(4,
                      br(),
                      h4(strong("Dataset example")),
                      offset=1
               ),
               column(4,     
                      tableOutput("codeexample"),
                      offset = 2
               ),
               column(4,      
                      tableOutput("dataexample")
               ),
               column(12, 
                      br(),
                      h4(strong("Reference")),
                      offset = 2
               ),
               
               
               
               column(8,
                      
                      h4("1. Wang Shi-Bo, Wen Yang-Jun, Ren Wen-Long, Ni Yuan-Li, Zhang Jin, Feng Jian-Ying, Zhang Yuan-Ming*.  
                         Mapping small-effect and linked quantitative trait loci for complex traits in backcross or DH populations via
                         a multi-locus GWAS methodology. Scientific Reports 2016,6:29951."),
                      h4("2. Wen Yang-Jun, Zhang Ya-Wen, Zhang Jin, Feng Jian-Ying, Jim M. Dunwell, Zhang Yuan-Ming*. An
                         efficient multi-locus mixed model framework for the detection of small and linked QTLs 
                         in F2.Briefings in Bioinformatics doi:10.1093/bib/bby058."),
                      h4("3. Zhang Ya-Wen, Jim M. Dunwell, Zhang Yuan-Ming*. QTL.gCIMapping.GUI v2.0: An R software for detecting small-effect and linked QTLs 
                         for quantitative traits in bi-parental segregation populations. Computational and Structural Biotechnology Journal, in revision."),
                      
                      br(),
                      h4("Authors: Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, Zhang Yuan-Ming"),
                      h4("Maintainer: Zhang Yuan-Ming (soyzhang at mail.hzau.edu.cn)"), 
                      h4("QTL.gCIMapping version 2.0, Realeased November 2019"),
                      offset = 2)
                      ),
      
      
      tabPanel(strong("Start"),
               
               titlePanel("QTL.gCIMapping.GUI"),
               sidebarLayout(
                 sidebarPanel(
                   radioButtons("dataformat", "Please select data format", choices = c("GCIM","WinQTLCart","QTLIciMapping"),selected="GCIM"),
                   fileInput("fileDataset", "Input dataset",multiple = TRUE),
                   conditionalPanel("input.dataformat == 'QTLIciMapping'",
                                    fileInput("fileCov", "Input covariate file",multiple = TRUE)
                   ),
                   selectInput("fileType","Show dataset:",choices=c("Genotype","Phenotype","Linkage map","Covariate")),
                   br(),
                   br(),
                   br(),
                   actionButton("parSe", label = "Parameter Settings",width=250,style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                   br(),
                   br(),
                   br(),
                   actionButton("drawPl", label = "Figure",width = 250,style="color: #fff; background-color: #337ab7; border-color: #2e6da4"), 
                   br(),
                   br(),
                   br(),
                   actionButton("manl", label = "User manual",width=250,style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                   uiOutput("manll")
                 ),
                 mainPanel(
                   tabsetPanel(id="inTabset",
                               tabPanel("Dataset",value = "DA",
                                        conditionalPanel("input.fileType == 'Genotype'",
                                                         dataTableOutput("genTable")),
                                        conditionalPanel( "input.fileType == 'Phenotype'",
                                                          dataTableOutput("pheTable")),
                                        conditionalPanel( "input.fileType == 'Linkage map'",
                                                         dataTableOutput("mapTable")),
                                        conditionalPanel( "input.fileType == 'Covariate'",
                                                         dataTableOutput("covTable")),
                                        conditionalPanel( "input.dataformat == 'WinQTLCart'",
                                                          verbatimTextOutput("TableWIN"))
                               ),
                               tabPanel("Parameter Settings",value = "PA",
                                        fluidRow(  
                                          br(), 
                                          br(),
                                          column(6,  
                                                 radioButtons("Pop", "Please select Population", choices = c("DH","RIL","BC1","BC2","F2"),selected="DH"),
                                                 radioButtons("model", "Please select Model", choices = c("Random model","Fixed model")),
                                                 conditionalPanel("input.Pop == 'F2'",
                                                 radioButtons("likelihood","Likelihood function (only for F2):",choices=c("REML","ML")),
                                                 radioButtons("Rqtl","Completing CIM in the neighborhood (only for F2)",choices=c("TRUE","FALSE"),selected="FALSE")
                                                ),
                                                 textInput("Walk", "Walk Speed for Genome-wide Scanning (cM):",value="1"),
                                                 textInput("Crilod", "Critical LOD score",value="2.5")
                                           ),
                                          
                                          column(6, 
                                                 textInput("trait", "Traits analyzed", value="1"),
                                                 textInput("SavePath", "Save path",value = "C:/Users/Administrator/Desktop"),
                                                 radioButtons("drawplot","Draw plot or not",choices=c("TRUE","FALSE")),
                                                 conditionalPanel("input.drawplot == 'TRUE'",
                                                 radioButtons("resolution","Resolution of plot",choices=c("General","High")),
                                                 radioButtons("Plotformat","Plot format",choices=c("*.png","*.tiff","*.jpeg","*.pdf"))
                                                 )
                                               ),
                                          
                                          
                                         column(12, 
                                                 br(),
                                                 br(),
                                                 actionButton("run", label = "Run",width=280, icon("paper-plane"), 
                                                              style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                          ),
                                          column(12,
                                                 br(),
                                                 br(),
                                                 uiOutput("result")
                                          )
                                        )
                               ),
                               
                               tabPanel("Figure",value = "FIG",
                                        h3("Genome-wide composite interval mapping (GCIM) figure"),
                                        fluidRow(
                                          column(12,
                                                 radioButtons("Mose", "", c("Parameter Settings", "Draw plot"),inline = TRUE)
                                          ),
                                          
                                          conditionalPanel("input.Mose == 'Parameter Settings'",
                                                           column(6,
                                                                  br(),
                                                                  radioButtons("Resolution2", "Select resolution of plot", c("General resolution", "High resolution", "Set by yourself")),
                                                                  textInput("widthGCIM","Width (px):",value="1500"),
                                                                  textInput("heightGCIM","Height (px):",value="600"),
                                                                  textInput("pointGCIM","Word resolution (1/72 inch, ppi):",value="12"),
                                                                  textInput("ppiGCIM","Figure resolution (ppi):",value="72"),
                                                                  textInput("lodGCIM","Critical LOD score:",value="2.5"),
                                                                  textInput("maGCIM","LOD line size:",value="1.0"),
                                                                  textInput("baGCIM","Size for -log10(P) curve:",value="0.5")
                                                                 ),
                                                          column(6,
                                                                  br(),
                                                                  textInput("leGCIM","Legend and tick marks:",value="1.0"),
                                                                  textInput("marGCIM","Margin space:",value="1.5"),
                                                                  textInput("axisGCIM","Space between tick marks and axis:",value="1.0"),
                                                                  textInput("logGCIM","Times for max{-log10(P)}:",value="1.5"),
                                                                  selectInput("co1GCIM","LOD line color:",choices=c("red","black","blue","yellow","green","pink","purple","gray50","brown")),
                                                                  selectInput("co2GCIM","-log10(P) curve color1:",choices=c("gray50","black","blue","yellow","green","pink","purple","red","brown")),
                                                                  selectInput("co3GCIM","-log10(P) curve color2 (only for F2):",choices=c("green","black","blue","yellow","red","pink","purple","gray50","brown"))
                                                           )
                                                           
                                                          
                                          ),
                                          conditionalPanel("input.Mose == 'Draw plot'",
                                                           
                                                           column(6,
                                                                  br(),
                                                                  
                                                                  selectInput("plotpopGCIM", "Select population",choices = c("DH BC1 BC2 RIL","F2"))   
                                                           ),                       
                                                           column(6,
                                                                  br(),
                                                                  
                                                                  
                                                                  fileInput("fileplotGCIM", "Input file to draw plot",multiple = TRUE,accept = ".xlsx")
                                                           ),   
                                                           
                                                           column(12,
                                                                  radioButtons("Plotformat2","Plot format",choices=c("*.png","*.tiff","*.jpeg","*.pdf"),inline = "TRUE")
                                                           ),
                                                           column(12,
                                                                  downloadButton("downloadplotGCIM", "Download plot",style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                                           ),
                                                           column(12,
                                                                  plotOutput("plotGCIM")
                                                           )
                                          )
                                        )
                               )
                   )
                 )
               )
      )
    )
  )
  
   shinyApp(ui<-ui,server<-server)
}

