add_stats <- function(dataset,plotx,means,ymin,ymax,predictor,p_colony,output,contrast.matrix,survival,p_queen=NA,lab_title){
  model <- output[["modellist"]][["24"]][["1"]]
  p_interaction <- output[["interaction_problist"]][["24"]][["1"]]
  ####reduce contrast.matrix
  comparisons <- data.frame(name=rownames(contrast.matrix),stringsAsFactors = F) ####get all the comparisons
  comparisons["term1"] <- unlist(lapply(comparisons$name,function(x)unlist(strsplit(x,split=" - "))[1]))
  comparisons["term2"] <- unlist(lapply(comparisons$name,function(x)unlist(strsplit(x,split=" - "))[2]))
  
  comparisons["treatment_index"] <- unlist(lapply(comparisons$term1,function(x)which(unlist(strsplit(x,split="_"))%in%treatments)))
  treatment_index <- unique(comparisons$treatment_index)
  
  comparisons["treatment1"] <- unlist(lapply(comparisons$term1,function(x)unlist(strsplit(x,split="_"))[treatment_index]))
  comparisons["treatment2"] <- unlist(lapply(comparisons$term2,function(x)unlist(strsplit(x,split="_"))[treatment_index]))
  
  comparisons["level1"] <- unlist(lapply(comparisons$term1,function(x)paste(unlist(strsplit(x,split="_"))[(treatment_index+1):length(unlist(strsplit(x,split="_")))],collapse="_")))
  comparisons["level2"] <- unlist(lapply(comparisons$term2,function(x)paste(unlist(strsplit(x,split="_"))[(treatment_index+1):length(unlist(strsplit(x,split="_")))],collapse="_")))
  
  comparisons["to_keep"] <- F
  comparisons[which(comparisons$level1==comparisons$level2),"to_keep"] <- T
  comparisons[which(comparisons$level1=="queen"),"to_keep"] <- F
  
  contrast.matrix <- contrast.matrix[which(comparisons$to_keep),]
  ####if more than 2 levels, then add individual change comparisaon statistics
  if (!is.na(p_interaction)&p_interaction <= 0.05){
    
    post_hoc <- summary(glht(model,contrast.matrix),test=adjusted("BH"))
    print("z value");print(post_hoc$test$tstat);print("Pr>|z|");print(post_hoc$test$pvalues);
    p_values <- as.numeric(post_hoc$test$pvalues);names(p_values) <- names(post_hoc$test$coefficients)
    if (!is.na(p_queen)){
      names(p_queen) <- "queen"
      p_values <- c(p_values,p_queen)
    }
    
    for (idx in 1:length(p_values)){
      p_value <- p_values[idx]
      if (names(p_value)=="queen"){level <- "queen"}else{level <- comparisons[which(comparisons$name==names(p_value)),"level1"]}
      print(p_value)
      if (p_value>0.05){p_cex <- par("cex") *inter_cex;adjust_line <- 0;fonty <- 1}else{p_cex <- par("cex") *max_cex*1.1;adjust_line <- -0.3; fonty <-  2}
      mtext(from_p_to_ptext(p_value),side=3,line=stat_line+adjust_line,at=mean(plotx[which(means$predictor==level)]),xpd=T,cex=p_cex,font=fonty)
    }
  }else if  (!is.na(p_interaction)&p_interaction > 0.05){
    anov <- Anova(model)
    if (!all(!grepl("moved",lab_title))){ ###for distance moved we are interested in the full interaction
      p_value <- p_interaction
    }else{
      p_value <- anov["Pr(>Chisq)"]["time:treatment","Pr(>Chisq)"]
    }
    if (p_value>0.05){p_cex <- par("cex") *inter_cex;adjust_line <- 0;fonty <- 1}else{p_cex <- par("cex") *max_cex*1.1;adjust_line <- -0.3; fonty <-  2}
    mtext(from_p_to_ptext(p_value),side=3,line=stat_line+adjust_line,xpd=T,cex=p_cex,font=fonty)
    
  }
}

apply_alpha <- function(coli,alphi){
  ###Convert to rgb
  coli <- as.vector(col2rgb(coli))
  ###Convert to CMYK
  coli <- rgb2cmyk(R=coli[1],G=coli[2],B=coli[3])
  ###Multiply K by alphi
  coli["K"] <- alphi*coli["K"]
  ###Convert back to rgb
  coli <- cmyk2rgb(C=coli["C"],M=coli["M"],Y=coli["Y"],K=coli["K"])
  ###Convert from RGB to hex
  rcol <- str_pad(as.hexmode(coli["R"]), width=2, side="left", pad="0")
  gcol <- str_pad(as.hexmode(coli["G"]), width=2, side="left", pad="0")
  bcol <- str_pad(as.hexmode(coli["B"]), width=2, side="left", pad="0")
  Colourz <- paste("#",rcol,gcol,bcol,sep="" )
  return(Colourz)
}

clean <- function(){
  rm(list=ls(envir = .GlobalEnv)[!ls(envir = .GlobalEnv)%in%to_keep], envir = .GlobalEnv)
  no_print <- gc(verbose=F)
}

closest_match <- function(x,y){
  return(min(which(abs(x-y)==min(abs(x-y))),na.rm=T))
}

colo <- function(x){
  return(paste("colony",paste(rep(0,3-nchar(x)),collapse=""),x,sep=""))
}

cmyk2rgb <- function(C,M,Y,K){
  R <- round(255 * (1-C) * (1-K))
  G <- round(255 * (1-M) * (1-K))
  B <- round(255 * (1-Y) * (1-K))
  outcol <- c(R,G,B)
  names(outcol) <- c("R","G","B")
  return(outcol)
}

partial_least_square_regression <- function(experiments){
  variables <- c("modularity","density","diameter","degree_mean","clustering","task_assortativity","efficiency","colony_size")
  names(variables) <- c("Modularity","Density","Diameter","Mean degree","Clustering","Task assortativity","Network efficiency","Colony size")
  desired_treatments  <- c("pathogen")
  
  #####Read qPCR data
  infection_data <- NULL
  for (experiment in experiments){
    temp <- read.table(paste(disk_path,experiment,"original_data/qPCR/qPCR_results.txt",sep="/"),header=T,stringsAsFactors=F)
    temp["experiment"] <- experiment
    if (!"age"%in%names(temp)){
      temp["age"] <- NA
    }
    temp <- temp[,order(names(temp))]
    infection_data <- rbind(infection_data,temp)
    rm(list="temp")
  }
  infection_data <- infection_data[c("experiment","colony","treatment","tag","status","measured_load_ng_per_uL","above_detection_threshold")]

  infection_data <-infection_data[which(infection_data$treatment%in%desired_treatments),]
  infection_data <-infection_data[which(infection_data$status=="untreated"),]
  infection_data <- infection_data[which(!grepl("brood",infection_data$tag)),]
  
  names(infection_data)[ names(infection_data)=="above_detection_threshold"]<- "contaminated"
  infection_data$contaminated <- as.numeric(infection_data$contaminated)
  replac_val  <- min(infection_data$measured_load_ng_per_uL[infection_data$measured_load_ng_per_uL!=0],na.rm=T)/sqrt(2)
  
  #####Read network data
  network <- NULL
  for (experiment in experiments){
    temp <- read.table(paste(disk_path,experiment,"processed_data/network_properties/post_treatment/network_properties_observed.txt",sep="/"),header=T,stringsAsFactors = F)
    temp["experiment"] <- experiment
    temp <- temp[,order(names(temp))]
    network <- rbind(network,temp)
  }
  network <- network[which(network$time_hours>=0),]
  network <- network[which(network$treatment%in%desired_treatments),]
  
  prevalence          <- aggregate(contaminated~experiment+colony,FUN=mean,data=infection_data)
  names(prevalence)    <- c("experiment","colony","prevalence")
  intensity           <- aggregate(log10(measured_load_ng_per_uL+replac_val)~experiment+colony,FUN=mean,data=infection_data)
  names(intensity)    <- c("experiment","colony","mean_received_load")
  
  network <- aggregate(cbind(modularity, density,diameter,degree_mean,clustering,task_assortativity,efficiency)~experiment+colony+treatment+colony_size,FUN=mean,data=network)
  
  network <- merge(network,prevalence)
  network <- merge(network,intensity)
  
  network$contaminated            <- as.numeric(network$prevalence)
  network$colony_size             <- as.numeric(network$colony_size )
  network$mean_received_load      <- as.numeric(network$mean_received_load )

    ####is the lmer model acceptable?
    model <- lm(mean_received_load ~ modularity+ density+diameter+degree_mean+clustering+task_assortativity+efficiency+colony_size,data=network)
    print("Properties' VIF:")
    print(vif(model))
    

    for (variable in c("mean_received_load","prevalence")){
      network["variable"] <- network[,variable]
     
      ###determine number of components (number of dimensions with lowest cross validation error) 
      pls.model <- plsr(variable ~ modularity+density+diameter+degree_mean+clustering+task_assortativity+efficiency+colony_size, ncomp = 8, data = network, validation = "CV",scale=T)
      cv = RMSEP(pls.model)
      best.dims = which.min(cv$val[estimate = "adjCV", , ][2:length(cv$val[estimate = "adjCV", , ])]) 
     
      ###refit pls with desired number of components 
      pls.model <- plsr(variable ~ modularity+density+diameter+degree_mean+clustering+task_assortativity+efficiency+colony_size, ncomp = as.numeric(best.dims), data = network,scale=T)
      
      ###extract model coefficients
      coefficients = coef(pls.model)
      
      ###sort coefficients from lowest to highest
      coefficients = sort(coefficients[, 1 , 1])
      
      ###define colours: negative coefficients in black, positive coefficients in white
      colorz <- rep(NA,length(coefficients))
      colorz[coefficients<0] <- "black"
      colorz[coefficients>0] <- "white"
      
      ###Plot coefficients
      par_mar_ori <- par("mar")
      par_mgp_ori <- par("mgp")
      par(mar=par("mar")+c(1.5,0,0,0))
      par(mgp=par("mgp")+c(0.1,0,0))
      par(cex.lab=inter_cex)
      barplot(coefficients,ylab="Partial least square regression coefficients",cex.axis = min_cex,cex.names = min_cex,las=2,names.arg=names(variables[ match(names(coefficients),variables)]),col=colorz)
      abline(h=0)
      if (variable=="mean_received_load"){
        mtext(expression(bold(paste("Mean measured load (ng/", mu, "L)"))),line=1,side=3,cex=par("cex")*max_cex,font=2)
        ylab <- "Mean fitted load"
        xlab <- "Mean measured load"
      }else{
        mtext("Measured prevalence (Prop. untreated workers)",line=1,side=3,cex=par("cex")*max_cex,font=2)
        ylab <- "Fitted prevalence"
        xlab <- "Measured prevalence"
      }
      par(mar=par_mar_ori)
      par(mgp=par_mgp_ori)
      
      ###Plot fit
      fitted <- predict(pls.model, ncomp = best.dims, newdata = network)
      measured <- network$variable
      ymin <- min(floor(fitted*2)/2)
      ymax <- max(ceiling(fitted*2)/2)
      
      xmin <- min(floor(measured*2)/2)
      xmax <- max(ceiling(measured*2)/2)
      
      plot(fitted~measured,bty="l",xlim=c(xmin,xmax),ylim=c(ymin,ymax),ylab=ylab,xlab=xlab,cex.axis=min_cex,cex.lab=inter_cex,yaxt="n",xaxt="n",pch=16)
      where <- axisTicks(c(par("usr")[3],par("usr")[4]),log=F)
      where <- where[which(where==round(where))]
      axis(2,at=where,labels=format(10^(where),scientific=T),cex.lab=inter_cex,cex.axis=min_cex,lwd=0,lwd.ticks=1)
      where <- axisTicks(c(par("usr")[1],par("usr")[2]),log=F)
      where <- where[which(where==round(where))]
      axis(1,at=where,labels=format(10^(where),scientific=T),cex.lab=inter_cex,cex.axis=min_cex,lwd=0,lwd.ticks=1)
      model <- lm(fitted~measured)
      abline(a=coef(model)["(Intercept)"],b=coef(model)["measured"],xpd=F,col="red")
      pval <- Anova(model)["Pr(>F)"]["measured","Pr(>F)"]
      title(from_p_to_ptext(pval))  
    }
}

create_diff <- function(dataset,predictor,type,form_stat=NULL,collective,plot_untransformed=F,diff_type="absolute_difference"){
  if (is.null(dataset)){
    return(NULL)
  }else{
    #################################################################################
    ####plot
    #################################################################################
    befores <- dataset[dataset$time=="Before",];afters <-dataset[dataset$time=="After",]
    # print(plot_untransformed)
    if (!plot_untransformed){
      names(befores)[names(befores)=="variable"] <- "variable_before";names(afters)[names(afters)=="variable"] <- "variable_after"
    }else{
      names(befores)[names(befores)=="untransformed_variable"] <- "variable_before";names(afters)[names(afters)=="untransformed_variable"] <- "variable_after"
    }
    if ((!all(!grepl("time:treatment:predictor",as.character(form_stat))))|(!all(!grepl("time \\* treatment \\* predictor",as.character(form_stat))))){
      befores <- aggregate(na.rm=T,na.action="na.pass",variable_before~colony_size+colony+predictor+treatment+antid+time_of_day+time_of_day_bis,FUN=mean,data=befores) 
      afters <- aggregate(na.rm=T,na.action="na.pass",variable_after~colony_size+colony+predictor+treatment+antid+time_of_day+time_of_day_bis,FUN=mean,data=afters) 
    }else{
      if (!grepl("age",root_path)){
        if (collective|type!="individual"){
          befores <- aggregate(na.rm=T,na.action="na.pass",variable_before~colony_size+colony+predictor+time_of_day+time_of_day_bis,FUN=mean,data=befores) 
          afters <- aggregate(na.rm=T,na.action="na.pass",variable_after~colony_size+colony+predictor+time_of_day+time_of_day_bis,FUN=mean,data=afters) 
        }else{
          befores <- aggregate(na.rm=T,na.action="na.pass",variable_before~colony_size+colony+predictor+antid+time_of_day+time_of_day_bis,FUN=mean,data=befores) 
          afters <- aggregate(na.rm=T,na.action="na.pass",variable_after~colony_size+colony+predictor+antid+time_of_day+time_of_day_bis,FUN=mean,data=afters) 
        }  
      }else{
        if (collective){
          befores <- aggregate(na.rm=T,na.action="na.pass",variable_before~colony_size+colony+treatment+time_of_day+time_of_day_bis,FUN=mean,data=befores) 
          afters <- aggregate(na.rm=T,na.action="na.pass",variable_after~colony_size+colony+treatment+time_of_day+time_of_day_bis,FUN=mean,data=afters) 
        }else{
          if(type!="individual"){
            befores <- aggregate(na.rm=T,na.action="na.pass",variable_before~colony_size+colony+predictor+treatment+time_of_day+time_of_day_bis,FUN=mean,data=befores) 
            afters <- aggregate(na.rm=T,na.action="na.pass",variable_after~colony_size+colony+predictor+treatment+time_of_day+time_of_day_bis,FUN=mean,data=afters) 
          }else{
            befores <- aggregate(na.rm=T,na.action="na.pass",variable_before~colony_size+colony+predictor+treatment+antid+time_of_day+time_of_day_bis,FUN=mean,data=befores) 
            afters <- aggregate(na.rm=T,na.action="na.pass",variable_after~colony_size+colony+predictor+treatment+antid+time_of_day+time_of_day_bis,FUN=mean,data=afters) 
          }
        }  
      }
    }
    befores["average_before"] <- mean(befores$variable_before,na.rm=T)
    diff <- merge(befores,afters,all=T)
    if (diff_type=="absolute_difference"){
      diff["variable"] <- diff$variable_after-diff$variable_before;diff["time"] <- "diff"
    }
    if (diff_type=="normalised_by_average_before"){
      diff["variable"] <- diff$variable_after-(diff$variable_before-diff$average_before);diff["time"] <- "diff"
    }
    if (diff_type=="relative_difference"){
      diff["variable"] <- diff$variable_after-diff$variable_before; diff["time"] <- "diff"
      if ("treatment"%in%names(diff)){
        diff_control <- diff[which(as.character(diff$treatment)=="control"),]
        diff <- diff[which(!as.character(diff$treatment)=="control"),]
        diff_control <- aggregate(na.rm=T,na.action="na.pass",variable~predictor,FUN=get(relative_function),data=diff_control)
        names(diff_control)[names(diff_control)=="variable"] <- "mean_diff_control"
        diff <- merge(diff,diff_control,all.x=T)
        diff$variable <- diff$variable-diff$mean_diff_control
        diff <- diff[,which(names(diff)!="mean_diff_control")]
      }else{
        diff_control <- diff[which(as.character(diff$predictor)=="control"),]
        diff <- diff[which(!as.character(diff$predictor)=="control"),]
        
        diff_control <- aggregate(na.rm=T,na.action="na.pass",variable~1,FUN=get(relative_function),data=diff_control)
        names(diff_control)[names(diff_control)=="variable"] <- "mean_diff_control"
        
        diff <- merge(diff,diff_control,all.x=T)
        diff$variable <- diff$variable-diff$mean_diff_control
        diff <- diff[,which(names(diff)!="mean_diff_control")]
        
      }
      
      
    }
    if (!grepl("age",root_path)){diff$predictor <- factor(diff$predictor)}else{diff$treatment <- factor(diff$treatment)}
    
    return(diff)
  }
}

from_p_to_ptext <- function(pvalue){
  if (is.na(pvalue)){
    pvaluetext <- "p = NA"
  }else{
    if 
    # (pvalue< 0.00001){
    #   pvaluetext <- "*****"
    # }else if (pvalue< 0.0001){
    #   pvaluetext <- "****"
    # }else if 
    (pvalue< 0.0001){
      pvaluetext <- "***"
    }else if (pvalue< 0.005){
      pvaluetext <- "**"
    }else if (pvalue< 0.05){
      pvaluetext <- "*"
    }else if (pvalue< 0.1){        
      pvaluetext <- paste("p=",sprintf("%.3f", pvalue),sep="")
    }else if (pvalue< 1){ 
      pvaluetext <- paste("p=",sprintf("%.2f", pvalue),sep="")
    }else{
      pvaluetext <- "p = 1"
    }#if (pvalue< 10^-6)
  }
  return(pvaluetext)
}

GetColorHex <- function(color){
  clr <- col2rgb(color)
  hex_and_col <- sprintf("#%02X%02X%02X %3d %3d %3d", clr[1],clr[2],clr[3], clr[1], clr[2], clr[3])
  hex <- unlist(strsplit(hex_and_col,split=" "))[1]
  return(hex)
}

get_posthoc_groups <- function(model,matrix,levels,randy,dataset){
  post_hoc <- summary(glht(model,matrix),test=adjusted("BH"))
  print("z value");print(post_hoc$test$tstat);print("Pr>|z|");print(post_hoc$test$pvalues);
  p_values <- as.numeric(post_hoc$test$pvalues);names(p_values) <- names(post_hoc$test$coefficients)
  
  post_hoc_levels <- names(levels)
  post_hoc_mat <- matrix(NA,nrow=length(post_hoc_levels)-1,ncol=length(post_hoc_levels)-1)
  rownames(post_hoc_mat) <- post_hoc_levels[2:length(post_hoc_levels)]
  colnames(post_hoc_mat) <- post_hoc_levels[1:(length(post_hoc_levels)-1)]
  for (i in 1:nrow(post_hoc_mat)){
    for (j in 1:i){
      if (!is.null(randy)){
        post_hoc_mat[i,j] <- as.logical(as.numeric(p_values[paste(randy," ",colnames(post_hoc_mat)[j]," minus ",randy," ",rownames(post_hoc_mat)[i],sep="")])>0.05)
      }else{
        post_hoc_mat[i,j] <- as.logical(as.numeric(p_values[paste(colnames(post_hoc_mat)[j]," - ",rownames(post_hoc_mat)[i],sep="")])>0.05)
      }
    }
  }
  g <- post_hoc_mat
  g <- cbind(rbind(NA, g), NA)
  g <- replace(g, is.na(g), FALSE)
  g <- g + t(g)
  diag(g) <- 1
  n <- length(post_hoc_levels)
  rownames(g) <- 1:n
  colnames(g) <- 1:n
  #g
  same <- which(g==1)
  topology <- data.frame(N1=((same-1) %% n) + 1, N2=((same-1) %/% n) + 1)
  topology <- topology[order(topology[[1]]),] # Get rid of loops and ensure right naming of vertices
  g3 <- simplify(graph.data.frame(topology,directed = FALSE))
  #get.data.frame(g3)
  #plot(g3)
  res <- maximal.cliques(g3)
  
  # Reorder given the smallest level
  clique_value <- NULL
  means <- aggregate(variable~treatment+predictor,FUN=mean,data=dataset)
  means$predictor <- names(levels)
  for (i in 1:length(res)){
    clique_value <- c(clique_value,mean(means[as.numeric(unlist(res[[i]])),"variable"]))
  }
  res <- res[order(clique_value)]
  
  # Get group letters
  lab.txt <- vector(mode="list", n)
  lab <- letters[seq(res)]
  for(i in seq(res)){
    for(j in res[[i]]){
      lab.txt[[j]] <- paste0(lab.txt[[j]], lab[i])
    }
  }
  post_hoc_groups <- unlist(lab.txt); names(post_hoc_groups) <- levels[post_hoc_levels]
  return(post_hoc_groups)
}

is.even <- function(x) {
  return(x %% 2 == 0)
}

log_transf <- function(x){
  if (all(x>0)){
    replac_val <- 0
  }else if (all(x>=0)){
    replac_val <- (min(x[x!=0],na.rm=T))/sqrt(2)
  }else{
    replac_val_1 <- -min(x,na.rm=T)
    y <- x+replac_val_1
    replac_val <- replac_val_1 + (min(y[y!=0],na.rm=T))/sqrt(2)
  }
  return(log10(x+replac_val))
}

meta_analysis <- function(p_values,effects,std.errors){
  ####Get effect signs
  effect_signs <- sign(effects)
  
  ########## Create a vector of one-sided p-values for each possible effect direction each side
  p_values_side1 <- p_values ###one-sided p testing whether observed is lower than random 
  p_values_side2 <- 1-p_values ###one-sided p testing whether observed is higher than random 
  
  ############# Meanp method
  ####Check which of the two effect direction is the one to keep
  best_idx <- which(c(meanp(p_values_side1)$p,meanp(p_values_side2)$p)== min(meanp(p_values_side1)$p,meanp(p_values_side2)$p))
  ####Get one-sided combined p-value and test statistics
  p_value <-c(meanp(p_values_side1)$p,meanp(p_values_side2)$p)[best_idx]
  z <- c(meanp(p_values_side1)$z,meanp(p_values_side2)$z)[best_idx]
  ####convert back to 2-sided pvalue
  two_sided_p_value <-  p_value*2
  
  p_values_meta_analysis <- data.frame(meta_statistic = z,one_sided_p=p_value,two_sided_p=two_sided_p_value)
  return(p_values_meta_analysis)
}

mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}

mysquare <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   squares=2*size, add=TRUE, inches=FALSE)
         })
}

normalize_to_range <- function(Input, Min, Max){
  Range  <- max(Input) - min(Input)
  Input  <- (Input - min(Input)) / Range
  Range2 <- Max - Min
  Input  <- (Input * Range2) + -1
  return(Input)
}

perform_barplot_analysis <- function(root_path,collective=F,dataset,lab_title,type="individual",excluded=NULL,survival=F,pool=F,violin_params=NULL,pool_plot=F,adjust_title_line=0){
  par_mar_ori <- par()$mar
  for (tab in c("dataset","excluded")){
    table <- get(tab)
    if (!is.null(table)){
      if (!grepl("age",root_path)){predictor <- "treatment";table["predictor"]<-table$treatment;inter <- "time*predictor";col_vector <- treatment_colours}else{predictor <- "status";table["predictor"]<-table$status;inter <- "time";col_vector <- statuses_colours}
      assign(tab,table)
    }
  }
  original_dataset <- dataset
  
  #######simple barplot, comparison with full last 24h period, meaning there are some times that are compared with minus 48h
  dataset <- original_dataset
  if(!exists("prepare_stats")){
    prepare_stats <- prepare_stats_1(collective,dataset,type,predictor,inter,survival); form_stat <- prepare_stats[["form_stat"]]; p_colony <- prepare_stats[["p_colony"]];contrast.matrix <- prepare_stats[["contrast.matrix"]]
  }
  temp <- transform_dataset(dataset,cut=F,predictor=predictor,form_stat,excluded); dataset <- temp[["dataset"]]; form_stat <- temp[["form_stat"]];excluded <- temp[["excluded"]]; rm(list=c("temp"))
  output <- prepare_stats_2(dataset,form_stat,survival=survival)
  diff <- create_diff(dataset,predictor,type,collective=collective,diff_type="relative_difference")
  
  formula_plot <- update(form_stat,.~.-(1 | colony)-(1 | time_of_day_bis)-(1|time_of_day)-(1|antid)-colony_size)
  diff["predictor_plot"] <- as.character(diff$predictor)
  diff$predictor <- as.character(diff$predictor)
  
  par(mar=par_mar_ori+c(-0.8,-2+0.1+0.5,0,0.5))
  plot_diff(diff,lab_title,col_vector,predictor,form_stat,p_colony,contrast.matrix,output,dataset,collective=collective,violin_params = violin_params,adjust_title_line=adjust_title_line)
  
  
  
  
  
  par(mar=par_mar_ori)
}

perform_barplot_analysis_refined <- function(root_path,collective=NULL,dataset,lab_title,type=NULL,survival=F,excluded=NULL,pool=F,violin_params=NULL,pool_plot=F,plot_untransformed=F,aligned=F,adjust_title_line=0){
  par_mar_ori <- par()$mar
  col_vector <- statuses_colours
  ###sort predictor levels
  levels <- sort(unique(dataset$predictor))
  levels <- levels[order(match(levels,status_order))]
  dataset$predictor <- factor(dataset$predictor,levels = levels)
  levels_treatment <- unique(as.character(dataset$treatment))
  levels_treatment <- levels_treatment[order(match(levels_treatment,status_order))]
  dataset$treatment <- factor(dataset$treatment,levels = levels_treatment)
  
  
  queens <- dataset[which(dataset$predictor=="queen"),]
  if(nrow(queens)>0){
    print("Queens")
    if (!grepl("age",root_path)){predictory <- "treatment";queens["predictor"]<-queens$treatment;inter <- "time*predictor";col_vector <- statuses_colours}else{predictory <- "status";queens["predictor"]<-queens$status;inter <- "time";col_vector <- statuses_colours}
    
    prepare_stats <- prepare_stats_1(collective,queens,type="individual",predictor=predictory,inter="time*predictor",survival); form_stat <- prepare_stats[["form_stat"]]; p_colony <- prepare_stats[["p_colony"]];contrast.matrix <- prepare_stats[["contrast.matrix"]]
    temp <- transform_dataset(queens,cut=F,predictor=predictory,form_stat,excluded); queens <- temp[["dataset"]]; form_stat <- temp[["form_stat"]];excluded <- temp[["excluded"]]; rm(list=c("temp"))
    reduced_contrast_matrix <- rbind(contrast.matrix[1,])
    row.names(reduced_contrast_matrix) <- row.names(contrast.matrix)[1]
    output <- prepare_stats_2(queens,form_stat,survival,T,reduced_contrast_matrix)
    p_queen <- as.numeric(output$interaction_problist[["24"]])
  }else{
    p_queen <- NA
  }
  
  prepare_stats <- prepare_stats_3(dataset,predictor,survival); form_stat <- prepare_stats[["form_stat"]]; p_colony <- prepare_stats[["p_colony"]];contrast.matrix <- prepare_stats[["contrast.matrix"]]
  output <- prepare_stats_4(dataset,form_stat,survival)
  dataset["time_of_day_bis"] <- dataset$time_of_day
  diff <- create_diff(dataset,predictor,type="individual",form_stat=form_stat,collective=collective,plot_untransformed=plot_untransformed,diff_type="relative_difference")
  diff["time"] <- "Delta"
  ylab <- lab_title
  par(mar=par_mar_ori+c(-0.3,-2+0.1+0.5,0,0.5))
  if (pool_plot){
    diff <- aggregate(na.rm=T,na.action="na.pass",variable~predictor+colony_size+colony+treatment+antid+time,FUN=mean,data=diff)
  }
  plot_refined(diff,lab_title,col_vector,predictor,output,contrast.matrix,survival,dataset,p_queen,violin_params=violin_params,aligned=aligned,adjust_title_line=adjust_title_line)
  
  par(mar=par_mar_ori)
}

perform_barplot_analysis_simple <- function(root_path,collective=NULL,dataset,lab_title,type=NULL,survival=F,excluded=NULL,pool=F,violin_params=NULL,pool_plot=F,plot_untransformed=F,aligned=F,adjust_title_line=0){
  par_mar_ori <- par()$mar
  col_vector <- statuses_colours
  ###sort predictor levels
  levels <- sort(unique(dataset$predictor))
  levels <- levels[order(match(levels,status_order))]
  dataset$predictor <- factor(dataset$predictor,levels = levels)
  levels_treatment <- unique(as.character(dataset$treatment))
  levels_treatment <- levels_treatment[order(match(levels_treatment,status_order))]
  dataset$treatment <- factor(dataset$treatment,levels = levels_treatment)
  print(par_mar_ori)
  par(mar=par_mar_ori+c(-0.3,-0.1,0,0))
  if (!is.null(violin_params)){
    violin_params <- as.numeric(unlist(violin_params))
    ##read violin param
    range <- violin_params[1]
    ylim_fac1 <- violin_params[2]
    ylim_fac2 <- violin_params[3]
    wex <- violin_params[4]
    h <- violin_params[5]
    
  }
  
  level_names <- expand.grid(levels(dataset$treatment),levels(dataset$predictor))
  level_names <- within(level_names,Var3<-paste(Var1,Var2,sep="."))$Var3
  names(level_names) <- level_names
  
  if (length(level_names)==6){
    contrast.matrix <- rbind(
      "1 - 2"=c(0,-1,0,0,0,0,0),
      "1 - 3"=c(0,0,-1,0,0,0,0),
      "1 - 4"=c(0,-1,-1,0,0,-1,0),
      "1 - 5"=c(0,0,0,-1,0,0,0),
      "1 - 6"=c(0,-1,0,-1,0,0,-1),
      "2 - 3"=c(0,1,-1,0,0,0,0),
      "2 - 4"=c(0,0,-1,0,0,-1,0),
      "2 - 5"=c(0,1,0,-1,0,0,0),
      "2 - 6"=c(0,0,0,-1,0,0,-1),
      "3 - 4"=c(0,-1,0,0,0,-1,0),
      "3 - 5"=c(0,0,1,-1,0,0,0),
      "3 - 6"=c(0,-1,1,-1,0,0,-1),
      "4 - 5"=c(0,1,1,-1,0,1,0),
      "4 - 6"=c(0,0,1,-1,0,1,-1),
      "5 - 6"=c(0,-1,0,0,0,0,-1)
    )
  }else if (length(level_names)==4){
    contrast.matrix <- rbind(
      "1 - 2" = c(0,-1,0,0,0),
      "1 - 3" = c(0,0,-1,0,0),
      "1 - 4" = c(0,-1,-1,0,-1),
      "2 - 3" = c(0,1,-1,0,0),
      "2 - 4" = c(0,0,-1,0,-1),
      "3 - 4" = c(0,-1,0,0,-1)
    )
  }
  for (i in 1:length(level_names)){
    rownames(contrast.matrix) <- gsub(i,level_names[i],rownames(contrast.matrix))
  }
  if ("antid"%in%names(dataset)){
    form_stat <- as.formula(paste("variable~", paste(c("treatment*predictor","colony_size","(1|colony)","(1|antid)"), collapse= "+")))
  }else{
    form_stat <- as.formula(paste("variable~", paste(c("treatment*predictor","colony_size","(1|colony)"), collapse= "+")))
  }
  model <- do.call(lmer, list(formula=form_stat, data=dataset))
  test_norm(residuals(model))
  anov <- Anova(model)
  print(anov)
  posthoc_groups <- get_posthoc_groups(model=model,matrix=contrast.matrix,levels=level_names,randy=NULL,dataset=dataset)
  
  for (randy in  levels_treatment){
    temp <- posthoc_groups[grepl(randy,names(posthoc_groups))]; names(temp) <- gsub(paste(randy,"\\.",sep=""),"",names(temp))
    assign(paste("post_hoc_",randy,sep=""),temp)
  }
  
  means <- aggregate(na.rm=T,na.action="na.pass",variable~treatment+predictor,FUN="mean",data=dataset);ses <- aggregate(na.rm=T,na.action="na.pass",variable~treatment+predictor,FUN="std.error",data=dataset);
  names(means)[names(means)=="variable"] <- "mean";names(ses)[names(ses)=="variable"] <- "se";means <- merge(means,ses)
  means <- means[order(match(means$treatment,levels(dataset$treatment)),match(means$predictor,levels(dataset$predictor))),]
  to_plot <- unique(means[c("treatment","predictor")])
  
  
  means[is.na(means$se),"se"] <- 0
  ymin <- min(c(means$mean-means$se),na.rm=T);ymax<- max(c(means$mean+means$se),na.rm=T)
  if (ymin>0){ymin <- 0};if (ymax<0){ymax <- 0}
  if (grepl("point",plot_type)|grepl("boxplot",plot_type)|grepl("violinplot",plot_type)){
    rangy <- max(dataset$variable,na.rm=T)-min(dataset$variable,na.rm=T)
    ymin <- min(ymin, min(dataset$variable,na.rm=T)-0.1*rangy)
    ymax <- max(ymax, max(dataset$variable,na.rm=T)+0.1*rangy)
    
  }else{
    rangy <- ymax-ymin
    ymax <- ymax+0.25*rangy
  }
  
  
  barwidth <- 0.5; barwidth_fac_within <- 0.5; barwidth_fac_between <- 2
  barspace <- rep(c(barwidth_fac_between,rep(barwidth_fac_within,(length(unique(means$predictor))-1))),length(unique(means$treatment)))
  col_vector <- statuses_colours
  
  means["alpha"] <- as.numeric(alphas[as.character(means$treatment)])
  means["full_col"] <- col_vector[as.character(means$predictor)]
  means[which(as.character(means$predictor)=="outdoor_ant"),"alpha"] <- means[which(as.character(means$predictor)=="outdoor_ant"),"alpha"]*10
  # for (colidx in 1:nrow(means)){
  #   means[colidx,"final_col"] <- apply_alpha(means[colidx,"full_col"],means[colidx,"alpha"])
  # }
  means["final_col"] <- means$full_col
  plotx <- barplot(means$mean,plot=F,width=barwidth,space=barspace,lwd=line_max)
  
  if (grepl("bars",plot_type)){
    plotx <- barplot(means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)- 0.75*(barwidth_fac_between*barwidth/2+barwidth/2),max(plotx)+ 0.75*(barwidth_fac_between*barwidth/2+barwidth/2)),xlab="",ylab=lab_title,bty="n",col=means$final_col,xaxs="i",yaxs="i",cex=0.5,cex.axis=min_cex,cex.lab=inter_cex,width=barwidth,space=barspace,xaxt="n",lwd=line_max)
    plot_arrows(means=means,plotx=plotx,plot_type=plot_type,LWD=line_max,LENGTH=0.025,colz=means$final_col)
  }else if (grepl("boxplot",plot_type)){
    ####empty plot
    plot(plotx,means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)- (barwidth_fac_between*barwidth/2+barwidth/2),max(plotx)+ (barwidth_fac_between*barwidth/2+barwidth/2)),xlab="",ylab=lab_title,bty="n",xaxs="i",yaxs="i",type="n",lwd=line_min,cex.axis=min_cex,cex.lab=inter_cex,xaxt="n")
    par(bty="n")
    for (lvly in 1:nrow(means)){
      boxplot(dataset[which(dataset$treatment==means[lvly,"treatment"]&dataset$predictor==means[lvly,"predictor"]),"variable"],at=plotx[lvly],add=T,range=1.5,notch=T,names=F,col=means[lvly,"final_col"],xlab="",ylab="",xaxt="n",xaxs="i",yaxs="i",cex.axis=min_cex,cex.lab=inter_cex,xaxt="n",yaxt="n",medlwd=line_max,boxlwd=line_min+0.5*(line_max-line_min),whisklwd=line_min,whisklty=1,staplelwd=line_min,boxwex=0.7,pch=16)
    }
    
  }else if (grepl("violinplot",plot_type)){
    ####empty plot
    plot(plotx,means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)- (barwidth_fac_between*barwidth/2+barwidth/2),max(plotx)+ (barwidth_fac_between*barwidth/2+barwidth/2)),xlab="",ylab=lab_title,bty="n",xaxs="i",yaxs="i",type="n",lwd=line_min,cex.axis=min_cex,cex.lab=inter_cex,xaxt="n")
    par(bty="n")
    for (lvly in 1:nrow(means)){
      
      if (is.na(range)){
        VioPlot(na.omit(dataset[which(dataset$treatment==means[lvly,"treatment"]&dataset$predictor==means[lvly,"predictor"]),"variable"]),col=alpha(means[lvly,"final_col"],0.7), horizontal=F, at=plotx[lvly], add=TRUE,lty=1, rectCol="black",wex=wex,border=NA,lwd=line_min,mode="Mean")
      }else{
        VioPlot(na.omit(dataset[which(dataset$treatment==means[lvly,"treatment"]&dataset$predictor==means[lvly,"predictor"]),"variable"]),range=range, h=h,col=alpha(means[lvly,"final_col"],0.7), horizontal=F, at=plotx[lvly], add=TRUE,lty=1, rectCol="black",wex=wex,border=NA,lwd=line_min,mode="Mean")
      }
    }
  }else{
    ####empty plot
    plot(plotx,means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)- (barwidth_fac_between*barwidth/2+barwidth/2),max(plotx)+ (barwidth_fac_between*barwidth/2+barwidth/2)),xlab="",ylab=lab_title,bty="n",xaxs="i",yaxs="i",type="n",lwd=line_min,cex.axis=min_cex,cex.lab=inter_cex,xaxt="n")
    ####arrows
    plot_arrows(means=means,plotx=plotx,plot_type=plot_type,LWD=line_max,LENGTH=0.025,colz=means$final_col)
    ####points
    points(plotx,means$mean,col=means$final_col,cex=1.5*max_cex,pch=16,lwd=line_min)
  }
  
  par(xpd=T)
  if (all(nchar(full_statuses_names[as.character(means$predictor)])==1)){
    lassy <- 1
  }else{
    lassy <- 2
  }
  labels_even <- full_statuses_names[as.character(means$predictor)][2*(1:floor(length(full_statuses_names[as.character(means$predictor)])/2))]
  labels_odd  <- full_statuses_names[as.character(means$predictor)][(2*(1:ceiling(length(full_statuses_names[as.character(means$predictor)])/2)))-1]
  axis(1,at=plotx[2*(1:floor(length(plotx)/2))],labels=labels_even,tick=F,cex.axis=min_cex,las=lassy)
  axis(1,at=plotx[(2*(1:ceiling(length(plotx)/2)))-1],labels=labels_odd,tick=F,cex.axis=min_cex,las=lassy)
  par(xpd=F)
  par(xpd=F)
  abline(h=0)
  pvalue <- anov["treatment:predictor","Pr(>Chisq)"]
  for (i in 1:(-1+length(unique(means$treatment)))){
    abline(v= mean(plotx[c(max(which(means$treatment==unique(means$treatment)[i])), min(which(means$treatment==unique(means$treatment)[i+1])))]),lwd=line_max)
  }
  par(xpd=T)
  for (predy in unique(means$treatment)){
    textx <- mean(plotx[which(means$treatment==predy)])
    mtext(full_statuses_names[predy],side=3,line=stat_line,adj=0.5,at=textx,cex=par("cex") *inter_cex,col="black",font=1)
    
    for (idx in 1:length(get(paste("post_hoc_",predy,sep="")))){
      group <- get(paste("post_hoc_",predy,sep=""))[idx]
      mtext(group,side=3,line=stat_line-0.7,at=plotx[which(means$treatment==predy&means$predictor==names(get(paste("post_hoc_",predy,sep="")))[idx])],xpd=T,cex=par("cex") *inter_cex,font=2)
    }
    
  }
  
  
  par(mar=par_mar_ori)
}

plot_age_dol <- function(experiments){
  for (experiment in experiments){
    ####first plot time outside = f(age) ######
      if (grepl("age",experiment)){
        ######read time investment data
        data <- read.table(paste(disk_path,"/",experiment,"/processed_data/time_investment/time_investment.txt",sep=""),header=T,stringsAsFactors = F)
        ######modify time investment data to summarise behaviour over the entire 24-hour period
        data                   <- aggregate(na.rm=T,na.action="na.pass",cbind(outside,detected)~colony+colony_size+treatment+tag+age+status+period,FUN=sum,data=data)
        data$prop_time_outside <- data$outside/data$detected
        data$antid             <- as.character(interaction(data$colony,data$tag))
        ###remove queen
        data                   <- data[which(data$tag!=queenid),]
        ####list desired variables and transformations
        variable_list <- c("prop_time_outside")
        names(variable_list) <- c("Prop. of time outside")
        predictor_list <- c("age")
        names(predictor_list) <- c("Worker age (weeks)")
        transf_variable_list <- c("sqrt")
        transf_predictor_list <- c("none")
        
        analysis <- list(variable_list=variable_list,
                         predictor_list=predictor_list,
                         transf_variable_list=transf_variable_list,
                         transf_predictor_list=transf_predictor_list,              
                         violin_plot_param = list(c(1,-0.02,-0.02,3,0.1)))
        
        ####plot
        plot_regression(data=data,time_point="before",analysis=analysis,n_cat_horiz=15,n_cat_vertic=30,pool=c(F,F),collective=T,input_color=colour_palette_age)
      }
    #####second plot interaction frequencies, observed vs. random  ##############
    data <- read.table(paste(disk_path,"/",experiment,"/processed_data/collective_behaviour/random_vs_observed/interactions.dat",sep=""),header=T,stringsAsFactors = F)
    if (grepl("age",experiment)){
      variable_list <- c("slope_WW_contact_duration_f_age_diff","intra_caste_over_inter_caste_WW_contact_duration","slope_QW_contact_duration_f_W_age","QNurse_over_QForager_contact_duration")
      names(variable_list) <- c("Slope [Log(contact) = f(delta age)]","Within/Between-task contact","Slope [Log(contact) = f(W age)]","Q-N/Q-F contacts")
    }else{
      variable_list <- c("intra_caste_over_inter_caste_WW_contact_duration","QNurse_over_QForager_contact_duration")
      names(variable_list) <- c("Within/Between-task contact","Q-N/Q-F contacts")
    }
    plot_observed_vs_random(experiments=experiment,variables=variable_list,data_input=data)
    
    #####third plot assortativity ############
    if (grepl("age",experiment)){
      variable_list <- c("age_assortativity","task_assortativity")
      names(variable_list) <- c("Age assortativity","Task assortativity")
    }else{
      variable_list <- c("task_assortativity")
      names(variable_list) <- c("Task assortativity")
    }
    plot_observed_vs_random(experiments=experiment,variables=variable_list,pattern="network_properties",data_path="/processed_data/network_properties/random_vs_observed")
    
    #####fourth plot the community composition  ##############
    par_mar_ori <- par()$mar
    par(mar=par_mar_ori+c(0,0,0,0.5))
    
    
    overall_results <- read.table(paste(disk_path,"/",experiment,"/processed_data/network_properties/random_vs_observed/queen_community.dat",sep=""),header=T,stringsAsFactors = F)
    overall_results$treatment <- overall_results$randy; overall_results$time_hours <- -30; overall_results$time_of_day <- -6; 
    
    if (grepl("age",experiment)){
      variable_list <- c("age","proportion_of_foragers")
      names(variable_list) <- c("Mean W age (weeks)","Prop. foragers")
    }else{
      variable_list <- c("proportion_of_foragers")
      names(variable_list) <- c("Prop. foragers")
    }
    transf_variable_list <- c("none","none")
    predictor_list <- c("in_queen_comm","in_queen_comm")
    names(predictor_list) <- c("Community","Community")
    analysis <- list(variable_list=variable_list,
                     transf_variable_list=transf_variable_list,
                     predictor_list=predictor_list,
                     violin_plot_param = list(c(1,0,-0.02,0.08,0.008),c(1,-0.02,-0.02,0.08,0.008)))
    
    plot_regression(data=overall_results,time_point="comparison",analysis,n_cat_horiz=15,n_cat_vertic=30,pool=c(F,F),prepare=T,status="all",collective=F,pool_plot=F)
    
    
    par(mar=par_mar_ori)
  }
  ####plot dividing lines and titles#######
  fac <- sum(heits)
  par(xpd=NA)
  x_line <- grconvertX(2/3, from='ndc')
  y_line1 <- grconvertY(1, from='ndc')####lower left = 0-0;top_right=1-1
  y_line2 <- grconvertY(0, from='ndc')
  segments(x_line,y_line1,x_line,y_line2)
  
  ####write titles #######
  ##########MAIN TITLES #########
  x_text1 <- grconvertX(2/6, from='ndc')
  x_text2 <- grconvertX(2/3+1/6, from='ndc')
  y_text <- grconvertY(1- (til/fac)/2, from='ndc')
  
  text(x_text1,y_text,labels="Age experiment (n=11)",font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)
  text(x_text2,y_text,labels="Main experiment (n=22)",font=2, cex=max_cex,adj=c(0.5,0.5),xpd=NA)
  
  ##########LETTERS ###############
  x_text1 <- grconvertX(1/80, from='ndc')
  x_text2 <- grconvertX(1/3+1/80, from='ndc')
  x_text3 <- grconvertX(1/3+1/6+1/80, from='ndc')
  x_text4 <- grconvertX(2/3+1/80, from='ndc')
  
  X_text1 <- grconvertX(1/80, from='ndc')
  X_text2 <- grconvertX(1/3 +1/80, from='ndc')
  X_text3 <- grconvertX(2/3+1/80, from='ndc')
  X_text4 <- grconvertX(2/3+1/6+1/80, from='ndc')
  
  y_text1 <- grconvertY(1- (heits[1])/fac- 4*(heits[2]/fac)/5+1/80, from='ndc')
  y_text2 <- grconvertY(1- (sum(heits[c(1:3)]))/fac - 4*(heits[4]/fac)/5+1/80, from='ndc')
  y_text3 <- grconvertY(1- (sum(heits[c(1:5)]))/fac - 4*(heits[6]/fac)/5+1/80, from='ndc')
  
  text(x_text1,y_text1,labels=panel_casse("a"),font=2, cex=max_cex,adj=c(0.5,0),xpd=NA)
  text(x_text2,y_text1,labels=panel_casse("b"),font=2, cex=max_cex,adj=c(0.5,0),xpd=NA)
  text(x_text3,y_text1,labels=panel_casse("c"),font=2, cex=max_cex,adj=c(0.5,0),xpd=NA)
  text(x_text4,y_text1,labels=panel_casse("d"),font=2, cex=max_cex,adj=c(0.5,0),xpd=NA)
  text(X_text1,y_text2,labels=panel_casse("e"),font=2, cex=max_cex,adj=c(0.5,0),xpd=NA)
  text(X_text2,y_text2,labels=panel_casse("f"),font=2, cex=max_cex,adj=c(0.5,0),xpd=NA)
  text(X_text3,y_text2,labels=panel_casse("g"),font=2, cex=max_cex,adj=c(0.5,0),xpd=NA)
  text(X_text4,y_text2,labels=panel_casse("h"),font=2, cex=max_cex,adj=c(0.5,0),xpd=NA)
  text(X_text1,y_text3,labels=panel_casse("i"),font=2, cex=max_cex,adj=c(0.5,0),xpd=NA)
  text(X_text3,y_text3,labels=panel_casse("j"),font=2, cex=max_cex,adj=c(0.5,0),xpd=NA)
  
  ########INDIVIDUAL TITLES #########
  x_text1 <- grconvertX(1/6, from='ndc')
  x_text2 <- grconvertX(1/3+1/12, from='ndc')
  x_text3 <- grconvertX(1/3+1/6+1/12, from='ndc')
  x_text4 <- grconvertX(2/3+1/6, from='ndc')
  
  y_text1 <- grconvertY(1- (heits[1])/fac- 4*(heits[2]/fac)/5+1/160, from='ndc')
  y_text2 <- grconvertY(1- (sum(heits[c(1:3)]))/fac - 4*(heits[4]/fac)/5+1/160, from='ndc')
  y_text3 <- grconvertY(1- (sum(heits[c(1:5)]))/fac - 4*(heits[6]/fac)/5+1/160, from='ndc')
  
  text(x_text1,y_text1,labels="Age-related DoL",font=4, cex=inter_cex,adj=c(0.5,1),xpd=NA)
  text(x_text2,y_text1,labels="Community age",font=4, cex=inter_cex,adj=c(0.5,1),xpd=NA)
  text(x_text3,y_text1,labels="Community task",font=4, cex=inter_cex,adj=c(0.5,1),xpd=NA)
  text(x_text4,y_text1,labels="Community task",font=4, cex=inter_cex,adj=c(0.5,1),xpd=NA)
  
  x_text1 <- grconvertX(1/6, from='ndc')
  x_text2 <- grconvertX(1/3+1/6, from='ndc')
  x_text3 <- grconvertX(2/3+1/12, from='ndc')
  x_text4 <- grconvertX(2/3+1/6+1/12, from='ndc')
  text(x_text1,y_text2,labels="W-W contacts",font=4, cex=inter_cex,adj=c(0.5,1),xpd=NA)
  text(x_text2,y_text2,labels="Q-W contacts",font=4, cex=inter_cex,adj=c(0.5,1),xpd=NA)
  text(x_text3,y_text2,labels="W-W contacts",font=4, cex=inter_cex,adj=c(0.5,1),xpd=NA)
  text(x_text4,y_text2,labels="Q-W contacts",font=4, cex=inter_cex,adj=c(0.5,1),xpd=NA)
  
  x_text1 <- grconvertX(1/6, from='ndc')
  x_text2 <- grconvertX(1/3+1/6, from='ndc')
  x_text3 <- grconvertX(2/3+1/6, from='ndc')
  
  text(x_text1,y_text3,labels="Age assortativity",font=4, cex=inter_cex,adj=c(0.5,1),xpd=NA)
  text(x_text2,y_text3,labels="Task assortativity",font=4, cex=inter_cex,adj=c(0.5,1),xpd=NA)
  text(x_text3,y_text3,labels="Task assortativity",font=4, cex=inter_cex,adj=c(0.5,1),xpd=NA)
  par(xpd=F)
  
}

plot_arrows <- function(means,plotx,plot_type,LWD,LENGTH,colz=NULL,direction="normal"){
  options(warn=-1)
  if (grepl("points",plot_type)){
    LENGTH <- LENGTH*2
  }
  if (grepl("bars",plot_type)){
    if (direction=="normal"){
      arrows_low <- means$mean-1*as.numeric(sign(means$mean)<=0)*means$se
      arrows_high <- means$mean+1*as.numeric(sign(means$mean)>=0)*means$se
    }else{
      arrows_low <- means$mean-1*as.numeric(sign(means$mean)>=0)*means$se
      arrows_high <- means$mean+1*as.numeric(sign(means$mean)<=0)*means$se
    }
    code1 <- which(arrows_high==means$mean&arrows_low<means$mean)
    code2 <- which(arrows_low==means$mean&arrows_high>means$mean)
    code3 <- which(arrows_low==means$mean&arrows_high==means$mean)
    
    arrows (plotx[code1],arrows_low[code1],plotx[code1],arrows_high[code1],code=1,angle=90,col="black",lwd=LWD,length=LENGTH)
    arrows (plotx[code2],arrows_low[code2],plotx[code2],arrows_high[code2],code=2,angle=90,col="black",lwd=LWD,length=LENGTH)
    arrows (plotx[code3],arrows_low[code3],plotx[code3],arrows_high[code3],code=3,angle=90,col="black",lwd=LWD,length=LENGTH)
  }else{
    arrows_low <- means$mean-means$se
    arrows_high <- means$mean+means$se
    arrows (plotx,arrows_low,plotx,arrows_high,code=3,angle=90,col=colz,lwd=1.5*LWD,length=1.5*LENGTH)
  }
  options(warn=0)
}

plot_before_vs_after <- function(root_path,data_path,analysis,status,collective,pool_plot=F,pattern="",plot_untransformed=F,boldy=F,aligned=F){
  setwd(data_path)
  file_list <- list.files(pattern=pattern)
  data <- NULL
  for (file in file_list){
    data <- rbind(data,read.table(file,header=T,stringsAsFactors = F))  
  }
  if (status=="untreated"&"status"%in%names(data)){
    data <- data[which(data$status=="untreated"),]
  }
  if (("task_group"%in%analysis[["predictor_list"]])&(!"task_group"%in%names(data))){
    task_groups <- read.table(paste(root_path,"original_data",task_group_file,sep="/"),header=T,stringsAsFactors = F)
    data        <- merge(data,task_groups,all.x=T,all.y=F)
  }
  plot_regression(data=data,time_point="comparison",analysis,n_cat_horiz=15,n_cat_vertic=30,pool=F,prepare=T,status=status,collective=collective,pool_plot=pool_plot,plot_untransformed=plot_untransformed,boldy=boldy,aligned=aligned)
}

plot_diff <- function(plot_dat,ytitle,col_vector,predictor,form_stat,p_colony,contrast.matrix,output,dataset,collective,violin_params=NULL,adjust_title_line){
  if (!is.null(violin_params)){
    violin_params <- as.numeric(unlist(violin_params))
    ##read violin param
    range <- violin_params[1]
    ylim_fac1 <- violin_params[2]
    ylim_fac2 <- violin_params[3]
    wex <- violin_params[4]
    h <- violin_params[5]
    
  }
  par_mgp_ori <- par()$mgp
  par_mar_ori <- par()$mar
  if (!((grepl("age",root_path))&collective)){
    means <- aggregate(na.rm=T,na.action="na.pass",variable~time+predictor,FUN="mean",data=plot_dat);ses <- aggregate(na.rm=T,na.action="na.pass",variable~time+predictor,FUN="std.error",data=plot_dat);
    names(means)[names(means)=="variable"] <- "mean";names(ses)[names(ses)=="variable"] <- "se";means <- merge(means,ses)
    means <- means[order(match(means$predictor,status_order),match(means$time,status_order)),]
    to_plot <- unique(means[c("time","predictor")])
  }else{
    means <- aggregate(na.rm=T,na.action="na.pass",variable~time,FUN="mean",data=plot_dat);ses <- aggregate(na.rm=T,na.action="na.pass",variable~time,FUN="std.error",data=plot_dat);
    names(means)[names(means)=="variable"] <- "mean";names(ses)[names(ses)=="variable"] <- "se";means <- merge(means,ses)
    means <- means[order(match(means$predictor,status_order)),]
    to_plot <- unique(means[c("time")])
  }
  if (grepl("\n",full_statuses_names[as.character(means$predictor)][1])){
    lassy <- 1
    mgpy <- par_mgp_ori+c(0,par_mgp_ori[1]-par_mgp_ori[2],0)
    fonty <- 3
    cexy <- min_cex
  }else{
    lassy <- 2
    mgpy <- c(0.8,0.05,0)
    fonty <- 1
    cexy <- inter_cex
  }
  means[is.na(means$se),"se"] <- 0
  ymin <- min(c(means$mean-means$se),na.rm=T);ymax<- max(c(means$mean+means$se),na.rm=T)
  ymin_ori <- ymin; ymax_ori <- ymax
  
  if (ymin>0){ymin <- 0};if (ymax<0){ymax <- 0}
  ####Now center on 0 roughly
  yminmax <- max(abs(ymin),abs(ymax))
  ymin <- -yminmax
  ymax <- yminmax
  
  ####Now get an idea of the spacing between ticks
  prospected_vals <- axisTicks(c(ymin,ymax),log=F)
  
  ####And correct ymin and yrange accordingly
  interval <- diff(prospected_vals)[1]
  while(min(prospected_vals)>ymin){
    prospected_vals <- c(min(prospected_vals)-interval,prospected_vals)
  }
  ymin <- min(prospected_vals)
  
  while(max(prospected_vals)<ymax){
    prospected_vals <- c(prospected_vals,max(prospected_vals)+interval)
  }
  ymax <- max(prospected_vals)
  
  if (ymax<interval){ymax <- interval}
  if (ymin>-interval){ymin <- -interval}
  
  ####and now center 0 roughly
  yminmax <- max(abs(ymin),abs(ymax))
  ymin <- -yminmax
  ymax <- yminmax
  
  
  
  if (grepl("point",plot_type)|grepl("boxplot",plot_type)|grepl("violinplot",plot_type)){
    rangy <- max(plot_dat$variable,na.rm=T)-min(plot_dat$variable,na.rm=T)
    ymin <- min(ymin, min(plot_dat$variable,na.rm=T)-0.1*rangy)
    ymax <- max(ymax, max(plot_dat$variable,na.rm=T)+0.1*rangy)
    yminmax <- max(abs(ymin),abs(ymax))
    ymin <- -yminmax
    ymax <- yminmax
    
  }
  
  
  # rangy <- (ymax-ymin)+(ymax-ymin)*abs(jitter(0,factor=5))
  # ymin <- ymin-0.05*rangy
  # ymax <- ymax+0.05*rangy
  #plotx <- unique(c(1:nrow(means))+sort(rep(1:(nrow(means)/2),2)-1))
  barwidth <- 0.5
  barspace <- 0.5
  arrowcodes <- c(1,2,3); names(arrowcodes) <- c("-1","1","0")
  if (!(collective&(grepl("age",root_path)))){
    plotx <- barplot(means$mean,plot=F,width=barwidth,space=barspace)  
    if (grepl("bars",plot_type)){
      ###The bars
      plotx <- barplot(means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)-barwidth,max(plotx)+barwidth),xlab="",ylab="",xaxt="n",bty="n",col=col_vector[as.character(means$predictor)],xaxs="i",yaxs="i",cex=0.5,cex.axis=min_cex,cex.lab=inter_cex,width=barwidth,space=barspace,xaxt="n",yaxt="n",xpd=F)  
      plot_arrows(means=means,plotx=plotx,plot_type=plot_type,LWD=line_max,LENGTH=0.05,colz=col_vector[as.character(means$predictor)])
      ###The points
      if (grepl("point",plot_type)){
        
        for (lvly in 1:nrow(means)){
          stripchart(plot_dat[which(plot_dat$time==means[lvly,"time"]&plot_dat$predictor==means[lvly,"predictor"]),"variable"],at=plotx[lvly],add=T,vertical=T,pch=16,col=alpha("black",0.5),method = 'jitter',ylim=c(ymin,ymax),cex=min_cex)
        }
        ymin_ori <- ymin
        ymax_ori <- ymax
        
      }
      
      vals <- axisTicks(c(ymin,ymax),log=F,nint=6);interval <- mean(diff(vals))
      while(min(vals)>ymin_ori){
        vals <- c(min(vals)-interval,vals)
      }
      while(max(vals)<ymax_ori){
        vals <- c(vals,max(vals)+interval)
      }
      labs <- formatC(vals)
      if ("0" %in%labs ){
        idx0 <- which(labs=="0")
        if (is.even(idx0)){
          labs[!is.even(1:length(labs))] <-""
        }else{
          labs[is.even(1:length(labs))] <-""
        }
      }
      axis(2,at=vals,labels=labs,cex.axis=min_cex)
      par(xpd=F)
      abline(h=0 ,col="black")
      
    }else if (grepl("boxplot",plot_type)){
      ####empty plot
      plot(plotx,means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)-barwidth,max(plotx)+barwidth),xlab="",ylab="",xaxt="n",bty="n",xaxs="i",yaxs="i",type="n",lwd=line_min,cex.axis=min_cex,cex.lab=inter_cex,xaxt="n",yaxt="n")  
      par(bty="n")
      for (lvly in 1:nrow(means)){
        boxplot(plot_dat[which(plot_dat$time==means[lvly,"time"]&plot_dat$predictor==means[lvly,"predictor"]),"variable"],at=plotx[lvly],add=T,range=1.5,notch=T,names=F,col=col_vector[as.character(means[lvly,"predictor"])],xlab="",ylab="",xaxt="n",xaxs="i",yaxs="i",cex.axis=min_cex,cex.lab=inter_cex,xaxt="n",yaxt="n",medlwd=line_max,boxlwd=line_min+0.5*(line_max-line_min),whisklwd=line_min,whisklty=1,staplelwd=line_min,boxwex=0.7,pch=16)
      }
      vals <- axisTicks(c(ymin,ymax),log=F,nint=6);interval <- mean(diff(vals))
      while(min(vals)>ymin){
        vals <- c(min(vals)-interval,vals)
      }
      while(max(vals)<ymax){
        vals <- c(vals,max(vals)+interval)
      }
      labs <- formatC(vals)
      if ("0" %in%labs ){
        idx0 <- which(labs=="0")
        if (is.even(idx0)){
          labs[!is.even(1:length(labs))] <-""
        }else{
          labs[is.even(1:length(labs))] <-""
        }
      }
      axis(2,at=vals,labels=labs,cex.axis=min_cex)
      par(xpd=F)
      abline(h=0 ,col="black")
      
      
    }else if (grepl("violinplot",plot_type)){
      ####empty plot
      plot(plotx,means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)-barwidth,max(plotx)+barwidth),xlab="",ylab="",xaxt="n",bty="n",xaxs="i",yaxs="i",type="n",lwd=line_min,cex.axis=min_cex,cex.lab=inter_cex,xaxt="n",yaxt="n")  
      par(bty="n")
      for (lvly in 1:nrow(means)){
        if (is.na(range)){
          VioPlot(na.omit(plot_dat[which(plot_dat$time==means[lvly,"time"]&plot_dat$predictor==means[lvly,"predictor"]),"variable"]),col=alpha(col_vector[as.character(means[lvly,"predictor"])],0.7), horizontal=F, at=plotx[lvly], add=TRUE,lty=1, rectCol="black",wex=wex,border=NA,lwd=line_min,mode="Mean")
        }else{
          VioPlot(na.omit(plot_dat[which(plot_dat$time==means[lvly,"time"]&plot_dat$predictor==means[lvly,"predictor"]),"variable"]),range=range, h=h,col=alpha(col_vector[as.character(means[lvly,"predictor"])],0.7), horizontal=F, at=plotx[lvly], add=TRUE,lty=1, rectCol="black",wex=wex,border=NA,lwd=line_min,mode="Mean")
        }
      }
      vals <- axisTicks(c(ymin,ymax),log=F,nint=6);interval <- mean(diff(vals))
      while(min(vals)>ymin){
        vals <- c(min(vals)-interval,vals)
      }
      while(max(vals)<ymax){
        vals <- c(vals,max(vals)+interval)
      }
      labs <- formatC(vals)
      if ("0" %in%labs ){
        idx0 <- which(labs=="0")
        if (is.even(idx0)){
          labs[!is.even(1:length(labs))] <-""
        }else{
          labs[is.even(1:length(labs))] <-""
        }
      }
      axis(2,at=vals,labels=labs,cex.axis=min_cex)
      par(xpd=F)
      abline(h=0 ,col="black")
      
      
    }else{
      ####empty plot
      plot(plotx,means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)-barwidth,max(plotx)+barwidth),xlab="",ylab="",xaxt="n",bty="n",xaxs="i",yaxs="i",type="n",lwd=line_min,cex.axis=min_cex,cex.lab=inter_cex,xaxt="n",yaxt="n")  
      ####arrows
      plot_arrows(means=means,plotx=plotx,plot_type=plot_type,LWD=line_max,LENGTH=0.05,colz=col_vector[as.character(means$predictor)])
      ####points
      points(plotx,means$mean,col=col_vector[as.character(means$predictor)],cex=1.5*max_cex,pch=16,lwd=line_min)  
      par(xpd=F)
      abline(h=0,lwd=line_max,col="black")
      vals <- axisTicks(c(ymin,ymax),log=F,nint=6);interval <- mean(diff(vals))
      while(min(vals)>ymin_ori){
        vals <- c(min(vals)-interval,vals)
      }
      while(max(vals)<ymax_ori){
        vals <- c(vals,max(vals)+interval)
      }
      labs <- formatC(vals)
      if ("0" %in%labs ){
        idx0 <- which(labs=="0")
        if (is.even(idx0)){
          labs[!is.even(1:length(labs))] <-""
        }else{
          labs[is.even(1:length(labs))] <-""
        }
      }
      
      axis(2,at=vals,labels=labs,lwd=0,lwd.ticks=line_inter,cex.axis=min_cex)
      abline(v=min(plotx)-barwidth,lwd=line_max,col="black")
    }
    # plot(plotx,means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)-1,max(plotx)+1),xlab="",ylab=ytitle,xaxt="n",bty="n",col=col_vector[as.character(means$predictor)],pch=19,xaxs="i",yaxs="i",cex=0.5,cex.axis=min_cex,cex.lab=inter_cex)
    # arrows (plotx,means$mean-means$se,plotx,means$mean+means$se,code=3,angle=90,col=col_vector[as.character(means$predictor)],lwd=line_max,length=0.15)
    par(xpd=T)
    axis(1,at=plotx,labels=full_statuses_names[as.character(means$predictor)],tick=F,cex.axis=cexy,las=lassy,mgp=mgpy,font=fonty)
    par(xpd=F)
  }else{
    plotx <- barplot(means$mean,plot=F,width=barwidth,space=barspace)  
    if (grepl("bars",plot_type)){
      plotx <- barplot(means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)-barwidth,max(plotx)+barwidth),xlab="",ylab="",bty="n",col=col_vector[as.character(means$predictor)],xaxs="i",yaxs="i",cex=0.5,cex.axis=min_cex,cex.lab=inter_cex,width=barwidth,space=barspace,xaxt="n",yaxt="n",xpd=F)  
      
      ###The points
      if (grepl("point",plot_type)){
        for (lvly in 1:nrow(means)){
          stripchart(plot_dat[which(plot_dat$time==means[lvly,"time"]),"variable"],at=plotx[lvly],add=T,vertical=T,pch=16,col=alpha("black",0.5),method = 'jitter',ylim=c(ymin,ymax),cex=min_cex)
        }
        ymin_ori <- ymin
        ymax_ori <- ymax
        
      }
      
      vals <- axisTicks(c(ymin,ymax),log=F,nint=6);interval <- mean(diff(vals))
      while(min(vals)>ymin_ori){
        vals <- c(min(vals)-interval,vals)
      }
      while(max(vals)<ymax_ori){
        vals <- c(vals,max(vals)+interval)
      }
      labs <- formatC(vals)
      if ("0" %in%labs ){
        idx0 <- which(labs=="0")
        if (is.even(idx0)){
          labs[!is.even(1:length(labs))] <-""
        }else{
          labs[is.even(1:length(labs))] <-""
        }
      }
      
      axis(2,at=vals,labels=labs,cex.axis=min_cex)
      par(xpd=F)
      abline(h=0 ,col="black")
      plot_arrows(means=means,plotx=plotx,plot_type=plot_type,LWD=line_max,LENGTH=0.05,colz=col_vector[as.character(means$predictor)])
      
    }else if (grepl("boxplot",plot_type)){
      ####empty plot
      plot(plotx,means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)-barwidth,max(plotx)+barwidth),xlab="",ylab="",xaxt="n",bty="n",xaxs="i",yaxs="i",type="n",lwd=line_min,cex.axis=min_cex,cex.lab=inter_cex,xaxt="n",yaxt="n")  
      par(bty="n")
      for (lvly in 1:nrow(means)){
        boxplot(plot_dat[which(plot_dat$time==means[lvly,"time"]),"variable"],at=plotx[lvly],add=T,range=1.5,notch=T,names=F,col=col_vector[as.character(means[lvly,"predictor"])],xlab="",ylab="",xaxt="n",xaxs="i",yaxs="i",cex.axis=min_cex,cex.lab=inter_cex,xaxt="n",yaxt="n",medlwd=line_max,boxlwd=line_min+0.5*(line_max-line_min),whisklwd=line_min,whisklty=1,staplelwd=line_min,boxwex=0.7,pch=16)
      }
      vals <- axisTicks(c(ymin,ymax),log=F,nint=6);interval <- mean(diff(vals))
      while(min(vals)>ymin){
        vals <- c(min(vals)-interval,vals)
      }
      while(max(vals)<ymax){
        vals <- c(vals,max(vals)+interval)
      }
      labs <- formatC(vals)
      if ("0" %in%labs ){
        idx0 <- which(labs=="0")
        if (is.even(idx0)){
          labs[!is.even(1:length(labs))] <-""
        }else{
          labs[is.even(1:length(labs))] <-""
        }
      }
      axis(2,at=vals,labels=labs,cex.axis=min_cex)
      par(xpd=F)
      abline(h=0 ,col="black")
      
    }else if (grepl("violinplot",plot_type)){
      ####empty plot
      plot(plotx,means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)-barwidth,max(plotx)+barwidth),xlab="",ylab="",xaxt="n",bty="n",xaxs="i",yaxs="i",type="n",lwd=line_min,cex.axis=min_cex,cex.lab=inter_cex,xaxt="n",yaxt="n")  
      par(bty="n")
      for (lvly in 1:nrow(means)){
        if (is.na(range)){
          VioPlot(na.omit(plot_dat[which(plot_dat$time==means[lvly,"time"]),"variable"]),col=alpha(col_vector[as.character(means[lvly,"predictor"])],0.7), horizontal=F, at=plotx[lvly], add=TRUE,lty=1, rectCol="black",wex=wex,border=NA,lwd=line_min,mode="Mean")
        }else{
          VioPlot(na.omit(plot_dat[which(plot_dat$time==means[lvly,"time"]),"variable"]),range=range, h=h,col=alpha(col_vector[as.character(means[lvly,"predictor"])],0.7), horizontal=F, at=plotx[lvly], add=TRUE,lty=1, rectCol="black",wex=wex,border=NA,lwd=line_min,mode="Mean")
        }
      }
      vals <- axisTicks(c(ymin,ymax),log=F,nint=6);interval <- mean(diff(vals))
      while(min(vals)>ymin){
        vals <- c(min(vals)-interval,vals)
      }
      while(max(vals)<ymax){
        vals <- c(vals,max(vals)+interval)
      }
      labs <- formatC(vals)
      if ("0" %in%labs ){
        idx0 <- which(labs=="0")
        if (is.even(idx0)){
          labs[!is.even(1:length(labs))] <-""
        }else{
          labs[is.even(1:length(labs))] <-""
        }
      }
      axis(2,at=vals,labels=labs,cex.axis=min_cex)
      par(xpd=F)
      abline(h=0 ,col="black")
      
    }else{
      ####empty plot
      plot(plotx,means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)-barwidth,max(plotx)+barwidth),xlab="",ylab="",bty="n",xaxs="i",yaxs="i",type="n",cex.axis=min_cex,cex.lab=inter_cex,lwd=line_min,xaxt="n",yaxt="n")  
      ####arrows
      plot_arrows(means=means,plotx=plotx,plot_type=plot_type,LWD=line_max,LENGTH=0.05,colz=col_vector[as.character(means$predictor)])
      ####points
      points(plotx,means$mean,col=col_vector[as.character(means$predictor)],cex=1.5*max_cex,pch=16,lwd=line_min)  
      par(xpd=F)
      abline(h=0,lwd=line_max,col="black")
      vals <- axisTicks(c(ymin,ymax),log=F,nint=6);interval <- mean(diff(vals))
      while(min(vals)>ymin_ori){
        vals <- c(min(vals)-interval,vals)
      }
      while(max(vals)<ymax_ori){
        vals <- c(vals,max(vals)+interval)
      }
      labs <- formatC(vals)
      if ("0" %in%labs ){
        idx0 <- which(labs=="0")
        if (is.even(idx0)){
          labs[!is.even(1:length(labs))] <-""
        }else{
          labs[is.even(1:length(labs))] <-""
        }
      }
      
      axis(2,at=vals,labels=labs,lwd=0,lwd.ticks=line_inter,cex.axis=min_cex)
      abline(v=min(plotx)-barwidth,lwd=line_max,col="black")
    }
    
    # plot(plotx,means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)-1,max(plotx)+1),xlab="",ylab=ytitle,bty="n",col=col_vector[as.character(means$time)],pch=19,xaxs="i",yaxs="i",cex=0.5,cex.axis=min_cex,cex.lab=inter_cex)
    # arrows (plotx,means$mean-means$se,plotx,means$mean+means$se,code=3,angle=90,col=col_vector[as.character(means$time)],lwd=line_max,length=0.15)
  }
  ###Now plot ylab
  if (grepl("Modularity",ytitle)|grepl("transmission",ytitle)|grepl("with brood",ytitle)){
    ytil <- list(quote("Pathogen-induced changes"),quote("relative to sham-induced changes"))
    par(xpd=NA)
    mtext(side=2,do.call(expression,ytil),line=pars$mgp[1]*(rev(1:length(ytil))-1)+pars$mgp[1],cex=par("cex") *inter_cex,xpd=NA)
    par(xpd=F)
  }
  
  pvalue <- output[["interaction_problist"]][["24"]][["1"]]
  # print(pvalue)
  if (pvalue>0.05){p_cex <- inter_cex;adjust_line <- 0.2;fonty <- 1}else{p_cex <- max_cex*1.1;adjust_line <- -0.1; fonty <-  2}
  title(main=from_p_to_ptext(pvalue),cex.main=p_cex,font.main=fonty,line=stat_line+adjust_line,xpd=T)
  par(xpd=NA)
  title(main=ytitle,cex.main=inter_cex,font.main=2,line=1+adjust_title_line,xpd=NA)
  par(xpd=F)
}

plot_distribution <- function(experiments,desired_treatments){
  par_mar_ori <- par()$mar
  par(mar=par_mar_ori+c(0,0,0,0.5))
  if (experiments=="all"){
    experiments <- c("age_experiment","survival_experiment","main_experiment")
  }
  if (experiments =="both"){
    experiments <- c("age_experiment","main_experiment")
  }
  transf <- function(x){
    return(x^(1/2))
  }
  rev_transf <- function(x){
    return(x^2)
  }
  xlabel <- substitute(sqrt ( xlabely),list(xlabely="Simulated load"))
  ####read data
  infection_data <- NULL
  for (experiment in experiments){
    setwd(paste(disk_path,experiment,"transmission_simulations/pre_vs_post_treatment/experimentally_exposed_seeds",sep="/"))
    si_outcome <- read.table("individual_simulation_results_observed.txt",header=T,stringsAsFactors = T)
    si_outcome["ant_id"] <- as.character(interaction(experiment,si_outcome$colony,si_outcome$tag))
    
    for (time_point in c("before","after")){
      si_outcome_temp <- si_outcome[which(si_outcome$period==time_point),c("colony","treatment","tag","status","antid","simulated_load","transmission_rank")]
      names(si_outcome_temp)[names(si_outcome_temp)%in%c("simulated_load","transmission_rank")] <- paste(names(si_outcome_temp)[names(si_outcome_temp)%in%c("simulated_load","transmission_rank")],time_point,sep="_")
      assign(paste("si_outcome_",time_point,sep=""),si_outcome_temp)
    }
    si_outcome <- merge(si_outcome_before,si_outcome_after)
    infection_data <- rbind(infection_data,data.frame(experiment=experiment,si_outcome,stringsAsFactors = F))
  }
  
  ####fill in missing untreated and queen info
  infection_data[which(infection_data$tag==queenid),"tag"] <- "queen"
  
  ###modify data
  infection_data[,"colony"] <- as.character(interaction(infection_data$experiment,infection_data$colony))
  
  ####make new datasets for further analyses######
  infection_data$status <- as.character(infection_data$status)
  infection_data <- infection_data[infection_data$status!="treated",]###contains queens and untreated workers
  infection_data <- infection_data[infection_data$tag!="queen",]###contains untreated workers
  
  #####1. Make bins
  xmin <- min(c(transf(infection_data$simulated_load_before),transf(infection_data$simulated_load_after)),na.rm=T)
  xmin_bis <- min(
    c(transf(infection_data$simulated_load_before),transf(infection_data$simulated_load_after))
    [
      c(transf(infection_data$simulated_load_before),transf(infection_data$simulated_load_after))
      >0
      ]
    ,
    na.rm=T
  )
  xmax <- max(c(transf(infection_data$simulated_load_before),transf(infection_data$simulated_load_after)),na.rm=T)
  xmin <- min(c(floor((xmin)*100)/100))
  xmax <- ceiling((xmax)*100)/100
  breaks <- seq(from=xmin,to=xmax,length.out=23)
  bw <- 0.09
  xmax <- ceiling((xmax)*10)/10
  
  ###2. plot histograms plus densities
  infection_data <- infection_data[which(infection_data$treatment%in%desired_treatments),]
  
  afters_dens <- density(transf(infection_data$simulated_load_after),from=xmin,to=xmax,n=500,na.rm=T,bw=bw)
  befores_dens <- density(transf(infection_data$simulated_load_before),from=xmin,to=xmax,n=500,na.rm=T,bw=bw)
  
  percolony_density <- NULL
  for (colony in sort(unique(infection_data$colony[!is.na(infection_data$colony)]))){
    if (colony!="age_experiment.colony021"){
      subsett <- infection_data[which(infection_data$colony==colony),]
      
      afters <- density(transf(subsett$simulated_load_after),from=xmin,to=xmax,n=500,na.rm=T,bw=bw)
      afters$y <- afters$y/sum(afters$y)
      befores <- density(transf(subsett$simulated_load_before),from=xmin,to=xmax,n=500,na.rm=T,bw=bw)
      befores$y <- befores$y/sum(befores$y)
      percolony_density <- rbind(percolony_density,data.frame(colony=colony,xcoor=afters$x,density_after=afters$y,density_before=befores$y,density_diff=afters$y-befores$y))
    }
  }
  
  ####plot actual distributions
  forplot <- data.frame(as.matrix(aggregate(cbind(density_after,density_before,density_diff)~xcoor,function(x)cbind(mean(x),std.error(x)),data=percolony_density)))
  names(forplot) <- c("x","mean_after","std.error_after","mean_before","std.error_before","mean","std.error")
  forplot["lower_y"] <- forplot$mean-forplot$std.error
  forplot["top_y"] <- forplot$mean+forplot$std.error
  forplot <- forplot[order(forplot$x),]
  
  xshade <- c(forplot$x,rev(forplot$x))
  yshade <- c(forplot$lower_y,rev(forplot$top_y))
  
  ##get ylims for plots ####
  ######first get an idea of how the density plots will be distributed
  ymin_dens <- 2*min(forplot$lower_y)
  ymax_dens <- max(c(forplot$mean_after,forplot$mean_before))+0.1*(max(c(forplot$mean_after,forplot$mean_before))-min(forplot$mean))
  neg <- abs(ymin_dens)/abs(ymax_dens)
  
  
  ######second get an idea of how the histogram will be distributed
  afters <- hist(transf(infection_data$simulated_load_after),breaks=breaks,plot=F)
  befores <- hist(transf(infection_data$simulated_load_before),breaks=breaks,plot=F)
  ymax <- max(c(afters$density,befores$density))+0.1*max(c(afters$density,befores$density))
  #####...and deduce ymin
  ymin <- -neg*ymax
  
  ###then plot histograms; in frequency
  befores <- hist(transf(infection_data$simulated_load_before),breaks=breaks,plot=T,col=alpha("blue",0.4),prob=T,ylim=c(ymin,ymax),xlim=c(min(0,xmin),xmax),xaxt="n",yaxt="n",xaxs="i",yaxs="i",main="",bty="l",cex.axis=min_cex,cex.lab=inter_cex,xlab="",lwd=line_min/10,border="black")
  afters <- hist(transf(infection_data$simulated_load_after),breaks=breaks,col=alpha("red",0.3),add=T,plot=T,prob=T,lwd=line_min/10,border="black")
  prospected_ats <- axisTicks(c(min(0,xmin),xmax),log=F)
  prospected_ats <- c(prospected_ats,prospected_ats[length(prospected_ats)]+prospected_ats[length(prospected_ats)]-prospected_ats[-1+length(prospected_ats)])
  axis(1,at=prospected_ats,cex.axis=min_cex,cex.lab=inter_cex)
  axis(2,cex.axis=min_cex,cex.lab=inter_cex)
  title(xlab=xlabel,cex.axis=min_cex,cex.lab=inter_cex,mgp=par()$mgp+c(0.2,0,0))
  ###then add densities; with new axes
  par(new=T)
  #####get new ymin, ymax
  plot(forplot$x,forplot$mean, type="n", axes=F, xlab=NA, ylab=NA, cex=1.2,col=statuses_colours[desired_treatments],ylim=c(ymin_dens,ymax_dens),xaxs="i",yaxs="i",main="",bty="l",cex.axis=min_cex,cex.lab=inter_cex,xlim=c(min(0,xmin),xmax))
  polygon(xshade,yshade, border = alpha("yellow",0),col=alpha("yellow",0.6))
  lines(forplot$x,forplot$mean,lwd=line_max,col="black")
  
  ####write legend
  legend(x=par("usr")[2]+0.05*(par("usr")[2]-par("usr")[1]),y=par("usr")[4]-0.025*(par("usr")[4]-par("usr")[3]),xjust=1,yjust=1,legend=c("Pre-treatment","Post-treatment","Difference"),pt.bg=c(alpha("blue",0.4),alpha("red",0.3),alpha("yellow",0.6)),col=c("black","black",alpha("yellow",0)),bty='n',pch=22,lty=0,lwd=0,pt.lwd=1,pt.cex=1.5,text.col="white",cex=min_cex)
  legend(x=par("usr")[2]+0.05*(par("usr")[2]-par("usr")[1]),y=par("usr")[4]-0.025*(par("usr")[4]-par("usr")[3]),xjust=1,yjust=1,legend=c("Pre-treatment","Post-treatment","Difference"),col=c(alpha("blue",0),alpha("red",0),"black"),bty='n',lty=c(0,0,1),lwd=c(0,0,1),cex=min_cex)
  
  
  print("KS-test:")
  ks_test <- ks.test(transf(infection_data$simulated_load_after),transf(infection_data$simulated_load_before))
  print(ks_test)
  p_value <- ks_test$p.value
  
  where_to_print_stat <- median(c(transf(infection_data$simulated_load_after),transf(infection_data$simulated_load_before)))
  
  par(xpd=T) 
  mtext(full_statuses_names[desired_treatments],side=3,line=stat_line,adj=0.5,cex=par("cex") *inter_cex,font=2)
  mtext(from_p_to_ptext(p_value),side=3,line=stat_line-1,adj=0.5,cex=par("cex") *max_cex,font=2,at=where_to_print_stat)
  
  # print("Thresholds at which after becomes lower than before")
  forplot["positive"] <- forplot$mean>=0
  change <- min(which(diff(forplot$positive)==-1))
  threshold1 <- rev_transf(forplot[change,"x"])
  threshold2 <- rev_transf(forplot[change+1,"x"])
  if(!exists("threshold")){threshold <-round(((threshold1+threshold2)/2)*10000)/10000}
  
  par(xpd=F)
  lines(x=c(transf(threshold),transf(threshold)),y=c(0,ymin_dens),col="springgreen3",lty=1,lwd=2*line_max)
  ###Now write down "high level", "low_level"
  arrows(x0=transf(threshold),y0=1.7*min(forplot$lower_y),x1=par("usr")[2],y1=1.7*min(forplot$lower_y),col="springgreen4",code=3,length=0.025)
  text(x=mean(c(transf(threshold),par("usr")[2])),y=1.5*min(forplot$lower_y),labels="high load",adj=c(0.5,0),cex=min_cex,col="springgreen4",font=3)
  xmin <-min(c(transf(infection_data$simulated_load_before),transf(infection_data$simulated_load_after)),na.rm=T)
  arrows(x1=transf(threshold),y0=1.7*min(forplot$lower_y),x0=xmin_bis,y1=1.7*min(forplot$lower_y),col="springgreen2",code=3,length=0.025)
  text(x=mean(c(transf(threshold),xmin)),y=1.5*min(forplot$lower_y),labels="low load",adj=c(0.5,0),cex=min_cex,col="springgreen2",font=3)
  
  par(mar=par_mar_ori)
}

plot_observed_vs_random <- function(experiments,variables,data_input=NULL,data_path,pattern){
 if (is.null(data_input)) {
    ###read-in data###
    data <- NULL
    
    for (experiment in experiments){
      print(experiment)
      ### data files
      setwd(paste(disk_path,experiment,data_path,sep="/"))
      file_list <- list.files(pattern=pattern)
      temp <- NULL
      for (file in file_list){
        dat <- read.table(file,header=T,stringsAsFactors=F)
        dat <- dat[,which(names(dat)%in%c("randy","colony","treatment","period","time_hours","time_of_day",variables))]
        temp <- rbind(temp,dat)
        rm(list=c("dat"))
      }
      temp <- temp[,order(names(temp))]
      temp <- data.frame(experiment=experiment,temp,stringsAsFactors = F)
      if (!is.null(data)){
        if (!all(names(data)%in%names(temp))){
          temp[names(data)[which(!names(data)%in%names(temp))]] <- NA
          temp <- temp[,names(data)]
        }
      }
      
      data <- rbind(data,temp)
      rm(list=c("temp"))
      
    }
    
    ####modify period values to be simple and match what scatterplot function expects
    data["colony"] <- as.character(interaction(data$experiment,data$colony))
    
  }else{data <- data_input}
  
   for (variable in variables){
    data["variable"] <- data[,variable]
    data[which(!is.finite(data$variable)),"variable"] <- NA
    
    #####get random and observed mean
    randys <- aggregate(variable~colony+treatment+randy,FUN=mean,data=data[which(data$randy!="observed"),]);
    randys <- as.data.frame(as.matrix(randys));randys$variable <- as.numeric(as.character(randys$variable))
    #####then get mean, median and standard error
    randys <- aggregate(variable~colony+treatment,function(x)cbind(mean(x),median(x),std.error(x),length(x)),data=randys);
    randys <- as.data.frame(as.matrix(randys))
    names(randys)[names(randys)=="variable.1"] <- "random_mean";names(randys)[names(randys)=="variable.2"] <- "random_median";names(randys)[names(randys)=="variable.3"] <- "random_std.error";names(randys)[names(randys)=="variable.4"] <- "random_nb"
    
    ###do the same for observed
    observeds <- aggregate(variable~colony+treatment+randy,FUN=mean,data=data[which(data$randy=="observed"),])
    observeds <- as.data.frame(as.matrix(observeds));observeds$variable <- as.numeric(as.character(observeds$variable))
    observeds <- aggregate(variable~colony+treatment,function(x)cbind(mean(x),median(x),length(x)),data=observeds);
    observeds <- as.data.frame(as.matrix(observeds))
    names(observeds)[names(observeds)=="variable.1"] <- "observed_mean";names(observeds)[names(observeds)=="variable.2"] <- "observed_median";names(observeds)[names(observeds)=="variable.3"] <- "observed_nb"
    
    randys <- merge(randys,observeds);
    randys$colony <- as.character(randys$colony)
    randys$observed_mean <- as.numeric(as.character(randys$observed_mean))
    randys$random_mean <- as.numeric(as.character( randys$random_mean))
    randys$observed_median <- as.numeric(as.character(randys$observed_median))
    randys$random_median <- as.numeric(as.character( randys$random_median))
    randys$random_std.error <- as.numeric(as.character( randys$random_std.error))
    
    randys["deviation"] <- randys$observed_median-randys$random_median
    randys["relative_deviation"] <- (randys$observed_median-randys$random_median)/abs(randys$random_median)
    randys["p_value"] <- NA
    randys["variable"] <- variable
    randys["n_observed"] <- NA
    randys["n_random"] <- NA
    
    for (coli in unique(randys$colony)){
      random <- aggregate(variable~colony+treatment+randy,FUN=mean,data=data[which(data$colony==coli&data$randy!="observed"),])[,"variable"];
      observed <- aggregate(variable~colony+treatment+randy,FUN=mean,data=data[which(data$colony==coli&data$randy=="observed"),])[,"variable"];
      #####Get one-sided p: proportion of random values that are greater than observed values. Will be 0 if all values are lower and 1 if all values are greater
      one_sided_p <- length(which(random>observed))/length(random)
      
      randys[which(randys$colony==coli),"p_value"] <- one_sided_p
      randys[which(randys$colony==coli),"one_sided_p_value_obs_lower_than_rand"] <- 1-one_sided_p
      randys[which(randys$colony==coli),"one_sided_p_value_obs_greater_than_rand"] <- one_sided_p
      randys[which(randys$colony==coli),"effect_sign"] <- sign( randys[which(randys$colony==coli),"deviation"])
      randys[which(randys$colony==coli),"n_observed"] <- length(observed)
      randys[which(randys$colony==coli),"n_random"] <- length(random)
    }
    randys_toprint <- data.frame(variable=variable,randys[c("colony","treatment","deviation","relative_deviation","effect_sign","one_sided_p_value_obs_lower_than_rand","one_sided_p_value_obs_greater_than_rand")],stringsAsFactors = F)
    randys_toprint[which(randys_toprint$effect_sign==1),"effect_signs"] <- "+"
    randys_toprint[which(randys_toprint$effect_sign==-1),"effect_signs"] <- "-"
    randys_toprint[which(randys_toprint$effect_sign==0),"effect_signs"] <- "0"
    
    #####now the stats: meta-analysis
    p_values_meta_analysis <- meta_analysis(p_values = randys[,"p_value"],effects = randys[,"relative_deviation"],std.errors = randys[,"random_std.error"])
    
    #####modify randys for plot
    forplot <- data.frame(network="random",randys[c("colony","treatment","random_median")],stringsAsFactors=F); names(forplot)[grepl("median",names(forplot))] <- "median"
    forplot2 <- data.frame(network="observed",randys[c("colony","treatment","observed_median")],stringsAsFactors=F); names(forplot2)[grepl("median",names(forplot2))] <- "median"
    
    forplot <- aggregate(median~network+colony+treatment,FUN=mean,data=forplot)
    forplot2 <- aggregate(median~network+colony+treatment,FUN=mean,data=forplot2)
    forplot <- rbind(forplot,forplot2)
    forplot$network <- factor(forplot$network,levels = c("random","observed"))
    ###get colony ordering from observed
    col_medians <- aggregate(median~colony,FUN=median,data=forplot2)
    col_medians <- col_medians [order(col_medians$median),]
    col_list <- col_medians$colony
    colour_pal <- colorRampPalette(brewer.pal(11, "Spectral"))(length(col_list))
    par(bty="n",xaxt = "n")
    
    for (coli in rev(col_list)){
      if (coli==col_list[length(col_list)]){
        addy <- F
      }else{
        addy <- T
      }
      
      titl <- names(variables[variables==variable])
      if (grepl("delta",titl)){
        titl1 <- unlist(strsplit(titl,"delta"))[1]
        titl2 <- unlist(strsplit(titl,"delta"))[2]
        titl <- substitute(paste(labo,Delta,laby),list(labo=titl1,laby=titl2))
      }
      stripchart(median ~ network, data = forplot[forplot$colony==coli,],vertical = TRUE,pch=16,col=alpha(colour_pal[which(coli==col_list)],1),method = 'jitter', jitter = 0.3,ylim=c(min(c(0,forplot$median)),max(c(forplot$median))), main = "",ylab=titl,add=addy,bty="l",cex=0.5,cex.axis=min_cex,cex.lab=inter_cex)
    }
    ###make boxplot
    forplot3 <- data.frame(as.matrix(aggregate(median~network,function(x)cbind(mean(x),std.error(x)),data=forplot)),stringsAsFactors = F)
    names(forplot3) <- c("network","mean","se")
    
    
    boxplot(median ~ network, data = forplot, 
            outline = FALSE, notch=F,    ## avoid double-plotting outliers, if any
            main = "",yaxt="n",add=T,col=alpha("white",0),medlwd=line_max,boxlwd=line_min+0.5*(line_max-line_min),whisklwd=line_min,whisklty=1,staplelwd=line_min,boxwex=0.7,bty="l")
    
    par(xpd=T)
    
    ###add stat
    pval <- p_values_meta_analysis$two_sided_p
    one_sidedpval <- p_values_meta_analysis$one_sided_p
    statistic <- p_values_meta_analysis$meta_statistic
    
    print(paste(variable,": z=",statistic,"; one-sided p =",one_sidedpval,"; two-sided p =",pval))
    
    if (pval>0.05){p_cex <- inter_cex;adjust_line <- 0.3;fonty <- 1}else{p_cex <- max_cex*1.1;adjust_line <- 0; fonty <-  2}
    
    title(main=from_p_to_ptext(pval),cex.main=p_cex,font.main=fonty,line=stat_line+adjust_line,xpd=T)
    par(xpd=F)
    par(xaxt = "s")
    axis(side=1,at=c(1,2),labels=c(full_statuses_names["random"],""),tick=F,lty=0,cex.axis=inter_cex)
    axis(side=1,at=c(1,2),labels=c("",full_statuses_names["observed"]),tick=F,lty=0,cex.axis=inter_cex)
    
  }
  
}

plot_network <- function(case,which_to_draw,colours="task_group"){
  ####Define two new vertex shapes (circle and square) which will allow to control the width of the vertex frame
  add.vertex.shape("fcircle", clip=igraph.shape.noclip,
                   plot=mycircle, parameters=list(vertex.frame.color=1,
                                                  vertex.frame.width=1))
  
  add.vertex.shape("fsquare", clip=igraph.shape.noclip,
                   plot=mysquare, parameters=list(vertex.frame.color=1,
                                                  vertex.frame.width=1))
  
  ##read qPCR data
  qPCR <- read.table(paste(root_path,"original_data/qPCR/qPCR_results.txt",sep="/"),header=T,stringsAsFactors = F)
  ##read simulation data
  si_outcome <- read.table(paste(root_path,"/transmission_simulations/pre_vs_post_treatment/experimentally_exposed_seeds/individual_simulation_results_observed.txt",sep=""),header=T,stringsAsFactors = F)
  ##read task_group data
  task_groups <- read.table(paste(root_path,"original_data",task_group_file,sep="/"),header=T,stringsAsFactors = F)
  ##read experimentally treated list
  treated <- read.table(paste(root_path,"original_data","treated_worker_list.txt",sep="/"),header=T,stringsAsFactors = F)
  
  setwd(paste(root_path,"example_networks_for_figures",case,sep="/"))
  networks <- list.files()
  for (current in which_to_draw){
    network_file <- networks[which( (grepl(unlist(strsplit(current,split="_"))[1],networks)) & (grepl(unlist(strsplit(current,split="_"))[2],networks)) )]
    interactions <- read.table(network_file,header=T,stringsAsFactors = F)
    
    if (current=="PreTreatment_random"){
      network_title <- "Random network"
    }else if (current=="PreTreatment_observed"&"PostTreatment_observed"%in%which_to_draw){
      network_title <- "Pre-treatment network"
    }else if (current=="PreTreatment_observed"&"PreTreatment_random"%in%which_to_draw){
      network_title <- "Observed network"
    }else if (current=="PostTreatment_observed"){
      network_title <- "Post-treatment network"
    }
    
    ###Make graph object
    actors <- data.frame(name=as.character(unique(c(interactions$Tag1,interactions$Tag2))))
    net <- graph.data.frame(interactions[c("Tag1","Tag2")],directed=F,vertices=actors)
    E(net)$weight <- interactions[,"duration_frames"]
    net <- simplify(net,remove.multiple=TRUE,remove.loop=TRUE,edge.attr.comb="sum")
    
    ###Drawing network
    WEIGHTS <- E(net)$weight
    ########first change the weight values so it ranges from mini to 0, then define edge widths and colors
    mini <- -100
    normal <- (WEIGHTS-min(WEIGHTS)) / max(WEIGHTS-min(WEIGHTS)) ####ranging from 0 to 1
    normal <- (normal-max(normal))*(mini/min(normal-max(normal)))
    normal <- exp(normal/(-mini/(1/50)))
    normal <- (normal-min(normal)) / max(normal-min(normal))
    E(net)$ori_width <-  2*normal
    E(net)$width <-  0.5*line_min+(3*line_max-0.5*line_min)*((E(net)$ori_width-min(E(net)$ori_width))/max(E(net)$ori_width-min(E(net)$ori_width)))
    E(net)$color <-  grey( 8/10*( 1-normal ))
    ########define node size, shapes, and labels
    V(net)$size  <- 10
    V(net)$shape <- "fcircle"
    V(net)$label=""
    ## Define layout
    E(net)$Scaled_Contact_Weights           <- sqrt(E(net)$weight)    / max(sqrt(E(net)$weight))
    net2 <- net                                                                 ## get layout from edge=weight thresholded graph
    Edge_Cutoff_P <- quantile(E(net2)$Scaled_Contact_Weights, probs=0.25)
    net2 <- igraph::delete.edges(net2, E(net2) [E(net2)$Scaled_Contact_Weights < Edge_Cutoff_P ] )  ## delete edges below a threshold - leaves only the strongest edges indicative of high spatial correlatoin
    # spring-embeded layout using contact weights
    lay <-  data.frame(tag=V(net)$name, layout_with_fr (net2, weights=E(net2)$Scaled_Contact_Weights, niter=50000)) ; colnames(lay )[2:3] <- c("x","y")
    ## contract outlier nodes
    PLOT_CONTRACTIONS <- F
    N_CONTRACTIONS <- 4
    OUTLIER_QUANTILE <- 0.9
    if (PLOT_CONTRACTIONS==TRUE) {par(mfrow=c(2,N_CONTRACTIONS/2))}
    if (N_CONTRACTIONS>0)
    {
      for (IT in 1:N_CONTRACTIONS)
      {
        ## centre so the CoG of all the nodes is at 0,0
        lay$x <- lay$x - mean(lay$x)
        lay$y <- lay$y - mean(lay$y)
        ## get nearest neighbour distance for each node & define those above a threshold as outliers
        lay$nnd <- rowMeans(nndist(X=lay$x, Y=lay$y, k=1:2))
        Outliers <- which(lay$nnd > quantile(lay$nnd,probs=OUTLIER_QUANTILE))
        lay$Outliers <- 0 ;  lay$Outliers [Outliers] <- 1
        CoG <- colMeans(lay[lay$Outliers==0,c("x","y")])
        lay$x_zero <- lay$x - CoG[1]
        lay$y_zero <- lay$y - CoG[2]
        lay$corrections <- 1
        lay$corrections[lay$Outliers==1] <- lay$nnd[lay$Outliers==1]
        lay$corrections[lay$Outliers==1] <- lay$corrections[lay$Outliers==1] / max(lay$nnd[lay$Outliers==1])
        lay$corrections[lay$Outliers==1] <- 1 - lay$corrections[lay$Outliers==1]
        lay$corrections[lay$Outliers==1] <- lay$corrections[lay$Outliers==1] + mean(sqrt(((lay$x-CoG["x"])^2+((lay$x-CoG["y"])^2))))
        lay$corrections[lay$Outliers==1] <- lay$corrections[lay$Outliers==1] / max(lay$corrections[lay$Outliers==1])
        lay$X  <- NA; lay$Y  <- NA
        lay[c("X","Y")]            <- lay$corrections * lay [,c("x","y")]
        ## plot each contraction
        if (PLOT_CONTRACTIONS==TRUE)
        {
          plot(lay[,c("X","Y")], pch=lay[,"Outliers"]+1) ; abline(v=0, lty=2) ; abline(h=0, lty=2)
          points(CoG[1], CoG[2], pch=21, bg=2, cex=2)
          segments(x0 = lay$x[lay$Outliers==1], y0 = lay$y[lay$Outliers==1], x1 = lay$X[lay$Outliers==1], y1 = lay$Y[lay$Outliers==1], col=2)
        }
        ## update the starting  coords
        lay[c("x","y")] <- lay[c("X","Y")]
      }
    }
    ## normalize the corrected x,y coords to -1, 1
    colnames(lay ) [2:3]<- c("X","Y")
    lay[,c("X","Y")] <- apply(lay[,c("X","Y")], 2, normalize_to_range, Min=-1, Max=1)
    lay <- as.matrix(lay[match(V(net)$name,lay$tag),c("X","Y")])
    rownames(lay) <- V(net)$name
    
    ##rotate network to facilitate comparisons, using two reference ants, 665 and 329
    ref_vec <- c(1,-1)
    obs_vec <- c(lay["329","X"]-lay["665","X"],lay["329","Y"]-lay["665","Y"])
    theta <-  atan2(obs_vec[2],obs_vec[1]) - atan2(ref_vec[2],ref_vec[1])
    ####then get rotation center
    center <- c(X=mean(lay["665","X"]),Y=mean(lay["665","Y"]))
    
    ####Now do the rotation. First: make "center" the origin
    lay[,"X"] <- lay[,"X"]-center["X"]
    lay[,"Y"] <- lay[,"Y"]-center["Y"]
    ####Second make new lay vector and calculate new coordinates
    lay_bis <- lay
    lay_bis[,"X"] <- lay[,"X"]*cos(-theta)-lay[,"Y"]*sin(-theta)
    lay_bis[,"Y"] <- lay[,"Y"]*cos(-theta)+lay[,"X"]*sin(-theta)
    ###Third retranslate
    lay[,"X"] <- lay[,"X"]+center["X"]
    lay[,"Y"] <- lay[,"Y"]+center["Y"]
    lay_bis[,"X"] <- lay_bis[,"X"]+center["X"]
    lay_bis[,"Y"] <- lay_bis[,"Y"]+center["Y"]
    lay_bis[,c("X","Y")] <- apply(lay_bis[,c("X","Y")], 2, normalize_to_range, Min=-1, Max=1)
    
    ###define colours
    colony         <- unique(interactions$colony)
    colony_number  <- as.numeric(gsub("colony","",colony))
    colony_treated <- treated[which(treated$colony==colony_number),"tag"]
    for (colour in colours){
      if (colour=="task_group"){
        colony_task_groups <- task_groups[which(task_groups$colony==colony),]
        colony_task_groups["color"] <- statuses_colours[colony_task_groups[,"task_group"]]
        V(net)$color <- colony_task_groups[match(V(net)$name,colony_task_groups$tag),"color"]
        if (length(colours)>1){
          network_title <- "Task group"
        }
      }else if (colour=="measured_load"){
        qPCR <- qPCR[which(qPCR$colony==colony),c("tag","status","measured_load_ng_per_uL")]
        qPCR[qPCR$tag=="queen","tag"] <- queenid
        qPCR <- qPCR[which(!is.na(as.numeric(qPCR$tag))),]
        if (0%in%qPCR$measured_load_ng_per_uL){
          replace_val <- min(qPCR$measured_load_ng_per_uL[qPCR$measured_load_ng_per_uL!=0])/sqrt(2)
        }else{
          replace_val <- 0
        }
        
        qPCR$log_measured <- log10(qPCR$measured_load_ng_per_uL+replace_val)
        mini_val <- min(qPCR[,"log_measured"])
        maxi_val <- max(qPCR[,"log_measured"])
        
        qPCR$log_measured_normalised <- 0 + ((qPCR$log_measured-mini_val)/(maxi_val-mini_val))*(1-0)
        
        palette <- rev(brewer.pal(9,"YlGnBu"))
        colour_palette <- colorRampPalette(palette)(1001)
        col_threshold <- "red"
        
        colour_range <- c (mini_val,maxi_val)
        qPCR[,"colour"] <- colour_palette[1+ceiling(qPCR[,"log_measured_normalised"]*1000)]
        V(net)$color <- qPCR[match(V(net)$name,qPCR$tag),"colour"]
        
        high_loadzes                                     <- qPCR[which(qPCR[,"measured_load_ng_per_uL"]>translated_high_threshold),"tag"]
        V(net)$shape                                     <- "fcircle"
        V(net)$shape[which(V(net)$name%in%high_loadzes)] <- "fsquare"
        
        network_title <- "qPCR-measured pathogen load"
      }else if (colour=="simulated_load"){
        si_outcome <- si_outcome[which(si_outcome$colony==colony&si_outcome$period=="after"),c("tag","simulated_load")]
        
        high_loadzes <- si_outcome[which(si_outcome$simulated_load>high_threshold),"tag"]
        V(net)$shape                                     <- "fcircle"
        V(net)$shape[which(V(net)$name%in%high_loadzes)] <- "fsquare"
        
        si_outcome$log_load <- log10(si_outcome$simulated_load)
        mini_val <- min(si_outcome[,"log_load"])
        maxi_val <- max(si_outcome[,"log_load"])
        
        si_outcome$log_load_normalised <- 0 + ((si_outcome$log_load-mini_val)/(maxi_val-mini_val))*(1-0)
        colour_palette <- inferno(1001, alpha = 1, begin = 0, end = 1, direction = 1)
        colour_range <- c (mini_val,maxi_val)
        si_outcome[,"colour"] <- colour_palette[1+ceiling(si_outcome[,"log_load_normalised"]*1000)]
        V(net)$color <- si_outcome[match(V(net)$name,si_outcome$tag),"colour"]
        network_title <- "Simulated pathogen load"
      }
      
      ###Finally, plot net
      plot_net <- net
      ###To improve visibility, delete thinnest edges before plotting
      plot_net <- igraph::delete.edges(plot_net, E(plot_net) [(E(net)$ori_width/2)<0.02 ] )
      
      par_mar_ori <- par()$mar
      if("measured_load" %in% colours){
        par(mar=c(1.2,0,0.2,0))
      }else{
        par(mar=c(0.2,0,1.2,0))
      }
      
      plot(plot_net, layout=lay_bis, edge.arrow.size=0.5, main="",rescale=F, vertex.frame.color="black",vertex.frame.width=line_min )#,"nodes")
      points(lay_bis[which(V(plot_net)$name%in%colony_treated),"X"],lay_bis[which(V(plot_net)$name%in%colony_treated),"Y"],pch=16,col="black",cex=0.7)
      if("measured_load" %in% colours){
        title(sub=network_title,font.sub=2,cex.sub=inter_cex,line=0.2,xpd=NA)
      }else{
        title(main=network_title,font.main=2,cex.main=inter_cex,line=0.2,xpd=NA)
      }
      
      if (grepl("load",colour)){
        par(mar=c(2,1,1,1))
        par(mgp=c(1,0.1,0))
        if (colour=="measured_load"){
          image(y=1, x=1:length(colour_palette), matrix(1:length(colour_palette), ncol =1),cex.lab=min_cex, col= colour_palette, xlab="", ylab="", xaxt="n", yaxt="n", main="",font=2); box(lwd=line_inter)
          title(main=expression(paste("(ng/", mu, "L)")),cex.main=inter_cex)
        }else{
          image(y=1, x=1:length(colour_palette), matrix(1:length(colour_palette), ncol =1),cex.lab=min_cex, col= colour_palette, xlab="", ylab="", xaxt="n", yaxt="n", main="",font=2); box(lwd=line_inter)
          title(main="(Proportion of exposure dose)",cex.main=inter_cex,font.main=1)
        }
        colour_vals <- mini_val+c(0:1000)*(maxi_val-mini_val)/1000
        ats <- ceiling(colour_range[1]):floor(colour_range[2])
        
        labs <- format(10^(ats),scientific=T)
        for (at in 1:length(ats)){
          ats[at] <- closest_match(ats[at] ,colour_vals)
        }
        axis(side = 1, at=ats, labels = labs, las=1, cex.axis=min_cex,lend="square",tcl=-0.2,lwd=line_inter)  ## ensure to back-transform the labels!!
        if (colour=="simulated_load"){
          thresh <- closest_match(log10(high_threshold) ,colour_vals)
          
          par(xpd=F)
          abline(v=thresh,col="white",lwd=1.5)
          abline(v=thresh,col="springgreen2",lwd=1)
          
        }else{
          thresh <- closest_match(log10(translated_high_threshold) ,colour_vals)
          
          par(xpd=F)
          abline(v=thresh,col="white",lwd=1.5)
          abline(v=thresh,col=col_threshold,lwd=1)
        }
      }
      par(mar=par_mar_ori)
    }
  }
}

plot_qpcr <- function(experiments){
  plot1 <- function(data,experiment,ylim=ylim){
    par_mari <- par("mar")
    par(mar=par_mari-c(0,0.35,0,par_mari[4]))
    
    ###Plot (load)=f(caste)
    data["variable"] <- data$measured_load_ng_per_uL;data["predictor"] <- data$task_group
    replac_val <- min(data$variable[data$variable!=0],na.rm=T)/sqrt(2)
    
    ####transform 
    data <- aggregate(na.action=na.pass,na.rm=T,log10(variable+replac_val)~predictor+colony,FUN=mean,data=data)
    names(data)[grepl("variable",names(data))] <-"variable"
    
    means <- aggregate(na.rm=T,na.action="na.pass",variable~predictor,FUN="mean",data=data);ses <- aggregate(na.rm=T,na.action="na.pass",variable~predictor,FUN="std.error",data=data);
    
    names(means)[grepl("variable",names(means))] <- "mean";names(ses)[grepl("variable",names(ses))] <- "se";means <- merge(means,ses)
    means <- means[order(match(means$predictor,status_order)),]
    
    means[is.na(means$se),"se"] <- 0
    if (is.null(ylim)){
      ymin <- min(c(means$mean-means$se),na.rm=T);ymax<- max(c(means$mean+means$se),na.rm=T)
      ymin <- floor(ymin)
      ymax <- ceiling(ymax)
      rangy <- ymax-ymin
      ymax <- ymax+0.1*rangy
      yrange <- c(ymin,ymax)
    }else{
      ymin <- ylim[1]
      ymax <- ylim[2]
      yrange <- ymax-ymin
      yrange <- c(ymin-0.04*yrange,ymax+0.04*yrange)
    }
    
    barwidth <- 0.5; barwidth_fac_within <- 0.5; barwidth_fac_between <- 2
    barspace <- c(barwidth_fac_between,barwidth_fac_within,barwidth_fac_within)
    
    plotx <- barplot(means$mean,plot=F,width=barwidth,space=barspace)
    ####empty plot
    plot(plotx,means$mean,ylim=yrange,xlim=c(min(plotx)- 0.6*(barwidth_fac_between*barwidth/2+barwidth/2),max(plotx)+ 0.6*(barwidth_fac_between*barwidth/2+barwidth/2)),xlab="",ylab=expression(paste("Mean measured pathogen load (ng/", mu, "L)")),bty="n",xaxs="i",yaxs="i",type="n",cex.axis=min_cex,cex.lab=inter_cex,lwd=line_min,xaxt="n",yaxt="n")
    ####arrows
    plot_arrows(means=means,plotx=plotx,plot_type="means",LWD=line_max,LENGTH=0.025,colz=statuses_colours[as.character(means$predictor)])
    ####points
    points(plotx,means$mean,col=statuses_colours[as.character(means$predictor)],pch=16,cex=max_cex*0.8,lwd=line_min)
    ####Y-axis
    axis(2,at=ymin:floor(ymax),labels=format(10^(ymin:floor(ymax)),scientific=T),cex.axis=min_cex,cex.lab=inter_cex,lwd=0,lwd.ticks=1)
    
    par(xpd=T)
    labses <- full_statuses_names[as.character(means$predictor)]
    labses <- c(labses[1],rep("",length(labses)-1))
    axis(1,at=plotx,labels=labses,tick=F,cex.axis=inter_cex,las=1,mgp=c(0.8,0.8,0))
    
    for (labse in 2:length(labses)){
      print(par("mgp")[2] + (labse-1)*1)
      if (labse/2==round(labse/2)){
        mtext(full_statuses_names[as.character(means$predictor)][labse],side=1,at=plotx[labse],line=0+par("mgp")[1] + 1,cex=par("cex")*inter_cex)
      }else{
        mtext(full_statuses_names[as.character(means$predictor)][labse],side=1,at=plotx[labse],line=0+par("mgp")[1] ,cex=par("cex")*inter_cex)
      }
    }
    
    segments(x0=min(plotx)- 0.6*(barwidth_fac_between*barwidth/2+barwidth/2),y0=yrange[1],x1=max(plotx)+ 0.6*(barwidth_fac_between*barwidth/2+barwidth/2),lwd=line_max)
    segments(x0=min(plotx)- 0.6*(barwidth_fac_between*barwidth/2+barwidth/2),y0=yrange[1],y1=yrange[2],lwd=line_max)
    par(xpd=F)
    
    
    data$predictor <- factor(data$predictor,levels=means$predictor)
    
    model_lmer <- lmer(variable~predictor+(1|colony),data=data)
    test_norm(residuals(model_lmer))
    pvalue <- Anova(model_lmer)["predictor","Pr(>Chisq)"]
    print(Anova(model_lmer))
    # title(ylab=lab_title,cex.lab=1.5,mgp=c(5,1,0))
    par(xpd=F)
    abline(h=0)
    if (pvalue>0.05){p_cex <- inter_cex}else{p_cex <- max_cex}
    title(main=from_p_to_ptext(pvalue),cex.main=p_cex,font.main=2,line=stat_line,xpd=T)
    
    post_hoc <- summary(glht(model_lmer,linfct = mcp(predictor="Tukey")),test=adjusted("BH"))
    print(post_hoc)
    post_hoc_groups <- cld(post_hoc)$mcletters$Letters
    for (idx in 1:length(post_hoc_groups)){
      group <- post_hoc_groups[as.character(means[idx,"predictor"])]
      text(x=plotx[idx],y=ymax-0.1*(ymax-ymin),adj=c(0.5,0),labels=as.character(group),xpd=T,cex=inter_cex)
    }
    par(mar=par_mari)
  }
  
  all_qpcr_data <- NULL
  data_for_plot <- NULL
  for (experiment in experiments){
    ###read qpcr data
    data <- read.table(paste(disk_path,"/",experiment,"/original_data/qPCR/qPCR_results.txt",sep=""),header=T,stringsAsFactors = F)
    ###keep only ants
    data <- data[which(!is.na(as.numeric(data$tag))),]
    ###remove workers that died before the end
    data <- data[which(data$alive_at_sampling_time),]
    
    ###read task group
    task_groups <-  read.table(paste(disk_path,"/",experiment,"/original_data/",task_group_file,sep=""),header=T,stringsAsFactors = F)
    ### add task groups info to data
    data <- merge(data,task_groups,all.x=T,all.y=F)
    
    ###keep only pathogen treatments
    data <- data[which(data$treatment%in%c("pathogen")),]
    
    ###add metadata
    data <- data.frame(experiment=experiment,data,stringsAsFactors = F)
    data$period <- "after"
    ####list desired variables and transformations
    data$antid <- as.character(data$colony,data$tag)
    all_qpcr_data <- rbind(all_qpcr_data,data)
    data$colony <- as.character(interaction(data$experiment,data$colony))
    
    ###remove treated workers
    data <- data[which(data$status!="treated"),]
    data_for_plot <- rbind(data_for_plot,data)
  }
  
  all_sim_data <- NULL
  for (experiment in experiments){
    data <-read.table(paste(disk_path,experiment,"transmission_simulations/calibration","individual_simulation_results.dat",sep="/"),header=T,stringsAsFactors=F)
    all_sim_data <- rbind(all_sim_data,data.frame(experiment=experiment,data[names(data)!="age"],stringsAsFactors = F))
  }
  
  ###Now join qPCR and simulation data into single data frame
  all_qpcr_data <- all_qpcr_data[which(!is.na(as.numeric(all_qpcr_data$tag))),]
  all_qpcr_data <- all_qpcr_data[which(all_qpcr_data$alive_at_sampling_time),]
  
  all_sim_qpcr_data <- merge(all_qpcr_data[c("experiment","colony","treatment","tag","status","task_group","age","measured_load_ng_per_uL")],all_sim_data[c("experiment","colony","tag","simulated_load")])
  
  ###remove treated individuals
  all_sim_qpcr_data <- all_sim_qpcr_data[which(all_sim_qpcr_data$status!="treated"),]
  all_sim_qpcr_data["antid"] <- as.character(interaction(all_sim_qpcr_data$experiment,all_sim_qpcr_data$colony,all_sim_qpcr_data$tag))
  
  varb <- "simulated_load"
  variable_list <- c("measured_load_ng_per_uL")
  names(variable_list) <-  c("Measured pathogen load")
  predictor_list <- c( "simulated_load")
  names(predictor_list) <- c("Simulated pathogen load")
  transf_variable_list <- c("log")
  transf_predictor_list <- c("log")
  
  ymin <- floor(log10(min(all_sim_qpcr_data$measured_load_ng_per_uL[all_sim_qpcr_data$measured_load_ng_per_uL!=0])/sqrt(2)))
  ymax <- ceiling(log10(max(all_sim_qpcr_data$measured_load_ng_per_uL)))
  xmin <- floor(log10(min(all_sim_qpcr_data[all_sim_qpcr_data[,varb]!=0,varb])/sqrt(2)))
  xmax <- ceiling(log10(max(all_sim_qpcr_data[,varb])))
  
  analysis <- list(variable_list=variable_list,
                   predictor_list=predictor_list,
                   transf_variable_list=transf_variable_list,
                   transf_predictor_list=transf_predictor_list,
                   violin_plot_param = list(c(1,0,-0.025,0.2,0.2)))
  predicted_value <- plot_regression(data=all_sim_qpcr_data,time_point="after",analysis=analysis,n_cat_horiz=20,n_cat_vertic=11,pool=c(F,F),collective=T,input_color=colour_palette_age,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,point_cex=1.5,predict=high_threshold)
  par(xpd=NA)
  if (varb =="simulated_load"){
    title(sub=expression(italic("(Prop. exposure dose)")),cex.sub=min_cex,font.sub=1,mgp=c(1,0.1,0))
  }else if (varb =="predicted_measured_load_ng_per_uL_SI"){
    title(sub=expression(paste("(ng/", mu, "L)")),cex.sub=min_cex,font.sub=1,mgp=c(1,0.1,0))
  }
  par(xpd=F)
  
  if (!exists("ymin")){
    ylim <- NULL
  }else{
    ylim <- c(ymin,ymax)
  }
  plot1(data=data_for_plot,experiment="both",ylim=NULL)
  if (exists("predicted_value")){return(predicted_value)}else{return(NULL)} 
}

plot_qpcr_vs_distance_to_treated <- function(experiments){
  all_qpcr_data    <- NULL
  for (experiment in experiments){
    ###read qpcr data
    data <- read.table(paste(disk_path,"/",experiment,"/original_data/qPCR/qPCR_results.txt",sep=""),header=T,stringsAsFactors = F)
    ###keep only ants
    data <- data[which(!is.na(as.numeric(data$tag))),]
    ###remove workers that died before the end
    data <- data[which(data$alive_at_sampling_time),]
    
    ###read task group
    task_groups <-  read.table(paste(disk_path,"/",experiment,"/original_data/",task_group_file,sep=""),header=T,stringsAsFactors = F)
    ### add task groups info to data
    data <- merge(data,task_groups,all.x=T,all.y=F)
    
    ###keep only pathogen treatments
    data <- data[which(data$treatment%in%c("pathogen")),]
    
    ###add metadata
    data <- data.frame(experiment=experiment,data,stringsAsFactors = F)
    data$period <- "after"
    ####list desired variables and transformations
    data$antid <- as.character(data$colony,data$tag)
    all_qpcr_data <- rbind(all_qpcr_data,data)
  }
  
  all_network_data <- NULL
  for (experiment in experiments){
    ###read data
    data <-read.table(paste(disk_path,experiment,"processed_data/network_properties/post_treatment/node_properties_observed.txt",sep="/"),header=T,stringsAsFactors=F)
    ###read task group
    task_group <- read.table(paste(disk_path,experiment,"original_data",task_group_file,sep="/"),header=T,stringsAsFactors = F)
    data <- merge(data,task_group,all.x=T,all.y=F)
    ###keep only pathogen treatments
    data <- data[which(data$treatment%in%c("pathogen")),]
    
    all_network_data <- rbind(all_network_data,data.frame(experiment=experiment,data[names(data)!="age"],stringsAsFactors = F))
  }
  all_network_data <- aggregate((mean_aggregated_distance_to_treated)~experiment+colony+colony_size+treatment+tag+status+task_group,FUN=mean,data=all_network_data)
  names(all_network_data)[grepl("aggregated",names(all_network_data))] <- "mean_aggregated_distance_to_treated"
  
  all_sim_data <- NULL
  for (experiment in experiments){
    data <-read.table(paste(disk_path,experiment,"transmission_simulations/calibration","individual_simulation_results.dat",sep="/"),header=T,stringsAsFactors=F)
    ###keep only pathogen treatments
    data <- data[which(data$treatment%in%c("pathogen")),]
    all_sim_data <- rbind(all_sim_data,data.frame(experiment=experiment,data[names(data)!="age"],stringsAsFactors = F))
  }
  
  all_interaction_data <- NULL
  for (experiment in experiments){
    data <-read.table(paste(disk_path,experiment,"processed_data/individual_behaviour/post_treatment/interactions_with_treated.txt",sep="/"),header=T,stringsAsFactors=F)
    ###keep only pathogen treatments
    data <- data[which(data$treatment%in%c("pathogen")),]
    ###keep only post-treatment values
    data <- data[which(data$time_hours>=0),]
    ###sum interactions across entire post-treatment period
    data$duration_of_contact_with_treated_min <- as.numeric(data$duration_of_contact_with_treated_min)
    data <- aggregate(na.rm=T,na.action="na.pass",duration_of_contact_with_treated_min~colony+colony_size+treatment+tag+period,FUN=sum,data=data)
    all_interaction_data <- rbind(all_interaction_data,data.frame(experiment=experiment,data[names(data)!="age"],stringsAsFactors = F))
  }
  
  ###Join all data into single data frame
  all_qpcr_data <- all_qpcr_data[which(!is.na(as.numeric(all_qpcr_data$tag))),]
  all_qpcr_data <- all_qpcr_data[which(all_qpcr_data$alive_at_sampling_time),]
  
  overall_data  <- merge(all_network_data[c("experiment","colony","colony_size","tag","status","task_group","mean_aggregated_distance_to_treated")],all_qpcr_data[c("experiment","colony","tag","measured_load_ng_per_uL")],all.x=T)
  overall_data  <- merge(overall_data,all_sim_data[c("experiment","colony","tag","simulated_load")],all.x=T)
  overall_data  <- merge(overall_data,all_interaction_data[c("experiment","colony","tag","duration_of_contact_with_treated_min")],all.x=T)
  
  ###remove treated individuals
  overall_data <- overall_data[which(overall_data$status!="treated"),]
  overall_data["antid"] <- as.character(interaction(overall_data$experiment,overall_data$colony,overall_data$tag))
  
  ####define high, medium and low contact with treated
  overall_data["contact_with_treated"] <- NA
  overall_data[
    which(overall_data$duration_of_contact_with_treated_min<quantile(overall_data$duration_of_contact_with_treated_min,probs=1/3,na.rm=T))    
    ,
    "contact_with_treated"] <- "low"
  
  overall_data[
    which(overall_data$duration_of_contact_with_treated_min>=quantile(overall_data$duration_of_contact_with_treated_min,probs=2/3,na.rm=T))    
    ,
    "contact_with_treated"] <- "high"

  overall_data[
    which(
      overall_data$duration_of_contact_with_treated_min>=quantile(overall_data$duration_of_contact_with_treated_min,probs=1/3,na.rm=T)
      &
        overall_data$duration_of_contact_with_treated_min<quantile(overall_data$duration_of_contact_with_treated_min,probs=2/3,na.rm=T)
    )    
    ,
    "contact_with_treated"] <- "medium"
 
  overall_data$colony <- as.character(interaction(overall_data$experiment))
  ###Now perform analysis for all workers, nurses, and low contact with treated
  for (who in c("all","nurses","low_contact_with_treated")){
    if (who=="all"){subdata <- overall_data;NCAT <- 20}else if(who=="nurses"){subdata <- overall_data[which(overall_data$task_group=="nurse"),]; NCAT <- 6}else if(who=="low_contact_with_treated"){subdata <- overall_data[which(overall_data$contact_with_treated=="low"),]; NCAT <- 20}
    
    variable_list <- c("measured_load_ng_per_uL","simulated_load")
    names(variable_list) <-  c("Measured pathogen load","Simulated pathogen load")
    predictor_list <- c("mean_aggregated_distance_to_treated","mean_aggregated_distance_to_treated")
    names(predictor_list) <- c("Mean network distance to treated","Mean network distance to treated")
    transf_variable_list <- c("log","sqrt")
    transf_predictor_list <- c("none","none")
    analysis <- list(variable_list=variable_list,
                     predictor_list=predictor_list,
                     transf_variable_list=transf_variable_list,
                     transf_predictor_list=transf_predictor_list,
                     violin_plot_param = list(c(1,0,-0.025,4,0.2),c(1,0,-0.025,4,0.02)))
    
    xmin_temp <- floor(round(10*(min(overall_data[overall_data$mean_aggregated_distance_to_treated!=0,"mean_aggregated_distance_to_treated"])/sqrt(2)))/10)
    xmax_temp <- ceiling(round(10*(max(overall_data$mean_aggregated_distance_to_treated)))/10)
    
    plot_regression(data=subdata,time_point="after",analysis=analysis,n_cat_horiz=NCAT,n_cat_vertic=11,pool=c(F,F),collective=T,input_color=colour_palette_age,point_cex=1.5,predict=NULL,xmin=xmin_temp,xmax=xmax_temp)
  }
}

plot_refined <- function(diff,lab_title,col_vector,predictor,output,contrast.matrix,survival,dataset,p_queen=NA,violin_params=NULL,aligned=F,adjust_title_line){
  if (!is.null(violin_params)){
    violin_params <- as.numeric(unlist(violin_params))
    ##read violin param
    range <- violin_params[1]
    ylim_fac1 <- violin_params[2]
    ylim_fac2 <- violin_params[3]
    wex <- violin_params[4]
    h <- violin_params[5]
    
  }
  par_mgp_ori <- par()$mgp
  par_mar_ori <- par()$mar
  means <- aggregate(na.rm=T,na.action="na.pass",variable~time+treatment+predictor,FUN="mean",data=diff);ses <- aggregate(na.rm=T,na.action="na.pass",variable~time+treatment+predictor,FUN="std.error",data=diff);
  names(means)[names(means)=="variable"] <- "mean";names(ses)[names(ses)=="variable"] <- "se";means <- merge(means,ses)
  means <- means[order(match(means$predictor,levels(diff$predictor)),match(means$treatment,levels(diff$treatment))),]
  to_plot <- unique(means[c("time","treatment","predictor")])
  #plotx <- c(1:nrow(means))+sort(rep(1:length(unique(means$predictor)),(nrow(means)/length(unique(means$predictor))))-1)
  
  means[is.na(means$se),"se"] <- 0
  ymin <- min(c(means$mean-means$se),na.rm=T);ymax<- max(c(means$mean+means$se),na.rm=T)
  ymin_ori <- ymin; ymax_ori <- ymax
  
  if (ymin>0){ymin <- 0};if (ymax<0){ymax <- 0}
  
  # ####Now center on 0 roughly
  if (Extended|aligned){
    yminmax <- max(abs(ymin),abs(ymax))
    ymin <- -yminmax
    ymax <- yminmax
  }
  
  ####Now get an idea of the spacing between ticks
  prospected_vals <- axisTicks(c(ymin,ymax),log=F)
  
  ####And correct ymin and yrange accordingly
  interval <- diff(prospected_vals)[1]
  while(min(prospected_vals)>ymin){
    prospected_vals <- c(min(prospected_vals)-interval,prospected_vals)
  }
  ymin <- min(prospected_vals)
  
  while(max(prospected_vals)<ymax){
    prospected_vals <- c(prospected_vals,max(prospected_vals)+interval)
  }
  ymax <- max(prospected_vals)
  
  if (ymax<interval){ymax <- interval}
  if (ymin>-interval){ymin <- -interval}
  
  ####and now center 0 roughly
  if (Extended|aligned){
    yminmax <- max(abs(ymin),abs(ymax))
    ymin <- -yminmax
    ymax <- yminmax
  }
  
  
  
  if (grepl("point",plot_type)|grepl("boxplot",plot_type)|grepl("violinplot",plot_type)){
    rangy <- max(diff$variable,na.rm=T)-min(diff$variable,na.rm=T)
    ymin <- min(ymin, min(diff$variable,na.rm=T)-0.1*rangy)
    ymax <- max(ymax, max(diff$variable,na.rm=T)+0.1*rangy)
    yminmax <- max(abs(ymin),abs(ymax))
    ymin <- -yminmax
    ymax <- yminmax
    
  }
  
  # rangy <- (ymax-ymin)
  # ymax <- ymax+0.3*rangy
  # 
  # plot(plotx,means$mean,ylim=c(ymin-0.3*(ymax-ymin),ymax+0.3*(ymax-ymin)),xlim=c(min(plotx)-1,max(plotx)+1),xlab="",ylab="",xaxt="n",bty="n",col=col_vector[as.character(means$treatment)],pch=19,xaxs="i",yaxs="i",cex.axis=1.1,yaxt="n")
  # arrows (plotx,means$mean-means$se,plotx,means$mean+means$se,code=3,angle=90,col=col_vector[as.character(means$treatment)],lwd=2,length=0.15)
  # axis(1,at=plotx,labels=treatment_names[as.character(means$treatment)],tick=F,cex.axis=1.54)
  barwidth <- 0.5; barwidth_fac_within <- 0.5; barwidth_fac_between <- 2
  barspace <- rep(c(barwidth_fac_between,barwidth_fac_within),length(means$mean)/2)
  
  plotx <- barplot(means$mean,plot=F,width=barwidth,space=barspace)
  ###prepare colours
  means["alpha"] <- as.numeric(alphas[as.character(means$treatment)])
  means["full_col"] <- col_vector[as.character(means$predictor)]
  means[which(as.character(means$predictor)=="outdoor_ant"),"alpha"] <- means[which(as.character(means$predictor)=="outdoor_ant"),"alpha"]*10
  # for (colidx in 1:nrow(means)){
  #   means[colidx,"final_col"] <- apply_alpha(means[colidx,"full_col"],means[colidx,"alpha"])
  # }
  means["final_col"] <- means["full_col"]
  if (grepl("bars",plot_type)){
    offset <- 0
    direction <- "normal"
    plotx <- barplot(means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)- (barwidth_fac_between*barwidth/2+barwidth/2),max(plotx)+ (barwidth_fac_between*barwidth/2+barwidth/2)),xlab="",ylab="",bty="n",col=means$final_col,xaxs="i",yaxs="i",cex=0.5,cex.axis=min_cex,cex.lab=inter_cex,width=barwidth,space=barspace,xaxt="n",yaxt="n",xpd=F,offset=offset)  
    
    
    if (grepl("point",plot_type)){
      
      for (lvly in 1:nrow(means)){
        stripchart(diff[which(diff$time==means[lvly,"time"]&diff$treatment==means[lvly,"treatment"]&diff$predictor==means[lvly,"predictor"]),"variable"],at=plotx[lvly],add=T,vertical=T,pch=16,col=alpha("black",0.5),method = 'jitter',ylim=c(ymin,ymax),cex=min_cex)
      }
      ymin_ori <- ymin
      ymax_ori <- ymax
    }
    
    vals <- axisTicks(c(ymin,ymax),log=F,nint=6);interval <- mean(diff(vals))
    while(min(vals)>ymin_ori){
      vals <- c(min(vals)-interval,vals)
    }
    while(max(vals)<ymax_ori){
      vals <- c(vals,max(vals)+interval)
    }
    labs <- formatC(vals)
    #print(labs)
    if ("0" %in%labs ){
      idx0 <- which(labs=="0")
      if (is.even(idx0)){
        labs[!is.even(1:length(labs))] <-""
      }else{
        labs[is.even(1:length(labs))] <-""
      }
    }
    #print(labs)
    
    axis(2,at=vals,labels=labs,cex.axis=min_cex)
    
    means$mean <- means$mean + offset
    plot_arrows(means=means,plotx=plotx,plot_type=plot_type,LWD=line_max,LENGTH=0.025,colz=means$final_col,direction=direction)
    par(xpd=F)
    abline(h=0 ,col="black")
    
  }else if (grepl("boxplot",plot_type)){
    ####empty plot
    plot(plotx,means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)-barwidth,max(plotx)+barwidth),xlab="",ylab="",xaxt="n",bty="n",xaxs="i",yaxs="i",type="n",lwd=line_min,cex.axis=min_cex,cex.lab=inter_cex,xaxt="n",yaxt="n")  
    par(bty="n")
    for (lvly in 1:nrow(means)){
      boxplot(diff[which(diff$time==means[lvly,"time"]&diff$treatment==means[lvly,"treatment"]&diff$predictor==means[lvly,"predictor"]),"variable"],at=plotx[lvly],add=T,range=1.5,notch=T,names=F,col=means[lvly,"final_col"],xlab="",ylab="",xaxt="n",xaxs="i",yaxs="i",cex.axis=min_cex,cex.lab=inter_cex,xaxt="n",yaxt="n",medlwd=line_max,boxlwd=line_min+0.5*(line_max-line_min),whisklwd=line_min,whisklty=1,staplelwd=line_min,boxwex=0.7,pch=16)
    }
    vals <- axisTicks(c(ymin,ymax),log=F,nint=6);interval <- mean(diff(vals))
    while(min(vals)>ymin){
      vals <- c(min(vals)-interval,vals)
    }
    while(max(vals)<ymax){
      vals <- c(vals,max(vals)+interval)
    }
    labs <- formatC(vals)
    #print(labs)
    if ("0" %in%labs ){
      idx0 <- which(labs=="0")
      if (is.even(idx0)){
        labs[!is.even(1:length(labs))] <-""
      }else{
        labs[is.even(1:length(labs))] <-""
      }
    }
    #print(labs)
    
    axis(2,at=vals,labels=labs,cex.axis=min_cex)
    
    par(xpd=F)
    abline(h=0 ,col="black")
    
  }else if (grepl("violinplot",plot_type)){
    ####empty plot
    plot(plotx,means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)-barwidth,max(plotx)+barwidth),xlab="",ylab="",xaxt="n",bty="n",xaxs="i",yaxs="i",type="n",lwd=line_min,cex.axis=min_cex,cex.lab=inter_cex,xaxt="n",yaxt="n")  
    par(bty="n")
    for (lvly in 1:nrow(means)){
      
      if (is.na(range)){
        VioPlot(na.omit(diff[which(diff$time==means[lvly,"time"]&diff$treatment==means[lvly,"treatment"]&diff$predictor==means[lvly,"predictor"]),"variable"]),col=alpha(means[lvly,"final_col"],0.7), horizontal=F, at=plotx[lvly], add=TRUE,lty=1, rectCol="black",wex=wex,border=NA,lwd=line_min,mode="Mean")
      }else{
        VioPlot(na.omit(diff[which(diff$time==means[lvly,"time"]&diff$treatment==means[lvly,"treatment"]&diff$predictor==means[lvly,"predictor"]),"variable"]),range=range, h=h,col=alpha(means[lvly,"final_col"],0.7), horizontal=F, at=plotx[lvly], add=TRUE,lty=1, rectCol="black",wex=wex,border=NA,lwd=line_min,mode="Mean")
      }
    }
    vals <- axisTicks(c(ymin,ymax),log=F,nint=6);interval <- mean(diff(vals))
    while(min(vals)>ymin){
      vals <- c(min(vals)-interval,vals)
    }
    while(max(vals)<ymax){
      vals <- c(vals,max(vals)+interval)
    }
    labs <- formatC(vals)
    #print(labs)
    if ("0" %in%labs ){
      idx0 <- which(labs=="0")
      if (is.even(idx0)){
        labs[!is.even(1:length(labs))] <-""
      }else{
        labs[is.even(1:length(labs))] <-""
      }
    }
    #print(labs)
    
    axis(2,at=vals,labels=labs,cex.axis=min_cex)
    
    par(xpd=F)
    abline(h=0 ,col="black")
    
  }else{
    ####empty plot
    plot(plotx,means$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)- (barwidth_fac_between*barwidth/2+barwidth/2),max(plotx)+ (barwidth_fac_between*barwidth/2+barwidth/2)),xlab="",ylab="",bty="n",xaxs="i",yaxs="i",type="n",lwd=line_min,cex.axis=min_cex,cex.lab=inter_cex,xaxt="n",yaxt="n")  
    ####arrows
    plot_arrows(means=means,plotx=plotx,plot_type=plot_type,LWD=line_max,LENGTH=0.025,colz=means$final_col)
    ####points
    points(plotx,means$mean,col=means$final_col,cex=1.5*max_cex,pch=16,lwd=line_min)  
    par(xpd=F)
    abline(h=0,lwd=line_max,col="black")
    vals <- axisTicks(c(ymin,ymax),log=F,nint=6);interval <- mean(diff(vals))
    while(min(vals)>ymin_ori){
      vals <- c(min(vals)-interval,vals)
    }
    while(max(vals)<ymax_ori){
      vals <- c(vals,max(vals)+interval)
    }
    labs <- formatC(vals)
    if ("0" %in%labs ){
      idx0 <- which(labs=="0")
      if (is.even(idx0)){
        labs[!is.even(1:length(labs))] <-""
      }else{
        labs[is.even(1:length(labs))] <-""
      }
    }
    
    axis(2,at=vals,labels=labs,lwd=0,lwd.ticks=line_inter,cex.axis=min_cex)
    abline(v=min(plotx)- (barwidth_fac_between*barwidth/2+barwidth/2),lwd=line_max,col="black")
  }
  if ((!(aligned&grepl(" low load",lab_title)))&((grepl(" load",lab_title)&!Extended)|!all(!grepl("time outside",lab_title))|!all(!grepl("with brood",lab_title))|!all(!grepl("Degree",lab_title))|!all(!grepl("time active",lab_title))|!all(!grepl("in the nest",lab_title)))){
    ytil <- list(quote("Pathogen-induced changes"),quote("relative to sham-induced changes"))
    par(xpd=NA)
    mtext(side=2,do.call(expression,ytil),line=pars$mgp[1]*(rev(1:length(ytil))-1)+pars$mgp[1],cex=par("cex") *inter_cex,xpd=NA)
    par(xpd=F)
  }
  if (grepl("\n",full_statuses_names[as.character(means$treatment)][1])){
    lassy <- 1
    mgpy <- par_mgp_ori+c(0,par_mgp_ori[1]-par_mgp_ori[2],0)
    fonty <- 3
    cexy <- min_cex
  }else{
    lassy <- 2
    mgpy <- c(0.8,0.05,0)
    fonty <- 1
    cexy <- inter_cex
  }
  
  par(xpd=T)
  axis(1,at=plotx,labels=full_statuses_names[as.character(means$treatment)],tick=F,cex.axis=cexy,las=lassy,mgp=mgpy,font=fonty)
  par(xpd=F)
  
  # title(ylab=lab_title,cex.lab=1.5,mgp=c(5,1,0))
  pvalue <- output[["interaction_problist"]][["24"]][["1"]]
  #title(main=from_p_to_ptext(pvalue),cex.main=max_cex,font.main=2,line=stat_line,xpd=T)
  # for (i in 1:(-1+length(unique(means$predictor)))){
  #   abline(v= mean(plotx[c(max(which(means$predictor==unique(means$predictor)[i])), min(which(means$predictor==unique(means$predictor)[i+1])))]),lwd=line_max)
  # }
  par(xpd=T)
  for (predy in unique(means$predictor)){
    textx <- mean(plotx[which(means$predictor==predy)])
    mtext(full_statuses_names[predy],side=1,line=0.5,adj=0.5,at=textx,cex=par("cex") *min_cex,col="black",font=1)
  }
  add_stats(dataset,plotx,means,ymin,ymax,predictor,p_colony,output,contrast.matrix,survival,p_queen,lab_title);
  par(xpd=NA)
  title(main=lab_title,cex.main=inter_cex,font.main=2,line=1+adjust_title_line,xpd=NA)
  par(xpd=F)
  
}

plot_regression <- function(data,time_point,analysis,n_cat_horiz,n_cat_vertic,pool=F,prepare=F,status=NULL,collective=NULL,pool_plot=F,input_color=NULL,plot_untransformed=F,boldy=F,aligned=F,ymin=NULL,ymax=NULL,xmin=NULL,xmax=NULL,point_cex=NULL,adjust_title_line=0,predict=NULL){
  adjust_title_line_ori <- adjust_title_line
  data_ori <- data
  #### plot regression for each desired combination of variable and predictor
  for (i in 1:length(analysis[["variable_list"]])){
    data <- data_ori
    adjust_title_line <- adjust_title_line_ori
    ####get variable and predictor
    variable <- analysis[["variable_list"]][i]
    
    ####if necessary: convert pixels to mm
    if (grepl("changetomm",names(variable))){
      if (grepl("changetomm2",names(variable))){
        data[,variable] <- data[,variable]*pix_to_mm_ratio*pix_to_mm_ratio
        names(variable) <- gsub("_changetomm2","",names(variable))
      }else{
        data[,variable] <- data[,variable]*pix_to_mm_ratio
        names(variable) <- gsub("_changetomm","",names(variable))
      }
    }
    
    print(variable)
    
    predictor <- analysis[["predictor_list"]][i]
    transf_variable <- analysis[["transf_variable_list"]][i]
    transf_predictor <- analysis[["transf_predictor_list"]][i]
    pooli <- pool[i]
    
    ####if necessary: include queen and treated into predictor function
    if (!is.null(predictor)){
      if ((!collective&refine!=""&predictor!="")&("tag"%in%names(data))){
        ####first change the content of predictor column 
        data["predictor"] <- data[,predictor]
        ###second add treated
        data[which(data$status=="treated"),"predictor"] <- "treated"
        ####fourth copy back into predictor column
        data[,predictor] <- data[,"predictor"]
        ####fifth if necessary remove queen
        if (length(unique(data$predictor))>1){
          if (!queen){
            data <- data[which(data$tag!=queenid),]
          }
          if (!treated){
            data <- data[which(data$predictor!="treated"),]
          }
          if (!nurses){
            data <- data[which(data$predictor!="nurse"),]
          }
          if (!foragers){
            data <- data[which(data$predictor!="forager"),]
          }
        }
      }else if (predictor!=""){
        data["predictor"] <- data[,predictor]
      }  
    }
    
    ####if necessary: apply prepare dataset function
    if (prepare){
      data <- prepare_dataset(data,variable)
    }else{
      ####process variable
      data["variable"] <- data[,variable]
      data[which(!is.finite(data$variable)),"variable"] <- NA
    }
    data["untransformed_variable"] <- data$variable
    ####transform variable
    ylabel <- names(variable)
    ylabel <- capitalize(ylabel)
    if (transf_variable=="log"){
      print("Logging variable...")
      data[!is.na(data$variable),"variable"] <- log_transf(data[!is.na(data$variable),"variable"] )
      if (!plot_untransformed){
        
        if(boldy){
          ylabel <- substitute(paste(bold(ylabel),bolditalic(" (log)")),list(ylabel=ylabel))
          adjust_title_line <- 0.17
        }else{
          ylabel <- substitute(paste(ylabel,italic(" (log)")),list(ylabel=ylabel))
          adjust_title_line <- 0.17
        }
        
      }
    }else if (grepl("power",transf_variable)){
      data[!is.na(data$variable),"variable"]  <- (data[!is.na(data$variable),"variable"] )^as.numeric(gsub("power","",transf_variable))
      if (!plot_untransformed){
        
        if(boldy){
          ylabel <- substitute(paste(bold(ylabel),bolditalic(" (") ^pow,bolditalic(")")),list(ylabel=ylabel,pow=as.numeric(gsub("power","",transf_variable))))
          adjust_title_line <- 0
        }else{
          ylabel <- substitute(paste(ylabel,italic(" (") ^pow,italic(")")),list(ylabel=ylabel,pow=as.numeric(gsub("power","",transf_variable))))
          adjust_title_line <- 0
        }
        
      }
    }else if (transf_variable=="sqrt"){
      data[!is.na(data$variable),"variable"]  <- sqrt_transf(data[!is.na(data$variable),"variable"] )
      if (!plot_untransformed){
        
        if (boldy){
          ylabel <- substitute(paste(bold(ylabel),bolditalic(" ("),sqrt(bolditalic(")"))),list(ylabel=ylabel))
          adjust_title_line <- 0.24
        }else{
          ylabel <- substitute(paste(ylabel,italic(" ("),sqrt(italic(")"))),list(ylabel=ylabel))
          adjust_title_line <- 0.24
        }
      }
      
    }
    ####for comparison: use perform_analysis_combined_function
    if (time_point=="comparison"){
      if ("randy"%in%names(data)){
        case <- "case3"
      }else{
        
        if (is.null(predictor)){
          case <- "case1"
        }else if (length(unique(data$predictor))==1){
          case <- "case1"
        }else if (!((!collective&refine!=""&predictor!=""))){
          case <- "case1"
        }else{
          case <- "case2"
        }
        
      }
      if (case=="case1"){
        perform_barplot_analysis(root_path=root_path,collective=collective,dataset=data,lab_title=ylabel,excluded=NULL,pool=pool,violin_params=analysis[["violin_plot_param"]][i],pool_plot=pool_plot,adjust_title_line=adjust_title_line)
      }else if (case=="case2"){
        perform_barplot_analysis_refined(root_path=root_path,collective=collective,dataset=data,lab_title=ylabel,excluded=NULL,pool=pool,violin_params=analysis[["violin_plot_param"]][i],pool_plot=pool_plot,plot_untransformed=plot_untransformed,aligned=aligned,adjust_title_line=adjust_title_line)
      }else{
        perform_barplot_analysis_simple(root_path=root_path,collective=collective,dataset=data,lab_title=ylabel,excluded=NULL,pool=pool,violin_params=analysis[["violin_plot_param"]][i],pool_plot=pool_plot,plot_untransformed=plot_untransformed,aligned=aligned,adjust_title_line=adjust_title_line)
      }
    }else{
      ####process predictor
      data["predictor"] <- data[,predictor]
      if (is.numeric(data$predictor)){
        data[which(!is.finite(data$predictor)),"predictor"] <- NA
        if (transf_predictor=="log"){
          if (!is.null(predict)){
            transfi <- log_transf(c(data[!is.na(data$predictor),"predictor"] , predict))
            predict <- transfi[length(transfi)]
            data[!is.na(data$predictor),"predictor"]  <- transfi[1:(length(transfi)-1)]
          }else{
            data[!is.na(data$predictor),"predictor"]  <- log_transf(data[!is.na(data$predictor),"predictor"] )
          }
          xlabel<- paste("Log(", names(predictor),")",sep="")
        }else if (transf_predictor=="power2"){
          data[!is.na(data$predictor),"predictor"]  <- (data[!is.na(data$predictor),"predictor"] )^2
          if (!is.null(predict)){predict <- predict^2}
          xlabel <- substitute( xlabely ^2,list(xlabely=names(predictor)))
        }else if (transf_predictor=="sqrt"){
          if (!is.null(predict)){
            transfi <- sqrt_transf(c(data[!is.na(data$predictor),"predictor"] , predict))
            predict <- transfi[length(transfi)]
            data[!is.na(data$predictor),"predictor"]  <- transfi[1:(length(transfi)-1)]
          }else{
            data[!is.na(data$predictor),"predictor"]  <- sqrt_transf(data[!is.na(data$predictor),"predictor"] )
          }
          xlabel <- substitute(sqrt ( xlabely),list(xlabely=names(predictor)))
        }else{
          xlabel <- names(predictor)
        }
      }else{
        xlabel <- ""
      }
      
      title <- ""
      if(!"period"%in%names(data)){data["period"] <- data$time}
      
      ###process pool argument
      if (pooli){
        dat <- aggregate(na.rm=T,na.action="na.pass",variable~predictor+colony+antid+status+colony_size+period,FUN="mean",data=data)
      }else{
        dat <- data
      }
      
      ###process horizontal bin argument
      if (length(unique(dat[,"predictor"]))>n_cat_horiz){
        dat[which(!is.na(dat[,"predictor"])),"predictor_plot"] <- as.numeric(gsub("\\(","",unlist(strsplit(as.character(cut(dat$predictor,n_cat_horiz,include_lowest=T)),split=","))[grepl("\\(",unlist(strsplit(as.character(cut(dat$predictor,n_cat_horiz,include_lowest=T)),split=",")))]))
      }else{#if (length(unique(dat[,"predictor"]))>n_cat_horiz)
        dat[,"predictor_plot"] <- dat[,"predictor"]
      }##else
      
      ###define formula
      formula_stat <- as.formula(paste("variable"," ~ ", paste(c("predictor","colony_size","(1|colony)","(1|antid)"), collapse= "+")))
      if (length(unique(na.omit(aggregate(variable~antid,FUN=length,data=dat)$variable)))==1&unique(na.omit(aggregate(variable~antid,FUN=length,data=dat)$variable))[1]==1){
        formula_stat <- update(formula_stat,~.-(1|antid))
      }
      formula_plot <- as.formula(paste("variable"," ~ ", paste(c("predictor_plot","period"), collapse= "+")))
      ###plot
      predicted_value <- scatterplot_violin_forpaper(formula_stat=formula_stat,formula_plot=formula_plot,ylabel=ylabel,xlabel=xlabel,title=title,dat=dat,sorting="period",time_point=time_point,output=F,violin_params = analysis[["violin_plot_param"]][i],input_color=input_color,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,point_cex=point_cex,predict=predict)
      if (!is.null(predicted_value)){
        predicted_value <- data[,variable][closest_match(predicted_value,data$variable) ]
        return(predicted_value)
      } 
    }
  }
}

plot_seeds <- function(experiments,seeds,variables,transf,color_pal){
  collective_data <- NULL
  for (seed in seeds){
    for (experiment in experiments){
      setwd(paste(disk_path,"/",experiment,"/transmission_simulations/random_vs_observed/",seed,sep=""))
      temp <- data.frame(experiment=experiment,seed=seed,read.table("collective_simulation_results_observed.txt",header=T,stringsAsFactors = F),stringsAsFactors = F)
      if (!is.null(collective_data)){
        collective_data <-collective_data[,names(collective_data)%in%names(temp)]
        temp <- temp[names(collective_data)]
      }
      collective_data <- rbind(collective_data,temp)
    }
  }
  collective_data["colony"] <- as.character(interaction(collective_data$experiment,collective_data$colony))
  
  for (variable in variables){
    collective_data["variable"] <- collective_data[,variable]
    dat <- collective_data[c("colony","seed","variable")]
    transf_variable <- transf[which(variable==variables)]
    ylabel <- names(variables[variables==variable])
    
    if (transf_variable=="log"){
      print("Logging variable...")
      dat[!is.na(dat$variable),"variable"] <- log_transf(dat[!is.na(dat$variable),"variable"] )
      ylabel <- substitute(paste(ylabel,italic(" (log)")),list(ylabel=ylabel))
      
    }else if (grepl("power",transf_variable)){
      dat[!is.na(dat$variable),"variable"]  <- (dat[!is.na(dat$variable),"variable"] )^as.numeric(gsub("power","",transf_variable))
      ylabel <- substitute(paste(ylabel,italic(" (") ^pow,italic(")")),list(ylabel=ylabel,pow=as.numeric(gsub("power","",transf_variable))))
    }else if (transf_variable=="sqrt"){
      dat[!is.na(dat$variable),"variable"]  <- sqrt_transf(dat[!is.na(dat$variable),"variable"] )
      ylabel <- substitute(paste(ylabel,italic(" ("),sqrt(italic(")"))),list(ylabel=ylabel))
    }
    
    
    
    
    
    dat$seed <- factor(dat$seed,levels=seeds)
    mod <- lmer(variable~seed+(1|colony),data=dat)
    
    pval <- Anova(mod)["seed","Pr(>Chisq)"]
    ###now plot
    forplot <- data.frame(as.matrix(aggregate(variable~seed,function(x)cbind(mean(x),std.error(x)),data=dat)),stringsAsFactors = F)
    names(forplot) <- c("seed","mean","se")
    forplot$mean <- as.numeric(forplot$mean);forplot$se <- as.numeric(forplot$se)
    
    ymin <- min(c(forplot$mean-forplot$se),na.rm=T);ymax<- max(c(forplot$mean+forplot$se),na.rm=T)
    if (ymin>0){ymin <- 0};if (ymax<0){ymax <- 0}
    if (grepl("point",plot_type)|grepl("boxplot",plot_type)|grepl("violinplot",plot_type)){
      rangy <- max(dat$variable,na.rm=T)-min(dat$variable,na.rm=T)
      ymin <- min(ymin, min(dat$variable,na.rm=T)-0.1*rangy)
      ymax <- max(ymax, max(dat$variable,na.rm=T)+0.1*rangy)
      
    }
    rangy <- ymax-ymin
    ymin <- ymin-0.005*rangy;ymax <- ymax+0.20*rangy
    barwidth <- 0.5
    barspace <- 0.5
    arrowcodes <- c(1,2,3); names(arrowcodes) <- c("-1","1","0")
    plotx <- barplot(forplot$mean,plot=F,width=barwidth,space=barspace)  
    if (grepl("bars",plot_type)){
      plotx <- barplot(forplot$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)-barwidth,max(plotx)+barwidth),xlab="",ylab=ylabel,xaxt="n",bty="n",xaxs="i",yaxs="i",cex=0.5,cex.axis=min_cex,cex.lab=inter_cex,width=barwidth,space=barspace,col=color_pal[forplot$seed],xaxt="n")  
      plot_arrows(means=forplot,plotx=plotx,plot_type=plot_type,LWD=line_max,LENGTH=0.05,colz=color_pal[forplot$seed])
    }else if (grepl("boxplot",plot_type)){
      ####empty plot
      plot(plotx,forplot$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)-barwidth,max(plotx)+barwidth),xlab="",ylab=ylabel,xaxt="n",bty="n",xaxs="i",yaxs="i",type="n",lwd=line_min,cex.axis=min_cex,cex.lab=inter_cex,xaxt="n")  
      par(bty="n")
      for (lvly in 1:nrow(forplot)){
        boxplot(dat[which(dat$seed==forplot[lvly,"seed"]),"variable"],at=plotx[lvly],add=T,range=1.5,notch=T,names=F,col=color_pal[forplot[lvly,"seed"]],xlab="",ylab="",xaxt="n",xaxs="i",yaxs="i",cex.axis=min_cex,cex.lab=inter_cex,xaxt="n",yaxt="n",medlwd=line_max,boxlwd=line_min+0.5*(line_max-line_min),whisklwd=line_min,whisklty=1,staplelwd=line_min,boxwex=0.7,pch=16)
      }
      
    }else{
      ####empty plot
      plot(plotx,forplot$mean,ylim=c(ymin,ymax),xlim=c(min(plotx)-barwidth,max(plotx)+barwidth),xlab="",ylab=ylabel,xaxt="n",bty="n",xaxs="i",yaxs="i",type="n",lwd=line_min,cex.axis=min_cex,cex.lab=inter_cex,xaxt="n")  
      ####arrows
      plot_arrows(means=forplot,plotx=plotx,plot_type=plot_type,LWD=line_max,LENGTH=0.05,colz=color_pal[forplot$seed])
      ####points
      points(plotx,forplot$mean,cex=1.5*max_cex,pch=16,lwd=line_min,col=color_pal[forplot$seed])  
    }
    
    xlabel <- names(seeds)
    ###if necessary cut xlab and ylab in 2 lines
    for (k in 1:length(xlabel)){
      temp <- xlabel[k]
      if (nchar(temp)>8){
        ###split label according to spaces
        temp <- unlist(strsplit(temp,split=" "))
        ###count nb of characters in each word
        temp_length <- unlist(lapply(temp,FUN=nchar))
        ###find roughly the middle
        cutsy <- sum(temp_length/2)
        cumul_length <- 0;cut_index <- NULL
        for (idx in c(1:length(temp))){
          cumul_length <- cumul_length + temp_length[idx]
          if(cumul_length>=cutsy & is.null(cut_index)){
            cut_index <- idx-1
          }#if(cumul_length>=cutsy & is.null(cut_index))
        }#for (idx in c(1:length(temp)))
        temp <- paste(paste(temp[c(1:cut_index)],collapse=" "),paste(temp[c((cut_index+1):length(temp))],collapse=" "),sep=" \n ")
        xlabel[k] <- temp
      }
    }#for (spec_lab in c("xlabel","ylabel"))
    
    par_mgp <- par()$mgp
    par(mgp=c(1.3,0.4,0))
    axis(1,at=plotx,labels=xlabel,tick=F,cex.axis=min_cex)
    title(xlab="Simulation seeds",cex.lab=inter_cex)
    par(mgp=par_mgp)
    par(xpd=F)
    abline(h=0)
    
    if (pval>0.05){p_cex <- inter_cex;adjust_line <- 0.4;fonty <- 1}else{p_cex <- max_cex*1.1;adjust_line <- 0.1; fonty <-  2}
    
    title(main=from_p_to_ptext(pval),cex.main=p_cex,font.main=fonty,line=stat_line+adjust_line,xpd=T)
    if (pval<0.05){
      post_hoc <- summary(glht(mod,linfct = mcp(seed="Tukey")),test=adjusted("BH"))
      post_hoc_groups <- cld(post_hoc)$mcletters$Letters
      post_hoc_groups <- post_hoc_groups[forplot$seed]
      for (i in 1:length(post_hoc_groups)){
        mtext(post_hoc_groups[i],side=3,line=stat_line-0.7,at=plotx[i],xpd=T,cex=par("cex") *inter_cex,font=2)
      }
    }
  }
}

points_plus_boxplot <- function(forplots,var,pred,levels,colour_pal,ylaby,pvalue,labels,option="boxplots",pool_plot){
  par_mar_ori <- par()$mar
  par(mar=par_mar_ori+c(0,0,0,0.5))
  
  if (pool_plot){
    jitter_param <- 0.15
  }else{
    jitter_param <- 0.3
  }
  forplot <- NULL; means <- NULL; ses <- NULL
  for (lvl in 1:length(forplots)){
    forplots[[lvl]]["var"] <- forplots[[lvl]][,var]
    forplots[[lvl]]["pred"] <- forplots[[lvl]][,pred]
    temp <- forplots[[lvl]]
    if (option=="barplots"){
      means <- rbind(means,aggregate(na.rm=T,na.action="na.pass",var~pred,FUN="mean",data=temp));ses <- rbind(ses,aggregate(na.rm=T,na.action="na.pass",var~pred,FUN="std.error",data=temp));
    }
    if (pool_plot){
      temp <- aggregate(var~colony+time+pred,FUN=mean,data=temp)
    }
    forplot <- rbind(forplot,temp)
  }
  
  forplot$pred <- factor(forplot$pred,levels = levels)  
  
  ###get colony ordering from observed
  if (length(unique(forplot$colony))==length(unique(forplots[[length(forplots)]]$colony))){
    col_means <- aggregate(var~colony,FUN=mean,data=forplots[[length(forplots)]])
  }else{
    col_means <- aggregate(var~colony,FUN=mean,data=forplot)
  }
  col_means <- col_means [order(col_means$var),]
  coli_list <- col_means$colony
  # colour_pal <- viridis(length(coli_list),alpha=0.5,begin=0.4,end=0.8,option="inferno")
  colour_pal <- colour_pal(length(coli_list))
  par(pars)  
  par(bty="n",xaxt = "n")
  
  for (coli in rev(coli_list)){
    if (coli==coli_list[length(coli_list)]){
      addy <- F
    }else{
      addy <- T
    }
    
    stripchart(var ~ pred, data = forplot[forplot$colony==coli,],vertical = TRUE,pch=16,col=alpha(colour_pal[which(coli==coli_list)],1),method = 'jitter', jitter = jitter_param,ylim=c(min(c(forplot$var)),max(c(forplot$var))), main = "",ylab=ylaby,add=addy,bty="l",cex=0.6,cex.axis=min_cex,cex.lab=inter_cex)
  }
  if (option=="boxplots"){
    ###make boxplot
    boxplot(var ~ pred, data = forplot, 
            outline = FALSE,     ## avoid double-plotting outliers, if any
            main = "",yaxt="n",add=T,col=alpha("white",0),medlwd=line_max,boxlwd=line_min+0.5*(line_max-line_min),whisklwd=line_min,whisklty=1,staplelwd=line_min,boxwex=0.7,bty="l")
  }else if (option=="barplots"){
    names(means)[names(means)=="var"] <- "mean";names(ses)[names(ses)=="var"] <- "se";means <- merge(means,ses)
    means <- means[match(means$pred,levels),]
    for (lvl in 1:length(levels)){
      rect(lvl - 0.008, means[lvl,"mean"]-means[lvl,"se"], lvl + 0.008, 
           means[lvl,"mean"]+means[lvl,"se"], col = "black")
      points(lvl, means[lvl,"mean"], pch = 21, col = "black",bg="white",cex=0.8)
    }
  }
  if (!is.na(pvalue)){
    par(xpd=T)
    title(main=from_p_to_ptext(pvalue),cex.main=max_cex,font.main=2,line=-0.25,xpd=T)
    par(xpd=F)
  }
  par(xaxt = "s")
  axis(side=1,at=c(1,2),labels=labels,tick=F,lty=0,cex.axis=inter_cex)
  par(mar=par_mar_ori)
}

prepare_dataset <- function(overall_results,variable,survival=F){
  if(is.null(overall_results)){
    dataset <- NULL
  }else{
    dataset <- overall_results
    dataset["variable"] <- dataset[variable]
    dataset <- dataset[(!is.na(dataset$variable))&(dataset$variable!=Inf),]
    names(dataset)[names(dataset)=="period"] <- "time"
    dataset$time <- as.character(dataset$time)
    dataset$time <- capitalize(dataset$time)
    dataset$time_of_day <- factor(dataset$time_of_day,levels=as.character(sort(unique(dataset$time_of_day))))
    dataset$treatment <- factor(dataset$treatment,levels=as.character(sort(unique(dataset$treatment))))
    if(survival){
      dataset["censor_variable"] <- dataset[paste("cens_",variable,sep="")]
      #####if survival, one must ensure that things are comparable between the 2 time windows
      #####so for each time of day get the smallest period considered
      min_time_considered <- aggregate(na.rm=T,na.action="na.pass",duration_considered_h~colony+time_of_day+tag,FUN=min,data=dataset)
      #####and then if recorded period is larger, replace it by the smaller value and censor it
      for (idx in 1:nrow(dataset)){
        if (dataset[idx,"variable"]>min_time_considered[min_time_considered$colony==dataset[idx,"colony"]&min_time_considered$time_of_day==dataset[idx,"time_of_day"]&min_time_considered$tag==dataset[idx,"tag"],"duration_considered_h"]){
          dataset[idx,"variable"] <- min_time_considered[min_time_considered$colony==dataset[idx,"colony"]&min_time_considered$time_of_day==dataset[idx,"time_of_day"]&min_time_considered$tag==dataset[idx,"tag"],"duration_considered_h"]
          dataset[idx,"censor_variable"] <- 0
        }
      }
    }
    if (exists("refine")){
      if ("tag"%in%names(dataset)){
        if (refine=="caste_health_care.txt"){
          dataset <- dataset[which(dataset$tag!=665),]
        }
      }
    }
  }
  if (("tag"%in%names(dataset))&(!"antid"%in%names(dataset))){
    dataset["antid"] <- as.character(interaction(dataset$colony,dataset$tag))
  }
  return(dataset)
}

prepare_stats_1 <- function(collective,dataset,type,predictor,inter,survival){
  if (survival){statista <- "Pr(>|Chi|)"}else{statista <- "Pr(>Chisq)" } 
  if(survival){first_term <- "Surv(variable,censor_variable) ~ "}else{first_term <- "variable ~ "}
  if(!collective){
    if (type=="individual"){
      form_stat <- as.formula(paste(first_term, paste(c("time*predictor","colony_size","(1|time_of_day)","(1|colony)","(1|antid)"), collapse= "+")))
      if (nrow(unique(dataset[c("time_of_day","time")]))==nrow(unique(dataset[c("time")]))){
        form_stat <- as.formula(paste(first_term, paste(c("time*predictor","colony_size","(1|colony)","(1|antid)"), collapse= "+")))
      }
      if (length(unique(dataset$antid))==1){
        form_stat <- update(form_stat,.~.-(1|antid))        
      }
    }else{#if (type=="individual")
      form_stat <- as.formula(paste(first_term, paste(c("time*predictor","colony_size","(1|time_of_day)","(1|colony)"), collapse= "+")))
      if (nrow(unique(dataset[c("time_of_day","time")]))==nrow(unique(dataset[c("time")]))){
        form_stat <- as.formula(paste(first_term, paste(c("time*predictor","colony_size","(1|colony)"), collapse= "+")))
      }
    }#else#if (type=="individual")
  }else{#!collective
    form_stat <- as.formula(paste(first_term, paste(c(inter,"colony_size","(1|time_of_day)","(1|colony)"), collapse= "+")))
    if (nrow(unique(dataset[c("time_of_day","time")]))==nrow(unique(dataset[c("time")]))){
      form_stat <- as.formula(paste(first_term, paste(c(inter,"colony_size","(1|colony)"), collapse= "+")))
    }
  }#!collective
  ######statistics
  if(exists("model")){rm(list=c("model"))} 
  if(exists("p_colony")){rm(list=c("p_colony"))}
  if(!survival){
    try(model <- do.call(lmer, list(formula=form_stat, data=dataset)),silent=T)
    try(anov <- Anova(model),silent=T)
    try(p_colony <- anov[statista]["colony_size",statista],silent=T)
  }else{
    try(model <- do.call(coxme, list(formula=form_stat, data=dataset)),silent=T)
    try(p_colony <- printme_coxme(model)["colony_size"],silent=T)
  }
  ####test effect of colony size ###
  ######1/case when the model was not identifiable: remove colony size
  if((!exists("model"))|(!exists("p_colony"))){
    form_stat <- update(form_stat,.~.-colony_size)
    rm(list=c("model"))
    p_colony <- NA
  }else{####2/case where model was identifiable
    if (is.na(p_colony)|((!is.na(p_colony))&(p_colony>0.05))){
      form_stat <- update(form_stat,.~.-colony_size)
      rm(list=c("model"))
    }
  }
  
  if (length(unique(dataset$predictor))==3){
    if ("colony_size"%in%all.vars(form_stat)){
      if (grepl("age",root_path)){
        contrast.matrix <- rbind(  "change in untreated minus change in queen" =    c(0, 0, 0, 0, 0, 0,-1),
                                   "change in treated minus change in queen" =      c(0, 0, 0, 0, 0,-1, 0),
                                   "change in untreated minus change in treated" =  c(0, 0, 0, 0, 0, 1,-1),
                                   "queen after minus queen before" =               c(0,-1, 0, 0, 0, 0, 0),
                                   "treated after minus treated before" =           c(0,-1, 0, 0, 0,-1, 0),
                                   "untreated after minus untreated before" =       c(0,-1, 0, 0, 0, 0,-1)                                      
        )
      }else{
        contrast.matrix <- rbind(  "change in pathogen minus change in control" =  c(0, 0, 0, 0, 0, 0,-1),
                                   "change in killeds minus change in control" =  c(0, 0, 0, 0, 0,-1, 0),
                                   "change in pathogen minus change in killeds" =  c(0, 0, 0, 0, 0, 1,-1),
                                   "pathogen after minus pathogen before" =         c(0,-1, 0, 0, 0, 0,-1),
                                   "killeds after minus killeds before" =         c(0,-1, 0, 0, 0,-1, 0),                                     
                                   "control after minus control before" =         c(0,-1, 0, 0, 0, 0, 0)                                       
        )
        
      }
    }else{
      if (grepl("age",root_path)){
        contrast.matrix <- rbind(  "change in untreated minus change in queen" =    c(0, 0, 0, 0, 0,-1),
                                   "change in treated minus change in queen" =      c(0, 0, 0, 0,-1, 0),
                                   "change in untreated minus change in treated" =  c(0, 0, 0, 0, 1,-1),
                                   "queen after minus queen before" =               c(0,-1, 0, 0, 0, 0),
                                   "treated after minus treated before" =           c(0,-1, 0, 0,-1, 0),
                                   "untreated after minus untreated before" =       c(0,-1, 0, 0, 0,-1)                                      
        )
      }else{
        contrast.matrix <- rbind(  "change in pathogen minus change in control" =    c(0, 0, 0, 0, 0,-1),
                                   "change in killeds minus change in control" =    c(0, 0, 0, 0,-1, 0),
                                   "change in pathogen minus change in killeds" =    c(0, 0, 0, 0, 1,-1),
                                   "pathogen after minus pathogen before" =           c(0,-1, 0, 0, 0,-1),
                                   "killeds after minus killeds before" =           c(0,-1, 0, 0,-1, 0),                                     
                                   "control after minus control before" =           c(0,-1, 0, 0, 0, 0)                                       
        )
      }
    }
    
  }else{
    if ("colony_size"%in%all.vars(form_stat)){
      if (!grepl("age",root_path)){        contrast.matrix <- rbind(  "change in pathogen minus change in control" =  c(0, 0, 0, 0,-1),
                                                                      "pathogen after minus pathogen before" =         c(0,-1, 0, 0,-1),
                                                                      "control after minus control before" =         c(0,-1, 0, 0, 0)                                       
      )
      }else{
        contrast.matrix <- rbind(  "change in treated minus change in untreated" =  c(0, 0, 0, 0,1),
                                   "treated after minus treated before" =         c(0,-1, 0, 0,0),
                                   "untreated after minus untreated before" =         c(0,-1, 0, 0, -1)                                       
        )
      }
    }else{
      if (!grepl("age",root_path)){
        contrast.matrix <- rbind(  "change in pathogen minus change in control" =  c(0, 0, 0,-1),
                                   "pathogen after minus pathogen before" =         c(0,-1, 0,-1),
                                   "control after minus control before" =         c(0,-1, 0, 0)                                       
        )
      }else{
        contrast.matrix <- rbind(  "change in treated minus change in untreated" =  c(0, 0, 0, 1),
                                   "treated after minus treated before" =         c(0,-1, 0, 0),
                                   "untreated after minus untreated before" =         c(0,-1, 0, -1)                                       
        )
      }
    }
  }
  if(survival){contrast.matrix <- contrast.matrix[,2:ncol(contrast.matrix)]}
  return(list(form_stat=form_stat,p_colony=p_colony,contrast.matrix=contrast.matrix))
}

prepare_stats_2 <- function(dataset,form_stat,survival,is_queens=F,contrast.matrix=NULL){
  if (survival){statista <- "Pr(>|Chi|)"}else{statista <- "Pr(>Chisq)" } 
  interaction_problist <- vector("list", length(unique(c(24))))
  names(interaction_problist) <- as.character(sort(unique(c(24))))
  predictor_problist <- vector("list", length(unique(c(24))))
  names(predictor_problist) <- as.character(sort(unique(c(24))))
  time_problist <- vector("list", length(unique(c(24))))
  names(time_problist) <- as.character(sort(unique(c(24))))
  modellist <- vector("list", length(unique(c(24))))
  names(modellist) <- as.character(sort(unique(c(24))))
  for (time_window_ref in sort(unique(c(24)))){
    if (time_window_ref==unit){
      formi <- update(form_stat,.~.-(1|time_of_day_bis))
    }else{
      formi <- form_stat
    }
    nb_time_windows <- floor(24/time_window_ref)
    interaction_probs <- NULL;time_probs <- NULL;predictor_probs <- NULL
    models <- vector("list", nb_time_windows);names(models) <- as.character(1:nb_time_windows)
    for (time_point in 1:nb_time_windows){
      ###subset dataset
      subset <- dataset[(as.numeric(as.character(dataset$time_of_day_bis))>=(time_window_ref*(time_point-1)))&(as.numeric(as.character(dataset$time_of_day_bis))<(time_window_ref*time_point)),]
      if (exists("model")){rm=list=c("model")}
      gc()
      if (!survival){
        try(model <- do.call(lmer, list(formula=formi, data=subset)),silent=T)
        try(anov <- Anova(model),silent=T)
      }else{
        try(model <- do.call(coxme, list(formula=formi, data=subset)),silent=T)
        try(anov <- anova(model),silent=T)
      }
      
      if(exists("model")){
        if (time_window_ref==24){test_norm(residuals(model));print(anov)}
        
        models[[as.character(time_point)]] <- model
        gc()
        if(exists("pees")){rm(list=c("pees"))}
        try(pees <- anov[statista]["time",statista],silent=T)
        if (exists("pees")){
          p_time <- anov[statista]["time",statista]
          p_predictor <- anov[statista]["predictor",statista]
          p_interaction <- anov[statista]["time:predictor",statista]
          rm(list=c("pees"))
        }else{
          p_time <- NA
          p_predictor <- NA
          p_interaction <- NA
        }
        
        if (is_queens){
          post_hoc <- summary(glht(model,contrast.matrix),test=adjusted("BH"))
          print("z value");print(post_hoc$test$tstat);print("Pr>|z|");print(post_hoc$test$pvalues);
        }
        rm(list=c("model"))
      }else{
        models[[as.character(time_point)]] <- NULL
        p_time <- NA
        p_predictor <- NA
        p_interaction <- NA
      }
      names(p_interaction) <- time_point;interaction_probs <- c(interaction_probs,p_interaction)
      names(p_predictor) <- time_point;predictor_probs <- c(predictor_probs,p_predictor)
      names(p_time) <- time_point;time_probs <- c(time_probs,p_time)
    }
    interaction_probs <- p.adjust(interaction_probs,method="hochberg"); time_probs <- p.adjust(time_probs,method="hochberg"); predictor_probs <- p.adjust(predictor_probs,method="hochberg");
    interaction_problist[[as.character(time_window_ref)]] <- interaction_probs;    time_problist[[as.character(time_window_ref)]] <- time_probs;predictor_problist[[as.character(time_window_ref)]] <- predictor_probs
    modellist[[as.character(time_window_ref)]] <- models
  }
  return(list(interaction_problist=interaction_problist,time_problist=time_problist,predictor_problist=predictor_problist,modellist=modellist))
}

prepare_stats_3 <- function(dataset,predictor,survival){
  if (survival){statista <- "Pr(>|Chi|)"}else{statista <- "Pr(>Chisq)" } 
  if(survival){first_term <- "Surv(variable,censor_variable) ~ "}else{first_term <- "variable ~ "}
  form_stat <- as.formula(paste(first_term, paste(c("time*treatment*predictor","colony_size","(1|time_of_day)","(1|colony)","(1|antid)"), collapse= "+")))
  ######statistics
  if(exists("model")){rm(list=c("model"))} 
  if(exists("p_colony")){rm(list=c("p_colony"))}
  if(!survival){
    try(model <- do.call(lmer, list(formula=form_stat, data=dataset)),silent=T)
    try(anov <- Anova(model),silent=T)
    try(p_colony <- anov[statista]["colony_size",statista],silent=T)
  }else{
    try(model <- do.call(coxme, list(formula=form_stat, data=dataset)),silent=T)
    try(p_colony <- printme_coxme(model)["colony_size"],silent=T)
  }
  ####test effect of colony size ###
  ######1/case when the model was not identifiable: remove colony size
  if((!exists("model"))|(!exists("p_colony"))){
    form_stat <- update(form_stat,.~.-colony_size)
    rm(list=c("model"))
    p_colony <- NA
  }else{####2/case where model was identifiable
    if (is.na(p_colony)|((!is.na(p_colony))&(p_colony>0.05))){
      form_stat <- update(form_stat,.~.-colony_size)
      rm(list=c("model"))
    }
  }
  
  #####define contrast matrix 
  ###first get the number and names of the levels
  level_names <- levels(dataset$predictor)
  treatment_names <- levels(dataset$treatment)
  time_names <- unique(dataset$time)
  ###deduce the index of colony size in the parameter matrix 
  colony_index <- (1 + length(time_names)-1 + length(treatment_names)-1 + length(level_names)-1+ 1)
  if (length(level_names)==2){
    contrast.matrix <-   rbind(
      "Delta_treatment1_1 - Delta_treatment1_2" =       c(0,0,0,0,0,0,1,0,0), 
      "Delta_treatment1_1 - Delta_treatment2_1" =       c(0,0,0,0,0,1,0,0,0),
      "Delta_treatment1_1 - Delta_treatment2_2" =       c(0,0,0,0,0,1,1,0,1),
      "Delta_treatment1_2 - Delta_treatment2_1" =       c(0,0,0,0,0,1,-1,0,0),
      "Delta_treatment1_2 - Delta_treatment2_2" =       c(0,0,0,0,0,1,0,0,1),
      "Delta_treatment2_1 - Delta_treatment2_2" =       c(0,0,0,0,0,0,1,0,1)
    )
  }else if (length(level_names)==3){
    contrast.matrix <-   rbind(
      "Delta_treatment1_1 - Delta_treatment1_2" =       c(0,0,0,0,0,0,0,1,0,0,0,0,0),
      "Delta_treatment1_1 - Delta_treatment1_3" =       c(0,0,0,0,0,0,0,0,1,0,0,0,0),
      "Delta_treatment1_1 - Delta_treatment2_1" =       c(0,0,0,0,0,0,1,0,0,0,0,0,0),
      "Delta_treatment1_1 - Delta_treatment2_2" =       c(0,0,0,0,0,0,1,1,0,0,0,1,0),
      "Delta_treatment1_1 - Delta_treatment2_3" =       c(0,0,0,0,0,0,1,0,1,0,0,0,1),
      "Delta_treatment1_2 - Delta_treatment1_3" =       c(0,0,0,0,0,0,0,-1,1,0,0,0,0),
      "Delta_treatment1_2 - Delta_treatment2_1" =       c(0,0,0,0,0,0,1,-1,0,0,0,0,0),
      "Delta_treatment1_2 - Delta_treatment2_2" =       c(0,0,0,0,0,0,1,0,0,0,0,1,0),
      "Delta_treatment1_2 - Delta_treatment2_3" =       c(0,0,0,0,0,0,1,-1,1,0,0,0,1),
      "Delta_treatment1_3 - Delta_treatment2_1" =       c(0,0,0,0,0,0,1,0,-1,0,0,0,0),
      "Delta_treatment1_3 - Delta_treatment2_2" =       c(0,0,0,0,0,0,1,1,-1,0,0,1,0),
      "Delta_treatment1_3 - Delta_treatment2_3" =       c(0,0,0,0,0,0,1,0,0,0,0,0,1),
      "Delta_treatment2_1 - Delta_treatment2_2" =       c(0,0,0,0,0,0,0,1,0,0,0,1,0),
      "Delta_treatment2_1 - Delta_treatment2_3" =       c(0,0,0,0,0,0,0,0,1,0,0,0,1),
      "Delta_treatment2_2 - Delta_treatment2_3" =       c(0,0,0,0,0,0,0,-1,1,0,0,-1,1)
    )
  }else if (length(level_names)==4){
    contrast.matrix <-   rbind(
      "Delta_treatment1_1 - Delta_treatment1_2" =       c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
      "Delta_treatment1_1 - Delta_treatment1_3" =       c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),
      "Delta_treatment1_1 - Delta_treatment1_4" =       c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0),
      "Delta_treatment1_1 - Delta_treatment2_1" =       c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),
      "Delta_treatment1_1 - Delta_treatment2_2" =       c(0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0),
      "Delta_treatment1_1 - Delta_treatment2_3" =       c(0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0),
      "Delta_treatment1_1 - Delta_treatment2_4" =       c(0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1),
      "Delta_treatment1_2 - Delta_treatment1_3" =       c(0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0),
      "Delta_treatment1_2 - Delta_treatment1_4" =       c(0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0),
      "Delta_treatment1_2 - Delta_treatment2_1" =       c(0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0),
      "Delta_treatment1_2 - Delta_treatment2_2" =       c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0),
      "Delta_treatment1_2 - Delta_treatment2_3" =       c(0,0,0,0,0,0,0,1,-1,1,0,0,0,0,0,1,0),
      "Delta_treatment1_2 - Delta_treatment2_4" =       c(0,0,0,0,0,0,0,1,-1,0,1,0,0,0,0,0,1),
      "Delta_treatment1_3 - Delta_treatment1_4" =       c(0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0),
      "Delta_treatment1_3 - Delta_treatment2_1" =       c(0,0,0,0,0,0,0,1,0,-1,0,0,0,0,0,0,0),
      "Delta_treatment1_3 - Delta_treatment2_2" =       c(0,0,0,0,0,0,0,1,1,-1,0,0,0,0,1,0,0),
      "Delta_treatment1_3 - Delta_treatment2_3" =       c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0),
      "Delta_treatment1_3 - Delta_treatment2_4" =       c(0,0,0,0,0,0,0,1,0,-1,1,0,0,0,0,0,1),
      "Delta_treatment1_4 - Delta_treatment2_1" =       c(0,0,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0),
      "Delta_treatment1_4 - Delta_treatment2_2" =       c(0,0,0,0,0,0,0,1,1,0,-1,0,0,0,1,0,0),
      "Delta_treatment1_4 - Delta_treatment2_3" =       c(0,0,0,0,0,0,0,1,0,1,-1,0,0,0,0,1,0),
      "Delta_treatment1_4 - Delta_treatment2_4" =       c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1),
      "Delta_treatment2_1 - Delta_treatment2_2" =       c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0),
      "Delta_treatment2_1 - Delta_treatment2_3" =       c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0),
      "Delta_treatment2_1 - Delta_treatment2_4" =       c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1),
      "Delta_treatment2_2 - Delta_treatment2_3" =       c(0,0,0,0,0,0,0,0,-1,1,0,0,0,0,-1,1,0),
      "Delta_treatment2_2 - Delta_treatment2_4" =       c(0,0,0,0,0,0,0,0,-1,0,1,0,0,0,-1,0,1),
      "Delta_treatment2_3 - Delta_treatment2_4" =       c(0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,-1,1)
    )
  }
  ###then replace rownames matric with correct treatment
  for (treaty in 1:length(levels(dataset$treatment))){
    rownames(contrast.matrix) <- gsub(paste("treatment",treaty,sep=""),levels(dataset$treatment)[treaty],rownames(contrast.matrix))
  }
  ###next remove colony size (the last column) if not in model
  if (!"colony_size"%in%all.vars(form_stat)){ 
    new_contrast.matrix <- rbind(contrast.matrix[,-colony_index])
    rownames(new_contrast.matrix) <- rownames(contrast.matrix)
    contrast.matrix <- new_contrast.matrix
  }
  ###next replace names in contrast.matrix
  for (i in 1:length(level_names)){
    rownames(contrast.matrix) <- gsub(i,level_names[i],rownames(contrast.matrix))
  }
  
  ###next, if survival=T, remove intercept from the matrix
  if (survival){
    new_contrast.matrix <- rbind(contrast.matrix[,-1])
    rownames(new_contrast.matrix) <- rownames(contrast.matrix)
    contrast.matrix <- new_contrast.matrix
  }
  return(list(form_stat=form_stat,p_colony=p_colony,contrast.matrix=contrast.matrix))
}

prepare_stats_4 <- function(dataset,form_stat,survival){
  if (survival){statista <- "Pr(>|Chi|)"}else{statista <- "Pr(>Chisq)" } 
  interaction_problist <- vector("list", length(24))
  names(interaction_problist) <- as.character(sort(24))
  modellist <- vector("list", length(24))
  names(modellist) <- as.character(sort(24))
  for (time_window_ref in sort(24)){
    if (time_window_ref==unit){
      formi <- update(form_stat,.~.-(1|time_of_day))
    }else{
      formi <- form_stat
    }
    nb_time_windows <- floor(24/time_window_ref)
    interaction_probs <- NULL;
    models <- vector("list", nb_time_windows);names(models) <- as.character(1:nb_time_windows)
    for (time_point in 1:nb_time_windows){
      ###subset dataset
      subset <- dataset[(as.numeric(as.character(dataset$time_of_day))>=(time_window_ref*(time_point-1)))&(as.numeric(as.character(dataset$time_of_day))<(time_window_ref*time_point)),]
      if (exists("model")){rm=list=c("model")}
      gc()
      if (!survival){
        try(model <- do.call(lmer, list(formula=formi, data=subset)),silent=T)
        try(anov <- Anova(model),silent=T)
      }else{
        try(model <- do.call(coxme, list(formula=formi, data=subset)),silent=T)
        try(anov <- anova(model),silent=T)
      }
      
      if(exists("model")){
        if (time_window_ref==24){test_norm(residuals(model));print(anov)}
        models[[as.character(time_point)]] <- model
        rm(list=c("model"))
        gc()
        if(exists("pees")){rm(list=c("pees"))}
        try(pees <- anov[statista]["predictor",statista],silent=T)
        if (exists("pees")){
          p_interaction <- anov[statista]["time:treatment:predictor",statista]
          rm(list=c("pees"))
        }else{
          p_interaction <- NA
        }
      }else{
        models[[as.character(time_point)]] <- NULL
        p_interaction <- NA
      }
      names(p_interaction) <- time_point;interaction_probs <- c(interaction_probs,p_interaction)
    }
    interaction_probs <- p.adjust(interaction_probs,method="hochberg");
    interaction_problist[[as.character(time_window_ref)]] <- interaction_probs; 
    modellist[[as.character(time_window_ref)]] <- models
  }
  return(list(interaction_problist=interaction_problist,modellist=modellist))
}

printme_coxme <-function (x, rcoef = FALSE, digits = options()$digits, ...) 
{
  cat("Cox mixed-effects model fit by maximum likelihood\n")
  if (!is.null(x$call$data)) 
    #cat("  Data:", deparse(x$call$data))
    if (!is.null(x$call$subset)) {
      #cat(";  Subset:", deparse(x$call$subset), "\n")
    }
  else cat("\n")
  beta <- x$coefficients
  nvar <- length(beta)
  nfrail <- nrow(x$var) - nvar
  omit <- x$na.action
  cat("  events, n = ", x$n[1], ", ", x$n[2], sep = "")
  if (length(omit)) 
    cat(" (", naprint(omit), ")", sep = "")
  loglik <- x$loglik + c(0, 0, x$penalty)
  temp <- matrix(loglik, nrow = 1)
  cat("\n  Iterations=", x$iter, "\n")
  dimnames(temp) <- list("Log-likelihood", c("NULL", "Integrated", 
                                             "Fitted"))
  print(temp)
  cat("\n")
  chi1 <- 2 * diff(x$loglik[c(1, 2)])
  chi1 <- 2 * diff(loglik[1:2])
  chi2 <- 2 * diff(loglik[c(1, 3)])
  temp <- rbind(c(round(chi1, 2), round(x$df[1], 2), signif(1 - 
                                                              pchisq(chi1, x$df[1]), 5), round(chi1 - 2 * x$df[1], 
                                                                                               2), round(chi1 - log(x$n[1]) * x$df[1], 2)), c(round(chi2, 
                                                                                                                                                    2), round(x$df[2], 2), signif(1 - pchisq(chi2, x$df[2]), 
                                                                                                                                                                                  5), round(chi2 - 2 * x$df[2], 2), round(chi2 - log(x$n[1]) * 
                                                                                                                                                                                                                            x$df[2], 2)))
  dimnames(temp) <- list(c("Integrated loglik", " Penalized loglik"), 
                         c("Chisq", "df", "p", "AIC", "BIC"))
  print(temp, quote = F, digits = digits)
  cat("\nModel: ", deparse(x$call$formula), "\n")
  if (nvar > 0) {
    se <- sqrt(diag(x$var)[nfrail + 1:nvar])
    tmp <- cbind(beta, exp(beta), se, round(beta/se, 2), 
                 signif(1 - pchisq((beta/se)^2, 1), 2))
    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)", 
                                         "se(coef)", "z", "p"))
    to_return <- signif(1 - pchisq((beta/se)^2, 1), 2)
  }
  if (rcoef) {
    coef <- unlist(lapply(ranef(x), function(y) {
      if (is.matrix(y)) {
        z <- c(y)
        dd <- dimnames(y)
        names(z) <- c(outer(dd[[1]], dd[[2]], paste, 
                            sep = ":"))
        z
      }
      else y
    }))
    se <- sqrt(diag(x$var)[1:nfrail])
    rtmp <- cbind(coef, exp(coef), se)
    dimnames(rtmp) <- list(names(coef), c("coef", "exp(coef)", 
                                          "Penalized se"))
  }
  if (nvar > 0 && rcoef) {
    cat("Fixed and penalized coefficients\n")
    print(rbind(tmp, cbind(rtmp, NA, NA)), na.print = "", 
          digits = digits)
  }
  else if (rcoef) {
    cat("Penalized coefficients\n")
    print(rtmp, digits = digits)
  }
  else if (nvar > 0) {
    cat("Fixed coefficients\n")
    print(tmp, digits = digits)
  }
  cat("\nRandom effects\n")
  random <- VarCorr(x)
  nrow <- sapply(random, function(x) if (is.matrix(x)) 
    nrow(x)
    else length(x))
  maxcol <- max(sapply(random, function(x) if (is.matrix(x)) 1 + 
                         ncol(x) else 2))
  temp1 <- matrix(NA, nrow = sum(nrow), ncol = maxcol)
  indx <- 0
  for (term in random) {
    if (is.matrix(term)) {
      k <- nrow(term)
      nc <- ncol(term)
      for (j in 1:k) {
        temp1[j + indx, 1] <- sqrt(term[j, j])
        temp1[j + indx, 2] <- term[j, j]
        if (nc > j) {
          indx2 <- (j + 1):nc
          temp1[j + indx, 1 + indx2] <- term[j, indx2]
        }
      }
    }
    else {
      k <- length(term)
      temp1[1:k + indx, 1] <- sqrt(term)
      temp1[1:k + indx, 2] <- term
    }
    indx <- indx + k
  }
  indx <- cumsum(c(1, nrow))
  temp3 <- rep("", nrow(temp1))
  temp3[indx[-length(indx)]] <- names(random)
  xname <- unlist(lapply(random, function(x) if (is.matrix(x)) 
    dimnames(x)[[1]]
    else names(x)))
  temp <- cbind(temp3, xname, ifelse(is.na(temp1), "", format(temp1, 
                                                              digits = digits)))
  if (maxcol == 2) 
    temp4 <- c("Group", "Variable", "Std Dev", "Variance")
  else temp4 <- c("Group", "Variable", "Std Dev", "Variance", 
                  "Corr", rep("", maxcol - 3))
  dimnames(temp) <- list(rep("", nrow(temp)), temp4)
  print(temp, quote = F)
  invisible(x)
  return(to_return)
}

printcoef_coxme <-function (x, rcoef = FALSE, digits = options()$digits, ...) 
{
  cat("Cox mixed-effects model fit by maximum likelihood\n")
  if (!is.null(x$call$data)) 
    #cat("  Data:", deparse(x$call$data))
    if (!is.null(x$call$subset)) {
      #cat(";  Subset:", deparse(x$call$subset), "\n")
    }
  else cat("\n")
  beta <- x$coefficients
  nvar <- length(beta)
  nfrail <- nrow(x$var) - nvar
  omit <- x$na.action
  cat("  events, n = ", x$n[1], ", ", x$n[2], sep = "")
  if (length(omit)) 
    cat(" (", naprint(omit), ")", sep = "")
  loglik <- x$loglik + c(0, 0, x$penalty)
  temp <- matrix(loglik, nrow = 1)
  cat("\n  Iterations=", x$iter, "\n")
  dimnames(temp) <- list("Log-likelihood", c("NULL", "Integrated", 
                                             "Fitted"))
  print(temp)
  cat("\n")
  chi1 <- 2 * diff(x$loglik[c(1, 2)])
  chi1 <- 2 * diff(loglik[1:2])
  chi2 <- 2 * diff(loglik[c(1, 3)])
  temp <- rbind(c(round(chi1, 2), round(x$df[1], 2), signif(1 - 
                                                              pchisq(chi1, x$df[1]), 5), round(chi1 - 2 * x$df[1], 
                                                                                               2), round(chi1 - log(x$n[1]) * x$df[1], 2)), c(round(chi2, 
                                                                                                                                                    2), round(x$df[2], 2), signif(1 - pchisq(chi2, x$df[2]), 
                                                                                                                                                                                  5), round(chi2 - 2 * x$df[2], 2), round(chi2 - log(x$n[1]) * 
                                                                                                                                                                                                                            x$df[2], 2)))
  dimnames(temp) <- list(c("Integrated loglik", " Penalized loglik"), 
                         c("Chisq", "df", "p", "AIC", "BIC"))
  print(temp, quote = F, digits = digits)
  cat("\nModel: ", deparse(x$call$formula), "\n")
  if (nvar > 0) {
    se <- sqrt(diag(x$var)[nfrail + 1:nvar])
    tmp <- cbind(beta, exp(beta), se, round(beta/se, 2), 
                 signif(1 - pchisq((beta/se)^2, 1), 2))
    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)", 
                                         "se(coef)", "z", "p"))
    to_return <- tmp
  }
  if (rcoef) {
    coef <- unlist(lapply(ranef(x), function(y) {
      if (is.matrix(y)) {
        z <- c(y)
        dd <- dimnames(y)
        names(z) <- c(outer(dd[[1]], dd[[2]], paste, 
                            sep = ":"))
        z
      }
      else y
    }))
    se <- sqrt(diag(x$var)[1:nfrail])
    rtmp <- cbind(coef, exp(coef), se)
    dimnames(rtmp) <- list(names(coef), c("coef", "exp(coef)", 
                                          "Penalized se"))
    to_return <- rtmp
  }
  if (nvar > 0 && rcoef) {
    cat("Fixed and penalized coefficients\n")
    print(rbind(tmp, cbind(rtmp, NA, NA)), na.print = "", 
          digits = digits)
  }
  else if (rcoef) {
    cat("Penalized coefficients\n")
    print(rtmp, digits = digits)
  }
  else if (nvar > 0) {
    cat("Fixed coefficients\n")
    print(tmp, digits = digits)
  }
  cat("\nRandom effects\n")
  random <- VarCorr(x)
  nrow <- sapply(random, function(x) if (is.matrix(x)) 
    nrow(x)
    else length(x))
  maxcol <- max(sapply(random, function(x) if (is.matrix(x)) 1 + 
                         ncol(x) else 2))
  temp1 <- matrix(NA, nrow = sum(nrow), ncol = maxcol)
  indx <- 0
  for (term in random) {
    if (is.matrix(term)) {
      k <- nrow(term)
      nc <- ncol(term)
      for (j in 1:k) {
        temp1[j + indx, 1] <- sqrt(term[j, j])
        temp1[j + indx, 2] <- term[j, j]
        if (nc > j) {
          indx2 <- (j + 1):nc
          temp1[j + indx, 1 + indx2] <- term[j, indx2]
        }
      }
    }
    else {
      k <- length(term)
      temp1[1:k + indx, 1] <- sqrt(term)
      temp1[1:k + indx, 2] <- term
    }
    indx <- indx + k
  }
  indx <- cumsum(c(1, nrow))
  temp3 <- rep("", nrow(temp1))
  temp3[indx[-length(indx)]] <- names(random)
  xname <- unlist(lapply(random, function(x) if (is.matrix(x)) 
    dimnames(x)[[1]]
    else names(x)))
  temp <- cbind(temp3, xname, ifelse(is.na(temp1), "", format(temp1, 
                                                              digits = digits)))
  if (maxcol == 2) 
    temp4 <- c("Group", "Variable", "Std Dev", "Variance")
  else temp4 <- c("Group", "Variable", "Std Dev", "Variance", 
                  "Corr", rep("", maxcol - 3))
  dimnames(temp) <- list(rep("", nrow(temp)), temp4)
  print(temp, quote = F)
  invisible(x)
  return(to_return)
}

pwe <- function(data, timevar, deathvar, bounds) {
  # pwe: expands an S data frame for Piece-Wise Exponential survival
  # G. Rodriguez, Nov 29, 1992
  #
  # Check arguments: time and death must be variables in the data frame
  # and boundaries must be non-negative and strictly increasing
  if(!is.data.frame(data)) stop("First argument must be a data frame")
  if(is.na(match(tn <- deparse(substitute(timevar)), names(data))))
    stop(paste("\n\tSurvival time", tn, 
               "must be a variable in the data frame"))
  if(is.na(match(dn <- deparse(substitute(deathvar)), names(data))))
    stop(paste("\n\tDeath indicator", dn, 
               "must be a variable in the data frame"))
  width <- diff(bounds)
  if(any(bounds < 0) | any(width <= 0)) stop(paste(
    "Invalid interval boundaries in", deparse(substitute(
      bounds))))      #
  # Expand the data frame creating one pseudo-observation for each
  # interval visited, add interval number, events and exposure time
  # (existing variables with these names will be overwriten)
  n <- cut(data[, tn], bounds)
  data <- data[rep(seq(along = n), n),  ]
  i <- NULL
  for(k in 1:length(n))
    i <- c(i, 1:n[k])
  data$events <- ifelse(data[, tn] > bounds[i + 1], 0, data[, dn])
  data$exposure <- ifelse(data[, tn] > bounds[i + 1], width[i], data[, tn
                                                                     ] - bounds[i])
  data$interval <- i
  attr(data$interval, "levels") <- attr(n, "levels")
  data
}

quantil <- function(observed,random){
  below <- length(which(random<observed))/length(random)
  above <- length(which(random>observed))/length(random)
  quantil <- below + (1-below-above)/2
  return(quantil)
}

read.tag <- function(tagfile){
  tag <- read.table(tagfile,sep=",",comment.char="%",as.is=TRUE,fill=TRUE,stringsAsFactors=F)
  ##find out in which line the names are contained
  index_names <- match("#tag",tag[,1])
  ##remove the stuff above
  if (index_names > 1){
    tag <- data.frame(tag[index_names:nrow(tag),])
  }
  ##in case there was a problem with the reading, fix it
  if (ncol(tag)==1){
    ncols <- min(which(!is.na(as.numeric(as.character(tag[,1]))))) -1
    new_tag <- {}
    for (line in 1:(nrow(tag)/ncols)){
      new_tag <- rbind(new_tag,as.character(tag[(((line - 1) *ncols) +1 ):(((line - 1) *ncols) + ncols ),]))
    }
    tag <- data.frame(new_tag,stringsAsFactors=FALSE)
  }
  
  ##update the line in which the names are contained
  index_names <- match("#tag",tag[,1])
  ##get name list
  original_name_list <- as.character(tag[index_names,])
  names(tag) <- original_name_list
  tag <-tag[tag["#tag"]!="#tag",]
  names(tag)[names(tag)=="#tag"] <- "tag"
  row.names(tag) <- 1:nrow(tag)
  return(tag)
}

rgb2cmyk <- function(R,G,B){
  R <- R/255
  G <- G/255
  B <- B/255
  
  K <- 1-max(c(R, G, B))
  C <- (1-R-K) / (1-K)
  M <- (1-G-K) / (1-K)
  Y <- (1-B-K) / (1-K)
  
  outcol <- c(C,M,Y,K)
  names(outcol) <- c("C","M","Y","K")
  return(outcol)
}

scatterplot_violin_forpaper <- function(formula_stat,formula_plot,ylabel,xlabel,title,dat,ymin=NULL,ymax=NULL,xmin=NULL,xmax=NULL,sorting="status",time_point,IC=NULL,output=F,means=T,input_color=NULL,violin_params,point_cex=NULL,predict=NULL){
  violin_params <- as.numeric(unlist(violin_params))
  ##read violin param
  range <- violin_params[1]
  ylim_fac1 <- violin_params[2]
  ylim_fac2 <- violin_params[3]
  wex <- violin_params[4]
  h <- violin_params[5]
  
  if (all(dat$predictor==dat$colony_size)){
    formula_stat <- update(formula_stat,~.-colony_size)
  }
  dep_var <- row.names(attr(terms(formula_plot),"factors"))[1]
  pf <- parent.frame()
  dat["variable"] <- eval(parse(text=dep_var),dat,pf)
  
  ##########plotting########
  if (is.numeric(dat$predictor)){
    categories <- sort(unique(dat$predictor_plot))
  }else{
    categories <- unique(c(dat$predictor_plot))
    categories <- categories[order(match(categories,status_order))]
  }
  if (is.null(input_color)){
    if (sorting=="treatment"){
      colour_pal <- NULL
      for (category in categories){
        colour_pal <- c(colour_pal,get(paste(category,"_colour",sep="")))
      }
    }else{
      colour_pal <- rev(colorRampPalette(colour_palette_workers)(length(categories)))
    }
  }else{
    colour_pal <- rev(colorRampPalette(input_color)(length(categories)))
  }
  
  names(colour_pal) <- categories
  dat["colour"] <-  colour_pal[match(dat$predictor_plot,names(colour_pal))]
  forplot <- aggregate(na.rm=T,na.action="na.pass",colour~predictor_plot,FUN=unique,data=dat)
  forplot_med <- aggregate(na.rm=T,na.action="na.pass",variable~predictor_plot,FUN=median,data=dat);names(forplot_med)[names(forplot_med)=="variable"] <- "median"
  forplot_mean <- aggregate(na.rm=T,na.action="na.pass",variable~predictor_plot,FUN=mean,data=dat);names(forplot_mean)[names(forplot_mean)=="variable"] <- "mean"
  forplot <- merge(merge(forplot,forplot_med),forplot_mean)
  if (is.null(ymin)){
    ymin <- min(dat$variable,na.rm=T) - ylim_fac1*(max(dat$variable,na.rm=T)-min(dat$variable,na.rm=T))
    # ymin <- min(dat$variable,na.rm=T)
  }
  if (is.null(ymax)){
    ymax <- max(dat$variable,na.rm=T) + ylim_fac2*(max(dat$variable,na.rm=T)-min(dat$variable,na.rm=T))
    # ymax <- max(dat$variable,na.rm=T) 
  }
  
  ###prepare plot shell
  par_mar_ori <- par()$mar
  if(is.character(dat$predictor)){
    forplot ["pred"] <- as.numeric(factor(forplot$predictor_plot,levels=categories))/2
    dat["pred"] <- as.numeric(factor(dat$predictor_plot,levels=categories))/2
    par(bty="n",xaxt = "n")
  }else{
    forplot ["pred"] <- forplot$predictor_plot
    dat["pred"] <- dat$predictor_plot
    par(bty='l')
  }
  par(mar=par_mar_ori+c(0,0,0,0.5))
  values <- sort(unique(forplot$pred))
  if (is.null(xmin)){
    xlim <- c(min(forplot$pred),max(forplot$pred))+mean(diff(values,lag=1))*c(-0.5,0.5)
  }else{
    xlim <- c(xmin,xmax)
  }
  
  plot(dat$pred,dat$variable,xlab="",ylab="",xaxt="n",yaxt="n",cex.main=inter_cex,cex.lab=inter_cex,font.lab=1,cex.axis=min_cex,xlim=xlim,ylim=c(ymin,ymax),pch=21,type="n")
  if (!all(!grepl("Log",xlabel))){
    where <- axisTicks(c(par("usr")[1],par("usr")[2]),log=F)
    where <- where[which(where==round(where))]
    axis(1,at=where,labels=format(10^(where),scientific=T),cex.lab=inter_cex,cex.axis=min_cex)
    xlab <- substr(gsub("Log\\(","",xlabel),1,-1+nchar(gsub("Log\\(","",xlabel)))
    if (xlab=="Measured pathogen load"){
      title(xlab=expression(paste("Measured pathogen load (ng/", mu, "L)")),cex.lab=inter_cex,mgp=par()$mgp+c(0.2,0,0))
    }else{
      title(xlab=xlab,cex.lab=inter_cex,mgp=par()$mgp+c(0.2,0,0))
    }
  }else{
    axis(1,cex.lab=inter_cex,cex.axis=min_cex)
    title(xlab=xlabel,cex.lab=inter_cex)
  }
  
  if(all(!grepl("log",ylabel))){
    axis(2,cex.lab=inter_cex,cex.axis=min_cex)
    title(ylab=ylabel,cex.lab=inter_cex)
  }else{
    where <- axisTicks(c(par("usr")[3],par("usr")[4]),log=F)
    where <- where[which(where==round(where))]
    axis(2,at=where,labels=format(10^(where),scientific=T),cex.lab=inter_cex,cex.axis=min_cex)
    ylab <- as.character(ylabel[2])
    if (ylab=="Measured pathogen load"){
      title(ylab=expression(paste("Measured pathogen load (ng/", mu, "L)")),cex.lab=inter_cex,mgp=par()$mgp+c(0.1,0,0))
    }else{
      title(ylab=ylab,cex.lab=inter_cex,mgp=par()$mgp+c(0.1,0,0))
    }
    
    
  }
  
  if(is.character(dat$predictor)){
    par(xaxt = "s")
    axis(side=1,at=sort(unique(forplot$pred)),labels=full_statuses_names[categories],tick=F,lty=0,cex.axis=inter_cex)
    ###plus add grid lines
    par(xpd=F)
    # for (idx in 1:nrow(forplot)){
    #   abline(h=forplot[idx,"median"],lwd=line_min,lty=3,col="black")
    # }
    abline(h=0,lwd=line_min,lty=3,col="black")
  }
  if(!is.null(title)){title(main=title,cex.main=inter_cex,line=2.5)}
  ###if input, add 95%IC
  if(!is.null(IC)){
    polygon(c(par("usr")[c(1,2)],par("usr")[c(2,1)]),c(IC[1],IC[1],IC[2],IC[2]),border=NA,col=alpha("orange",0.1))
  }
  par_cex_ori <- par()$cex
  par(cex=0.3)
  if (is.null(point_cex)){
    cexMed <- min_cex
  }else{
    cexMed <- point_cex
  }
  ###add violins
  for (i in 1:nrow(forplot)){
    subset <- dat[which(dat$pred==forplot[i,"pred"]),"variable"]
    if (is.na(range)){
      VioPlot(na.omit(subset),col=alpha(forplot[i,"colour"],0.7), horizontal=F, at=forplot[i,"pred"], add=TRUE,lty=1, rectCol="black",wex=wex,border=NA,lwd=line_min,mode="Median",cexMed=cexMed)
    }else{
      VioPlot(na.omit(subset),range=range, h=h,col=alpha(forplot[i,"colour"],0.7), horizontal=F, at=forplot[i,"pred"], add=TRUE,lty=1, rectCol="black",wex=wex,border=NA,lwd=line_min,mode="Median",cexMed=cexMed)
    }
  }
  par(cex=par_cex_ori)
  
  
  ##########stats############
  ###initialise
  to_plot <- data.frame(intercept=as.numeric(as.character()),slope=as.numeric(as.character()),colour=as.character(),stringsAsFactors=F)
  to_mtext <- data.frame(effect=as.character(),pvalue=as.numeric(as.character()),stringsAsFactors=F)
  pred<- attr(terms(formula_stat),"term.labels")
  pred <- pred[!grepl("\\|",pred)]
  
  ###get final stats formula (step by step in case of a multiple terms model); i.e., reduce model
  try(model_temp <- lmer(formula_stat,data=dat),silent=T)
  if (exists("model_temp")){
    coeff <- Anova(model_temp,type="III")
    if ("Pr(>Chisq)"%in%colnames(coeff)){
      ###first test if colony size is significant; if not remove it
      if ("colony_size" %in%pred){
        if (coeff["colony_size","Pr(>Chisq)"]>0.05){
          formula_stat <- update(formula_stat,~.-colony_size)
          model_temp <- lmer(formula_stat,data=dat)
          coeff <- Anova(model_temp,type="III")
        }#if (coeff["colony_size","Pr(>Chisq)"]>0.05)
        ###and now that it is dealt with, remove it from pred
        pred <- pred[pred!="colony_size"]
      }#("colony_size" %in%pred)
      
      if (time_point=="comparison"){
        for (effect in pred){
          to_mtext <- rbind(to_mtext,data.frame(effect=effect,pvalue=coeff[effect,"Pr(>Chisq)"],stringsAsFactors=F))
        }
        interaction_effect <- pred[grepl("\\:",pred)]
        p_interaction <- coeff[interaction_effect,"Pr(>Chisq)"]
        if (p_interaction>0.05){
          if (sorting=="status"){
            formula_stat <- update(formula_stat,~.-predictor:status)
          }
          if (sorting=="period"){
            formula_stat <- update(formula_stat,~.-predictor:period)
          }
          pred<- attr(terms(formula_stat),"term.labels")
          pred <- pred[!grepl("\\|",pred)&!grepl("colony_size",pred)]
          try(model_bis <- lmer(formula_stat,data=dat),silent=T)
          if (exists("model_bis")){
            coeff <- Anova(model_bis,type="III")
            pvalue <-   coeff[pred[!grepl("predictor",pred)],"Pr(>Chisq)"]
            to_mtext[to_mtext$effect==pred[!grepl("predictor",pred)],"pvalue"] <- pvalue
            if (pvalue>0.05){
              if (sorting=="status"){
                formula_stat <- update(formula_stat,~.-status)
              }
              if (sorting=="period"){
                formula_stat <- update(formula_stat,~.-period)
              }
              pred<- attr(terms(formula_stat),"term.labels")
              pred <- pred[!grepl("\\|",pred)&!grepl("colony_size",pred)]
            }
          }
        }#if (p_interaction>0.05)
      }#if (time_point=="comparison")
    }#if ("Pr(>Chisq)"%in%colnames(coeff))
    
    rm(list=c("model_temp"))
  }#(exists("model_temp"))
  
  
  ####get final pvalues based on final model
  ########get names of the predictors in the table for extracting coefficients
  formula_simple <- update(formula_stat,~.-(1|colony)-(1|antid_1)-(1|antid_2)-(1|antid))
  pred2 <- Names(  formula_simple,dat);pred2 <- pred2[!grepl("Intercept",pred2)&!grepl("colony_size",pred2)];pred2 <- pred2[!grepl("\\|",pred2)]
  try(model_final <- lmer(formula_stat,data=dat),silent=T)
  #########testing normality of residuals
  resids <- residuals(model_final)
  test_norm(resids)
  
  if (exists("model_final")){
    coeff <- Anova(model_final,type="III")
    print(coeff)
    coeff2 <- summary(model_final)$coefficients
    if ("Pr(>Chisq)"%in%colnames(coeff)){
      if (length(pred)==1){
        pvalue <- coeff[pred,"Pr(>Chisq)"]
        to_output <- list(coeff2[pred2,"Estimate"])
        names(to_output) <- time_point
        if (nrow(to_mtext)==0){
          to_mtext <- rbind(to_mtext,data.frame(effect=pred,pvalue=pvalue,stringsAsFactors=F))
        }else{
          to_mtext[to_mtext$effect==pred,"pvalue"] <- pvalue
        }
        if ( is.numeric(dat$predictor)){
          if (pvalue < 0.05){
            if ("colony_size"%in%Names(formula_simple,dat)){
              to_plot <- rbind(to_plot,data.frame(
                intercept = coeff2["(Intercept)","Estimate"]+coeff2["colony_size","Estimate"]*mean(dat$colony_size),
                slope = coeff2[pred2,"Estimate"],
                #colour = statuses_colours[time_point],stringsAsFactors=F))
                colour = "black",stringsAsFactors=F))
              
            }else{
              to_plot <- rbind(to_plot,data.frame(intercept = coeff2["(Intercept)","Estimate"],slope = coeff2[pred2,"Estimate"],
                                                  #colour = statuses_colours[time_point],stringsAsFactors=F))
                                                  colour = "black",stringsAsFactors=F))
            }
          }#if (pvalue < 0.05)         
        }
      }#if (length(pred)==1)
    }#if ("Pr(>Chisq)"%in%colnames(coeff))
  }#if (exists("model_final"))
  
  ###plot ablines
  if (!is.null(predict)){
    predicted_value <- to_plot[1,"intercept"] + to_plot[1,"slope"]*predict
    segments(x0=predict,y0=ymin-0.1*(ymax-ymin),y1=predicted_value,col="springgreen2",xpd=F,lty=2)
    segments(x0=xlim[1]-0.1*(xlim[2]-xlim[1]),y0=predicted_value,x1=predict,col="red",xpd=F,lty=2)
  }
  par(xpd=F)
  if (nrow(to_plot)>=1){
    for (i in 1:nrow(to_plot)){
      abline(a=to_plot[i,"intercept"],b=to_plot[i,"slope"],col=to_plot[i,"colour"],lwd=line_inter)
    }# i
  }#(nrow(to_plot)>=1)
  ###plot mtext
  if (nrow(to_mtext)>=1){
    for (i in 1:nrow(to_mtext)){
      pvalue <- to_mtext[i,"pvalue"];effect <- to_mtext[i,"effect"]
      if(grepl("\\:",effect)){effect_ <- "Interaction: "}
      if(!grepl(sorting,effect)){
        effect_ <- paste(ylabel[1],": ",sep="")
        if (nchar(effect_)>30){
          effect_ <- "main: "
        }
      }
      if (sorting=="status"){
        if(!grepl("predictor",effect)){effect_ <- "treated vs. nestmates: "}
      }
      if (sorting=="period"){
        if(!grepl("predictor",effect)){effect_ <- "before vs. after: "}
      }
      p_value <- from_p_to_ptext(pvalue)
      if (pvalue>0.05){p_cex <- inter_cex}else{p_cex <- max_cex}
      par(xpd=T)
      if (nrow(to_mtext)==1){
        title(main=p_value,cex.main=p_cex,font.main=2,line=stat_line-((i-1)*0.75),xpd=T)
      }else{
        title(main=paste(effect_,p_value,sep=""),cex.main=p_cex,font.main=2,line=stat_line-((i-1)*0.75),xpd=T)
      }
      par(xpd=T)
    } # i
  }#if (nrow(to_mtext)>=1)
  if(output){return(to_output)}
  
  par(mar=par_mar_ori)
  if(!is.null(predict)){return(predicted_value)}else{return(predict)}
}#scatterplot_bubbles_qpcr

sqrt_transf <- function(x){
  if(all(x>=0)){
    replac_val <- 0
  }else{
    replac_val <- -min(x,na.rm=T)
  }
  return(sqrt(x+replac_val))
}

survival_analysis <- function(experiment,which_to_plot="all"){
  
  ####plotting function
  plotty <- function(data_table,variable,survimin,legend_title){
    data_table["variable"] <- data_table[,variable]
    
    ###########PART 1 - by variable ##########
    model1 <- coxme(Surv(time=survival,event=censored_status)~1+(1|colony),data=data_table)
    survival_model <- coxme(Surv(time=survival,event=censored_status)~1+variable+(1|colony),data=data_table)
    survival_model_plot <- survfit(Surv(time=survival,event=censored_status)~variable,data=data_table)
    
    status_names <- gsub("variable=","",unlist(strsplit(names(survival_model_plot$strata),split="\\."))[grepl("\\=",unlist(strsplit(names(survival_model_plot$strata),split="\\.")))])
    
    if (which_to_plot=="all"){
      xmin <- -0.1
    }else{
      xmin <- -0.5
    }
    surviplot <- plot(survival_model_plot, mark.time=F,col=alpha(statuses_colours[status_names],0),yaxt="n",xaxt="n",lty=styles[status_names],lwd=widths[status_names],ylab="Proportion surviving",xlab="Days after treatment",cex.axis=min_cex,cex.lab=inter_cex,xlim=c(xmin,1+1.1*max(data_table$survival)),bty="n",xaxs="i",yaxs="i",ylim=c(survimin,1+0.35*(1-survimin)))
    axis(side=2,at=seq(from=survimin,to=1,by=0.05),tick=T,cex.axis=min_cex,xaxs="i",yaxs="i")
    ats <- 1+c(2*c(0:(floor((max(data_table$survival)-1)/2))));
    ats <- c(ats,max(ats)+2)
    labs <- ats+1;
    if (!xmin %in% ats){ats <- c(xmin,ats);labs <- c("",labs)}
    
    axis(side=1,
         at = ats
         ,
         labels = labs
         ,
         tick=T,cex.axis=min_cex
    )
    lines(survival_model_plot, mark.time=F,col=NA,yaxt="n",xaxt="n",lty=styles[status_names],lwd=widths[status_names],xlab="",cex.lab=inter_cex,cex.axis=min_cex,xlim=c(0,1.1*max(data_table$survival)),bty="n",xaxs="i",yaxs="i")
    strata <- levels(summary(survival_model_plot)$strata)
    
    final_surv <- c()
    for (stratum in rev(strata)){
      survival_values <- c(1,summary(survival_model_plot)$surv[summary(survival_model_plot)$strata==stratum])
      time_values <- summary(survival_model_plot)$time[summary(survival_model_plot)$strata==stratum]
      std_error_values <- c(0,summary(survival_model_plot)$std.err[summary(survival_model_plot)$strata==stratum])
      last_y <- rev(summary(survival_model_plot)$surv[summary(survival_model_plot)$strata==stratum])[1]
      y_lo <- survival_values-std_error_values
      y_hi <- survival_values+std_error_values
      
      xvalues <- c(0,rep(time_values,each=2),max(survival_model_plot$time)*(1+1/96))
      yvalues <- rep(c(survival_values),each=2)
      y_lovalues <- rep(y_lo,each=2)
      y_hivalues <- rep(y_hi,each=2)
      
      polygon(x=c(xvalues,rev(xvalues))
              ,
              y=c(y_lovalues,rev(y_hivalues))
              ,
              col=alpha(statuses_colours[gsub("variable=","",stratum)],0.2)
              ,
              border=NA
      )
      lines(xvalues,yvalues, col=statuses_colours[gsub("variable=","",stratum)],lty=styles[gsub("variable=","",stratum)],lwd=widths[gsub("variable=","",stratum)])
      names(last_y) <- stratum
      final_surv <- c(final_surv,last_y)
    }
    if (grepl("load",status_names[1])){
      legend("topright",xjust=1,yjust=1,legend=gsub(" predicted"," simulated",full_statuses_names[status_names]),pt.bg=alpha(statuses_colours[status_names],0.2),col=alpha(statuses_colours[status_names],0.2),bty='n',pch=15,lty=0,lwd=widths[status_names],pt.lwd=1,pt.cex=2,text.col="white",cex=min_cex,title=legend_title,y.intersp=1.1)
      legend("topright",xjust=1,yjust=1,legend=gsub(" predicted"," simulated",full_statuses_names[status_names]),pt.bg=alpha(statuses_colours[status_names],0.2),col=statuses_colours[status_names],bty='n',pch=NA,lty=styles[status_names],lwd=widths[status_names],pt.lwd=1,pt.cex=2,text.col="black",cex=min_cex,title=legend_title,y.intersp=1.1)
    }else{
      legend("topright",xjust=1,yjust=1,legend=capitalize(gsub("Untreated","",gsub("\n","",full_statuses_names[status_names]))),pt.bg=alpha(statuses_colours[status_names],0.2),col=alpha(statuses_colours[status_names],0.2),bty='n',pch=15,lty=0,lwd=widths[status_names],pt.lwd=1,pt.cex=2,text.col="white",cex=min_cex,title=legend_title,y.intersp=1.1)
      legend("topright",xjust=1,yjust=1,legend=capitalize(gsub("Untreated","",gsub("\n","",full_statuses_names[status_names]))),pt.bg=alpha(statuses_colours[status_names],0.2),col=statuses_colours[status_names],bty='n',pch=NA,lty=styles[status_names],lwd=widths[status_names],pt.lwd=1,pt.cex=2,text.col="black",cex=min_cex,title=legend_title,y.intersp=1.1)
    }
    
    arrow_x <- 0
    arrow_y0 <- 1
    arrow_y1 <- 1+0.05*(1-survimin)
    arrows(arrow_x,arrow_y1,arrow_x,arrow_y0+(1-survimin)/100,length=0.025,lwd=line_max,col="springgreen3")
    
    segments(x0=arrow_x,y0=survimin,y1=arrow_y0-0.00075,lty=3,lwd=line_inter,col="springgreen3")
    par(lheight=.75)
    mtext("time of load\nprediction",side=3,las=3,at=arrow_x,line=-1.75,cex=par('cex')*inter_cex,adj=0.5,col="springgreen3")
    par(lheight=1)
    
    top_surv_fit <- survfit(Surv(time=survival,event=censored_status)~1,data=data_table[which(data_table$variable==gsub("variable=","",names(final_surv[final_surv==max(final_surv)]))),])
    arrow_x <- top_surv_fit$time[closest_match((break_point-1),top_surv_fit$time)]
    arrow_y0 <- top_surv_fit$surv[closest_match((break_point-1),top_surv_fit$time)-1]
    arrow_y1 <- 1+0.05*(1-survimin)
    arrows(arrow_x,arrow_y1,arrow_x,arrow_y0+(1-survimin)/100,length=0.025,lwd=line_max,col="red")
    segments(x0=arrow_x,y0=survimin,y1=arrow_y0-0.00075,lty=3,lwd=line_inter,col="red")
    par(lheight=.75)
    mtext("increase\nin mortality\n(untreated)",side=3,las=3,at=arrow_x,line=-1.75,cex=par('cex')*inter_cex,adj=0.5,col="red")
    par(lheight=1)
    
    ####comparison of survival before break point
    subsett <- data_table
    subsett[which(subsett$survival>=(break_point-1)),"censored_status"] <- 0
    subsett[which(subsett$survival>=(break_point-1)),"survival"] <- break_point-1
    model1 <- coxme(Surv(time=survival,event=censored_status)~1+(1|colony),data=subsett)
    survival_model <- coxme(Surv(time=survival,event=censored_status)~1+variable+(1|colony),data=subsett)
    print("Between load prediction and inflexion")
    print(anova(model1,survival_model))
    p_before_break <- anova(model1,survival_model)["P(>|Chi|)"][2,"P(>|Chi|)"]
    
    if (gsub("variable=","",names(final_surv[final_surv==max(final_surv)]))==status_names[length(status_names)]){
      hazard_before_break <- sprintf("%.2f", exp(-printcoef_coxme(survival_model)[paste("variable",gsub("variable=","",names(final_surv[final_surv==max(final_surv)])),sep=""),"coef"]))
    }else{
      hazard_before_break <- sprintf("%.2f", exp(printcoef_coxme(survival_model)[paste("variable",gsub("variable=","",names(final_surv[final_surv==max(final_surv)])),sep=""),"coef"]))
    }
    
    # ####comparison of survival after break point
    subsett <- data_table[which((!is.na(data_table$tag)&data_table$focal_status=="untreated")),]
    subsett <- subsett[which(subsett$survival>(break_point-1)),]
    subsett$survival <- subsett$survival-(break_point-1)
    
    print("Between inflexion and end")
    model1                    <- coxme(Surv(time=survival,event=censored_status)~1+(1|colony),data=subsett)
    survival_model            <- coxme(Surv(time=survival,event=censored_status)~1+variable+(1|colony),data=subsett)
    print( anova(model1,survival_model))
    p_after_break <- anova(model1,survival_model)["P(>|Chi|)"][2,"P(>|Chi|)"]
    if (gsub("variable=","",names(final_surv[final_surv==max(final_surv)]))==status_names[length(status_names)]){
      hazard_after_break <- sprintf("%.2f", exp(-printcoef_coxme(survival_model)[paste("variable",gsub("variable=","",names(final_surv[final_surv==max(final_surv)])),sep=""),"coef"]))
    }else{
      hazard_after_break <- sprintf("%.2f", exp(printcoef_coxme(survival_model)[paste("variable",gsub("variable=","",names(final_surv[final_surv==max(final_surv)])),sep=""),"coef"]))
    }
    if (p_before_break<0.06){
      text(x=0.5*(break_point-1),y=survimin+(0.25-1/20)*(1-survimin),labels=paste("HR=",hazard_before_break,sep=""),cex=min_cex,col="black")
      if (p_before_break<0.05){
        p_cex <- max_cex*1.1;adjust_line <- +0.00035; fonty <-  2
      }else{
        p_cex <- inter_cex;adjust_line <- 0;fonty <- 1
      }
    }else{
      p_cex <- inter_cex;adjust_line <- 0;fonty <- 1
    }
    text(x=0.5*(break_point-1),y=survimin+0.25*(1-survimin)+adjust_line,labels=gsub(" ","",paste(from_p_to_ptext(p_before_break),sep="")),cex=p_cex,font=fonty)
    if (p_after_break<0.06){
      text(x= (break_point-1)+0.5*(max(data_table$survival)-(break_point-1)),y=survimin+(0.25-1/20)*(1-survimin),labels=paste("HR=",hazard_after_break,sep=""),cex=min_cex,col="black")
      if (p_after_break<0.05){
        p_cex <- max_cex*1.1;adjust_line <- +0.00035; fonty <-  2
      }else{
        p_cex <- inter_cex;adjust_line <- 0;fonty <- 1
      }
      
    }else{
      p_cex <- inter_cex;adjust_line <- 0;fonty <- 1
    }
    text(x=(break_point-1)+0.5*(max(data_table$survival)-(break_point-1)),y=survimin+0.25*(1-survimin)+adjust_line,labels=gsub(" ","",paste(from_p_to_ptext(p_after_break),sep="")),cex=p_cex,font=fonty)
  }
  
  ###parameters
  root_path <- paste(disk_path,experiment,sep="/")
  break_point <- 3.8
  styles <- c(1,2,1,2,1,2); names(styles) <- c("treated","untreated","high_predicted_load","low_predicted_load","forager","nurse")
  widths <- c(2*line_max,2*line_max,2*line_max,2*line_max,2*line_max,2*line_max); names(widths) <- c("treated","untreated","high_predicted_load","low_predicted_load","forager","nurse")
  
  ###read survival data
  survival_data <- read.table(paste(root_path,"/original_data/survival.dat",sep=""),header=T,stringsAsFactors = F)
  
  ###get rid of accidental deaths / escapes
  survival_data <- survival_data[which(survival_data$censored_status%in%c("alive","dead")),]
  
  ###update censored column
  survival_data[survival_data$censored_status!="dead","censored_status"] <- 0
  survival_data[survival_data$censored_status=="dead","censored_status"] <- 1
  survival_data$censored_status <- as.numeric(survival_data$censored_status)
  survival_data$colony <- as.character(survival_data$colony)
  
  #######add predicted load
  predicted_load <- read.table(paste(root_path,"/transmission_simulations/post_treatment/experimentally_exposed_seeds/individual_simulation_results_observed.txt",sep=""),header=T,stringsAsFactors = F)
  survival_data <- merge(survival_data,predicted_load[which(predicted_load$period=="after"),c("colony","tag","simulated_load")],all.x=T,all.y=F,sort=F)
  survival_data["final_load"] <- NA
  survival_data[which(survival_data$simulated_load>high_threshold),"final_load"] <- "high_predicted_load"
  survival_data[which(survival_data$simulated_load<=high_threshold),"final_load"] <- "low_predicted_load"
  survival_data$final_load <- factor(survival_data$final_load)
  
  #######add task group
  task_group                                                               <- read.table(paste(disk_path,experiment,"original_data",task_group_file,sep="/"),header=T,stringsAsFactors = F)
  survival_data                                                            <- merge(survival_data,task_group,all.x=T,all.y=F)
  survival_data[which(survival_data$focal_status=="treated"),"task_group"] <- "treated"
  survival_data[which(survival_data$tag==queenid),"task_group"]            <- "queen"
  
  #######add interaction with treated
  int_with_treated                                   <- read.table(paste(disk_path,experiment,"processed_data/individual_behaviour/post_treatment/interactions_with_treated.txt",sep="/"),header=T,stringsAsFactors = F)
  survival_data                                      <- merge(survival_data,int_with_treated[c("colony","tag","duration_of_contact_with_treated_min")],all.x=T,all.y=F)
  survival_data$duration_of_contact_with_treated_min <- as.numeric(survival_data$duration_of_contact_with_treated_min)
  Q1                                                 <- quantile(survival_data[survival_data$focal_status!="treated","duration_of_contact_with_treated_min"],probs=1/3,na.rm=T)
  Q2                                                 <- quantile(survival_data[survival_data$focal_status!="treated","duration_of_contact_with_treated_min"],probs=2/3,na.rm=T)
  survival_data["frequency_contact_with_treated"]                                                               <- NA
  survival_data[which(survival_data$duration_of_contact_with_treated_min<Q1),"frequency_contact_with_treated"]  <- "low"
  survival_data[which(survival_data$duration_of_contact_with_treated_min>=Q1
                      &
                      survival_data$duration_of_contact_with_treated_min<Q2),"frequency_contact_with_treated"]  <- "medium"
  survival_data[which(survival_data$duration_of_contact_with_treated_min>=Q2),"frequency_contact_with_treated"] <- "high"
  
  ###start at -3 days
  nb_days_before         <- 3
  survival_data          <- survival_data[survival_data$survival >=(-nb_days_before*24),]
  survival_data$survival <- survival_data$survival+(nb_days_before*24)
  ###modify survival time from hours to days
  survival_data$survival <- as.numeric(survival_data$survival)/24
  
  ####fit a simple model to the data in order to predict expected death rate
  ####first get mortality rate before treatment for untreated
  survival_data_untreated <- survival_data[which(survival_data$focal_status=="untreated"),]
  modified_survival       <- survival_data_untreated
  modified_survival[which(modified_survival$survival>=nb_days_before),"censored_status"] <- 0
  modified_survival[which(modified_survival$survival>=nb_days_before),"survival"]        <- nb_days_before
  wfit <- survreg(Surv(time=survival,event=censored_status)~1,data=modified_survival)
  pp <- 0:1000/10000
  wsurv <- predict(wfit, type='quantile', p=pp, 
                   newdata=modified_survival[1:2,])

  expected_survival_at_nb_days_before    <- 1-pp[closest_match(nb_days_before,t(wsurv)[,1])]
  observed_survival_at_nb_days_before    <- nrow(modified_survival[which(modified_survival$censored_status==0),])/nrow(modified_survival)
  expected_survival_at_break_point       <- 1-pp[closest_match(nb_days_before+break_point,t(wsurv)[,1])]
  observed_survival_at_break_point       <- nrow(survival_data_untreated[which(survival_data_untreated$survival>=nb_days_before+break_point),])/nrow(survival_data_untreated)
  
  
  end_of_experiment_in_days              <- max(survival_data_untreated$survival,na.rm=T)
  expected_survival_end_of_experiment    <- 1-pp[closest_match(end_of_experiment_in_days,t(wsurv)[,1])]
  actual_survival_at_end_of_experiment   <- nrow(survival_data_untreated[which(survival_data_untreated$censored_status==0),])/nrow(survival_data_untreated)

  ###plot 1: survival = f(status)
  model1 <- coxme(Surv(time=survival,event=censored_status)~1+(1|colony),data=survival_data)
  survival_model <- coxme(Surv(time=survival,event=censored_status)~1+focal_status+(1|colony),data=survival_data)
  survival_model_plot <- survfit(Surv(time=survival,event=censored_status)~focal_status,data=survival_data)
  
  status_names <- gsub("focal_status=","",unlist(strsplit(names(survival_model_plot$strata),split="\\."))[grepl("\\=",unlist(strsplit(names(survival_model_plot$strata),split="\\.")))])
  
  p_exposure <- anova(model1,survival_model)["P(>|Chi|)"][2,"P(>|Chi|)"]
  print("Survival of treated vs. untreated workers")
  print(anova(model1,survival_model))
  hazard_exposure <- sprintf("%.2f", exp(-printcoef_coxme(survival_model)["focal_statusuntreated","coef"]))
  
  survival_model_plot_untreated <- survfit(Surv(time=survival,event=censored_status)~1,data=survival_data[which(survival_data$focal_status=="untreated"),])
  
  par_mar_ori <-par()$mar
  par(mar=par()$mar+c(0,0,-0.75,0))
  
  if (which_to_plot == "detailed"){
    surviplot <- plot(survival_model_plot, mark.time=F,col=alpha(statuses_colours[status_names],0),yaxt="n",xaxt="n",lty=styles[status_names],lwd=widths[status_names],ylab="Proportion surviving",xlab="Days after treatment",cex.axis=min_cex,cex.lab=inter_cex,xlim=c(0,1.1*max(survival_data$survival)),bty="n",xaxs="i",yaxs="i",ylim=c(0,1.35))
    axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1),tick=T,cex.axis=min_cex,xaxs="i",yaxs="i")
    ats <- c(nb_days_before - rev(2*c(1:floor(nb_days_before/2))),nb_days_before+2*c(0:ceiling((max(survival_data$survival)-nb_days_before)/2)))
    labs <- ats-nb_days_before
    if (!0 %in% ats){ats <- c(0,ats);labs <- c("",labs)}
    axis(side=1,
         at = ats
         ,
         labels = labs
         ,
         tick=T,cex.axis=min_cex
    )
    lines(survival_model_plot, mark.time=F,col=NA,yaxt="n",xaxt="n",lty=styles[status_names],lwd=widths[status_names],xlab="",cex.lab=inter_cex,cex.axis=min_cex,xlim=c(0,1.1*max(survival_data$survival)),bty="n",xaxs="i",yaxs="i")
    
    strata <- levels(summary(survival_model_plot)$strata)
    for (stratum in rev(strata)){
      survival_values <- c(1,summary(survival_model_plot)$surv[summary(survival_model_plot)$strata==stratum])
      time_values <- summary(survival_model_plot)$time[summary(survival_model_plot)$strata==stratum]
      std_error_values <- c(0,summary(survival_model_plot)$std.err[summary(survival_model_plot)$strata==stratum])
      last_y <- rev(summary(survival_model_plot)$surv[summary(survival_model_plot)$strata==stratum])[1]
      y_lo <- survival_values-std_error_values
      y_hi <- survival_values+std_error_values
      
      xvalues <- c(0,rep(time_values,each=2),max(survival_model_plot$time)*(1+1/96))
      yvalues <- rep(c(survival_values),each=2)
      y_lovalues <- rep(y_lo,each=2)
      y_hivalues <- rep(y_hi,each=2)
      
      polygon(x=c(xvalues,rev(xvalues))
              ,
              y=c(y_lovalues,rev(y_hivalues))
              ,
              col=alpha(statuses_colours[gsub("focal_status=","",stratum)],0.2)
              ,
              border=NA
      )
      lines(xvalues,yvalues, col=statuses_colours[gsub("focal_status=","",stratum)],lty=styles[gsub("focal_status=","",stratum)],lwd=widths[gsub("focal_status=","",stratum)])
      
    }
    
    legend("topright",xjust=1,yjust=1,legend=gsub("\\\nworkers"," ",gsub("\\\nforagers"," ",full_statuses_names[status_names])),pt.bg=alpha(statuses_colours[status_names],0.2),col=alpha(statuses_colours[status_names],0.2),bty='n',pch=15,lty=0,lwd=widths[status_names],pt.lwd=1,pt.cex=2,text.col="white",cex=min_cex,title="Workers",y.intersp=1.2)
    legend("topright",xjust=1,yjust=1,legend=gsub("\\\nworkers"," ",gsub("\\\nforagers"," ",full_statuses_names[status_names])),pt.bg=alpha(statuses_colours[status_names],0.2),col=statuses_colours[status_names],bty='n',pch=NA,lty=styles[status_names],lwd=widths[status_names],pt.lwd=1,pt.cex=2,text.col="black",cex=min_cex,title="Workers",y.intersp=1.2)
    
    treated_y <- surviplot$y[c(1)]
    untreated_y <- surviplot$y[c(2)]
    
    par(xpd=T)
    if (p_exposure<0.06){
      text(x= (break_point+3)+0.5*(max(survival_data$survival)-(break_point+3)),y=(0.55-1/20),labels=paste("HR=",hazard_exposure,sep=""),cex=min_cex,col="black")
      if (p_exposure<0.05){
        p_cex <- max_cex*1.1;adjust_line <- +0.00035; fonty <-  2
      }else{
        p_cex <- inter_cex;adjust_line <- 0;fonty <- 1
      }
      
    }else{
      p_cex <- inter_cex;adjust_line <- 0;fonty <- 1
    }
    text(x=(break_point+3)+0.5*(max(survival_data$survival)-(break_point+3)),y=0.55+adjust_line,labels=gsub(" ","",paste(from_p_to_ptext(p_exposure),sep="")),cex=p_cex,font=fonty)

    arrow_x <- nb_days_before
    arrow_y0 <- 1
    arrow_y1 <- 1.05
    arrows(arrow_x,arrow_y1,arrow_x,arrow_y0+0.015,length=0.025,lwd=line_max,col=statuses_colours["pathogen"])
    segments(x0=arrow_x,y0=0,y1=arrow_y0-0.005,lty=3,lwd=line_inter,col=statuses_colours["pathogen"])
    par(lheight=.75)
    mtext("pathogen\nexposure",side=3,las=3,at=nb_days_before,line=-1.75,cex=par('cex')*inter_cex,adj=0.5,col=statuses_colours["pathogen"])
    par(lheight=1)
    
    #####Now study that break point for untreated workers
    ####trying to appl  y example from http://data.princeton.edu/wws509/stata/recidivismR.html
    ##We create an id variable  
    survival_data$id <- as.character(interaction(survival_data$colony,survival_data$tag))
    
    ##define breaks
    breaks <- c(0,nb_days_before,nb_days_before+break_point,ceiling(max(survival_data$survival)/0.5)*0.5)
    ##now we use the pwefunction to create new dataset for piecewise
    survival_datax <- pwe(survival_data,survival,censored_status,breaks)
    survival_datax$interval = factor(survival_datax$interval,labels=levels(survival_datax$interval))
    
    ##We are now ready to fit a proportional hazards model with a piecewise exponential baseline where the hazard changes from year to year. 
    ##We use the same model as Wooldridge(2002), involving ten predictors, all fixed covariates. 
    fit=glm(events~interval+offset(log(exposure)),data=survival_datax[which(survival_datax$focal_status=="untreated"),], family=poisson)
    summary(fit)
    coeffs <- summary(fit)$coefficients
    b = coef(fit)
    h = exp( b[1] + c(0,b[2:length(b)]) )
    H = cumsum( diff(breaks)*h)
    S = exp(-H)
    
    arrow_x <- survival_model_plot_untreated$time[closest_match((nb_days_before+1),survival_model_plot_untreated$time)]
    arrow_y0 <- survival_model_plot_untreated$surv[closest_match((nb_days_before+1),survival_model_plot_untreated$time)]
    arrow_y1 <- 1.05
    arrows(arrow_x,arrow_y1,arrow_x,arrow_y0+0.015,length=0.025,lwd=line_max,col="springgreen3")
    segments(x0=arrow_x,y0=0,y1=arrow_y0-0.005,lty=3,lwd=line_inter,col="springgreen3")
    par(lheight=.75)
    mtext("time of load\nprediction",side=3,las=3,at=arrow_x,line=-1.75,cex=par('cex')*inter_cex,adj=0.5,col="springgreen3")
    par(lheight=1)
    
    arrow_x <- survival_model_plot_untreated$time[closest_match((nb_days_before+break_point),survival_model_plot_untreated$time)]
    arrow_y0 <- survival_model_plot_untreated$surv[closest_match((nb_days_before+break_point),survival_model_plot_untreated$time)]
    arrow_y1 <- 1.05
    segments(x0=arrow_x,y0=0,y1=arrow_y0-0.005,lty=3,lwd=line_inter,col="red")
    arrows(arrow_x,arrow_y1,arrow_x,arrow_y0+0.015,length=0.025,lwd=line_max,col="red")
    par(lheight=.75)
    mtext("increase\nin mortality\n(untreated)",side=3,las=3,at=arrow_x,line=-1.75,cex=par('cex')*inter_cex,adj=0.5,col="red")
    par(lheight=1)
  }
  
  if (which_to_plot == "second_only"){
    ###plot 2: survival of untreated = f(predicted load)
    #####start 1 day after treatment, once the predicted load is known
    survival_data <- survival_data[survival_data$survival >=nb_days_before+1,]
    survival_data$survival <- survival_data$survival-(nb_days_before+1)
    survival_data <- survival_data[which(!is.na(survival_data$final_load)),]
    #####and remove treated workers
    survival_data <- survival_data[which(survival_data$focal_status=="untreated"),]
    
    ###prepare plot
      model1 <- coxme(Surv(time=survival,event=censored_status)~1+(1|colony),data=survival_data)
    survival_model <- coxme(Surv(time=survival,event=censored_status)~1+final_load+(1|colony),data=survival_data)
    survival_model_plot <- survfit(Surv(time=survival,event=censored_status)~final_load,data=survival_data)
    
    status_names <- gsub("final_load=","",unlist(strsplit(names(survival_model_plot$strata),split="\\."))[grepl("\\=",unlist(strsplit(names(survival_model_plot$strata),split="\\.")))])
    
    p_load <- anova(model1,survival_model)["P(>|Chi|)"][2,"P(>|Chi|)"]
    hazard_load <- sprintf("%.2f", exp(-printcoef_coxme(survival_model)["final_loadlow_predicted_load","coef"]))
    
    if (which_to_plot=="all"){
      xmin <- -0.1
    }else{
      xmin <- -0.5
    }
    surviplot <- plot(survival_model_plot, mark.time=F,col=alpha(statuses_colours[status_names],0),yaxt="n",xaxt="n",lty=styles[status_names],lwd=widths[status_names],ylab="Proportion surviving",xlab="Days after treatment",cex.axis=min_cex,cex.lab=inter_cex,xlim=c(xmin,1+1.1*max(survival_data$survival)),bty="n",xaxs="i",yaxs="i",ylim=c(0.85,1.0525))
    axis(side=2,at=c(0.85,0.9,0.95,1),tick=T,cex.axis=min_cex,xaxs="i",yaxs="i")
    ats <- 1+c(2*c(0:(floor((max(survival_data$survival)-1)/2))));
    ats <- c(ats,max(ats)+2)
    labs <- ats+1;
    if (!xmin %in% ats){ats <- c(xmin,ats);labs <- c("",labs)}
    
    axis(side=1,
         at = ats
         ,
         labels = labs
         ,
         tick=T,cex.axis=min_cex
    )
    lines(survival_model_plot, mark.time=F,col=NA,yaxt="n",xaxt="n",lty=styles[status_names],lwd=widths[status_names],xlab="",cex.lab=inter_cex,cex.axis=min_cex,xlim=c(0,1.1*max(survival_data$survival)),bty="n",xaxs="i",yaxs="i")
    strata <- levels(summary(survival_model_plot)$strata)
    
    for (stratum in rev(strata)){
      survival_values <- c(1,summary(survival_model_plot)$surv[summary(survival_model_plot)$strata==stratum])
      time_values <- summary(survival_model_plot)$time[summary(survival_model_plot)$strata==stratum]
      std_error_values <- c(0,summary(survival_model_plot)$std.err[summary(survival_model_plot)$strata==stratum])
      last_y <- rev(summary(survival_model_plot)$surv[summary(survival_model_plot)$strata==stratum])[1]
      y_lo <- survival_values-std_error_values
      y_hi <- survival_values+std_error_values
      
      xvalues <- c(0,rep(time_values,each=2),max(survival_model_plot$time)*(1+1/96))
      yvalues <- rep(c(survival_values),each=2)
      y_lovalues <- rep(y_lo,each=2)
      y_hivalues <- rep(y_hi,each=2)
      
      polygon(x=c(xvalues,rev(xvalues))
              ,
              y=c(y_lovalues,rev(y_hivalues))
              ,
              col=alpha(statuses_colours[gsub("final_load=","",stratum)],0.2)
              ,
              border=NA
      )
      lines(xvalues,yvalues, col=statuses_colours[gsub("final_load=","",stratum)],lty=styles[gsub("final_load=","",stratum)],lwd=widths[gsub("final_load=","",stratum)])
      
    }
    
    legend("topright",xjust=1,yjust=1,legend=gsub(" predicted"," simulated",full_statuses_names[status_names]),pt.bg=alpha(statuses_colours[status_names],0.2),col=alpha(statuses_colours[status_names],0.2),bty='n',pch=15,lty=0,lwd=widths[status_names],pt.lwd=1,pt.cex=2,text.col="white",cex=min_cex,title="Untreated workers",y.intersp=1.1)
    legend("topright",xjust=1,yjust=1,legend=gsub(" predicted"," simulated",full_statuses_names[status_names]),pt.bg=alpha(statuses_colours[status_names],0.2),col=statuses_colours[status_names],bty='n',pch=NA,lty=styles[status_names],lwd=widths[status_names],pt.lwd=1,pt.cex=2,text.col="black",cex=min_cex,title="Untreated workers",y.intersp=1.1)
    
    treated_y <- surviplot$y[c(1)]
    untreated_y <- surviplot$y[c(2)]
    
    arrow_x <- 0
    arrow_y0 <- 1
    arrow_y1 <- 1.0075
    arrows(arrow_x,arrow_y0+0.03,arrow_x,arrow_y0+0.00225,length=0.05,lwd=line_max,col="springgreen3")
    segments(x0=arrow_x,y0=0.85,y1=arrow_y0-0.00075,lty=3,lwd=line_inter,col="springgreen3")
    if (which_to_plot=="all"){
      mtext("time of load\nprediction",side=3,las=3,at=arrow_x,line=-1.75,cex=par('cex')*inter_cex,adj=0.5,col="springgreen3")
    }
    
    arrow_x <- survival_model_plot_untreated$time[closest_match((break_point-1),survival_model_plot_untreated$time)]
    arrow_y0 <- survival_model_plot_untreated$surv[closest_match((break_point-1),survival_model_plot_untreated$time)]
    arrow_y1 <- 1.0075
    arrows(arrow_x,arrow_y0+0.03,arrow_x,arrow_y0+0.00225,length=0.05,lwd=line_max,col="red")
    segments(x0=arrow_x,y0=0.85,y1=arrow_y0-0.00075,lty=3,lwd=line_inter,col="red")
    if (which_to_plot=="all"){
      mtext("increase\nin mortality\n(untreated)",side=3,las=3,at=arrow_x,line=-1.75,cex=par('cex')*inter_cex,adj=0.5,col="red")
    }
    
    ####comparison of survival before break point
    subsett <- survival_data[which(!is.na(survival_data$tag)&survival_data$focal_status=="untreated"),]
    subsett[which(subsett$survival>=(break_point-1)),"censored_status"] <- 0
    subsett[which(subsett$survival>=(break_point-1)),"survival"] <- break_point-1
    model1 <- coxme(Surv(time=survival,event=censored_status)~1+(1|colony),data=subsett)
    survival_model <- coxme(Surv(time=survival,event=censored_status)~1+final_load+(1|colony),data=subsett)
    print("Between load prediction and inflexion")
    print(anova(model1,survival_model))
    p_before_break <- anova(model1,survival_model)["P(>|Chi|)"][2,"P(>|Chi|)"]
    hazard_before_break <- sprintf("%.2f", exp(-printcoef_coxme(survival_model)["final_loadlow_predicted_load","coef"]))
    
    ####comparison of survival after break point
    print("Between inflexion and end")
    subsett <- survival_data[which(!is.na(survival_data$tag)&survival_data$focal_status=="untreated"),]
    subsett <- subsett[which(subsett$survival>(break_point-1)),]
    subsett$survival <- subsett$survival-(break_point-1)
    
    model1 <- coxme(Surv(time=survival,event=censored_status)~1+(1|colony),data=subsett)
    survival_model <- coxme(Surv(time=survival,event=censored_status)~1+final_load+(1|colony),data=subsett)
    print( anova(model1,survival_model))
    p_after_break <- anova(model1,survival_model)["P(>|Chi|)"][2,"P(>|Chi|)"]
    hazard_after_break <- sprintf("%.2f", exp(-printcoef_coxme(survival_model)["final_loadlow_predicted_load","coef"]))
    
    if (p_before_break<0.05){
      text(x=0.5*(break_point-1),y=0.87,labels=paste("HR=",hazard_before_break,sep=""),cex=min_cex,col="black")
      p_cex <- max_cex*1.1;adjust_line <- +0.00035; fonty <-  2
    }else{
      p_cex <- inter_cex;adjust_line <- 0;fonty <- 1
    }
    text(x=0.5*(break_point-1),y=0.88+adjust_line,labels=gsub(" ","",paste(from_p_to_ptext(p_before_break),sep="")),cex=p_cex,font=fonty)
    if (p_after_break<0.05){
      text(x= (break_point-1)+0.5*(max(survival_data$survival)-(break_point-1)),y=0.87,labels=paste("HR=",hazard_after_break,sep=""),cex=min_cex,col="black")
      p_cex <- max_cex*1.1;adjust_line <- +0.00035; fonty <-  2
    }else{
      p_cex <- inter_cex;adjust_line <- 0;fonty <- 1
    }
    text(x=(break_point-1)+0.5*(max(survival_data$survival)-(break_point-1)),y=0.88+adjust_line,labels=gsub(" ","",paste(from_p_to_ptext(p_after_break),sep="")),cex=p_cex,font=fonty)
  }
  
  if (which_to_plot == "detailed"){
    #####start 1 day after treatment, once the predicted load is known
    survival_data <- survival_data[survival_data$survival >=nb_days_before+1,]
    survival_data$survival <- survival_data$survival-(nb_days_before+1)
    survival_data <- survival_data[which(!is.na(survival_data$final_load)),]
    #####and remove treated workers
    survival_data <- survival_data[which(survival_data$focal_status=="untreated"),]
    #####and remove queen 
    survival_data <- survival_data[which(survival_data$tag!=queenid),]
    
    #########create new tables that have limited data
    survival_data_nurses <- survival_data[which(survival_data$task_group=="nurse"),]
    survival_data_foragers <- survival_data[which(survival_data$task_group=="forager"),]
    
    survival_data_low_contact_with_treated    <- survival_data[which(survival_data$frequency_contact_with_treated=="low"),]
    survival_data_medium_contact_with_treated <- survival_data[which(survival_data$frequency_contact_with_treated=="medium"),]
    survival_data_high_contact_with_treated   <- survival_data[which(survival_data$frequency_contact_with_treated=="high"),]
    
    ####Plot 1: survival = f(task_group)
    plotty(  data_table=survival_data[!is.na(survival_data$task_group),], variable="task_group",survimin=0.75,legend_title="Untreated workers")
    
    ####Plot 2: survival=f(load), nurses
    plotty(  data_table=survival_data_nurses, variable="final_load",survimin=0.75,legend_title="Nurses")
    
    ####Plot 3: survival=f(load), low contacts
    plotty(  data_table=survival_data_low_contact_with_treated, variable="final_load",survimin=0.75,legend_title="Low contact with treated") 
  }
  
  par(mar=par_mar_ori)
}

test_norm <- function(resids){
  print("Testing normality")
  if (length(resids)<=300){
    print("Fewer than 300 data points so performing Shapiro-Wilk's test")
    print(shapiro.test(resids))
  }else{
    print("More than 300 data points so using the skewness and kurtosis approach")
    print("Skewness should be between -3 and +3 (best around zero")
    print(skewness(resids))
    print("")
    print("Excess kurtosis (i.e. absolute kurtosis -3) should be less than 4; ideally around zero")
    print(kurtosis(resids))
  }
}

transform_dataset <- function(dataset,cut=T,predictor,form_stat,excluded){
  combis <- unique(dataset[c("colony","time_hours","time_of_day")],stringsAsFactors=F); combis$time_of_day <- as.numeric(as.character(combis$time_of_day));
  if (0 %in%combis$time_hours){
    combis <- combis[which(combis$time_hours==0),]
  }else{
    combis <- combis[which(combis$time_hours==-24),]
  }
  ######1/ get most frequent starting time
  time_ori <- aggregate (colony~time_of_day,FUN=length,data=combis);time_ori <- as.numeric(as.character(time_ori[which(time_ori$colony==max(time_ori$colony,na.rm=T)),"time_of_day"]))
  ######2/ for colonies that started before, cut beginning - if only one day kept
  if(cut){
    if (kept_days==1){
      lows <- combis[combis$time_of_day<time_ori,"colony"]
      if (length(lows)>0){
        for (coli in lows){
          to_cut <- (combis[combis$colony==coli,"time_of_day"]):(time_ori-1)
          dataset[which((dataset$colony==coli)&(dataset$time=="After")&(dataset$time_of_day%in%to_cut)),"variable"] <- NA
        }
      }
      ######3/ for colonies that started after, cut 
      highs <- combis[combis$time_of_day>time_ori,"colony"]
      if (length(highs)>0){
        for (coli in highs){
          to_cut <- ((time_ori):(combis[combis$colony==coli,"time_of_day"]-1))
          dataset[which((dataset$colony==coli)&(dataset$time=="After")&(dataset$time_of_day%in%to_cut)),"variable"] <- NA
        }
      }
    }#kept_days
  }
  dataset["time_of_day_bis"] <- as.numeric(as.character(dataset$time_of_day))-time_ori
  dataset$time_of_day_bis[dataset$time_of_day_bis<0] <- dataset$time_of_day_bis[dataset$time_of_day_bis<0] + 24
  form_stat <- update(form_stat,.~.-(1|time_of_day)+(1|time_of_day_bis))
  if (!is.null(excluded)){
    combis["offset"] <- combis$time_of_day-time_ori
    excluded <- merge(excluded,combis[c("colony","offset")])
    excluded$time_of_day <- as.numeric(as.character(excluded$time_of_day))
    excluded$time_of_day <- excluded$time_of_day - excluded$offset
    excluded["time_of_day_bis"] <-  excluded$time_of_day - time_ori
  }
  return(list(dataset=dataset,form_stat=form_stat,excluded=excluded))
}

VioPlot <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
                     horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                     lwd = 1, rectCol = "black", colMed = "black", pchMed = 21, bgMed = "white",
                     at, add = FALSE, wex = 1, drawRect = TRUE, mode="Median",cexMed=min_cex){
  datas <- list(x, ...)
  n <- length(datas)
  if (missing(at)){ 
    at <- 1:n
  }
  std_low <- vector(mode = "numeric", length = n)
  std_high <- vector(mode = "numeric", length = n)
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  meaN <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))){
    args <- c(args, h = h)
  } 
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    meaN[i] <- mean(data)
    std <- std.error(data)
    std_low[i] <- meaN[i]-std
    std_high[i] <- meaN[i]+std
    
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        if (mode=="Median"){
          lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
                lty = lty)
          rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
               q3[i], col = rectCol)
          points(at[i], med[i], pch = pchMed, col = colMed, bg=bgMed, cex=cexMed)
        }else if (mode=="Mean"){
          # rect(at[i] - boxwidth/2, std_low[i], at[i] + boxwidth/2, 
          # std_high[i], col = rectCol)
          plot_arrows(means=data.frame(mean=meaN[i],se=std),plotx=at[i],plot_type="violinplot",LWD=line_max,LENGTH=0.05,colz="black")
          points(at[i], meaN[i], pch = pchMed, col = colMed, bg=bgMed, cex=cexMed)
        }
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}