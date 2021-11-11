#v1.8 update list
#- Allow pie charts on merged layers
#- hapid=2 means no labels at all
#- text.color added
#- Proportion option
#- Trapped error if no sequences assigned to a middle level
TempNet<-function(file,ftype="fasta",ages=NA,mut_size=1.5,nohap_size=.5,layernm=NA,mut.len=1,color=NA,useprop=F,
	invert=F,planes=F,vcol=F,theta=0,phi=pi/6,hapid=F,mergelayer=list(),confirm=T,save=T,load=F,text.col="black") { #v1.8.1
  #ages must be specified as 1,2,3,...; 1 will be on the bottom
  #hapid: 0 # in group, 1 haplogroup id, 2 [empty]
    library(pegas); library(ape); library(igraph)
    path<-function(pt,net,cp=pt[1]) {
        #recursive function for finding the path between any two haplogroups
        cx=tail(cp,1)
        op=c(net[net[,2]==cx,1],net[net[,1]==cx,2])
        if(pt[2] %in% op) return(c(cp,pt[2]))
        for(i in op) if(!i%in%cp) { 
            x=path(pt,net,c(cp,i))
            if(x[1]!=-1) return(x)
        }
        return(-1)
    }
    children<-function(pt,net,not) {
        ix=c(net[net[,1]==pt,2],net[net[,2]==pt,1]); ix=ix[!(ix%in%not)]
        if(length(ix)==0) return
        c(pt,unlist(lapply(ix,children,net,not=c(not,pt))))
    }
    z0plane<-function(x) c((x[1]-x[2]*sin(theta)/sin(phi))/cos(theta),-x[2]/sin(phi))
    if(!load) {
      read.dna(file,format=ftype)->temp
  	  nseq=0; 
      if (is.na(ages[1])) {
    		ag <- dimnames(temp)[[1]]
    		if(is.null(ag)) { ag<-names(temp); nseq=length(temp) } else nseq=nrow(temp)
    		options(warn=-1); ages <- as.integer(sub("^.+[$]","",ag)); options(warn=0)
    		if (all(is.na(ages))) ages <- rep(1,nseq)
    		if(min(ages)==0) ages=ages+1
    		missinglvl=which(!(1:max(ages)%in%ages))
    		if(length(missinglvl)) warning(paste("No sequences assigned to level",missinglvl,collapse="\n"))
    	} else { 
    		ages<-as.integer(factor(ages))
    		length(ages)->nseq
    	}
    	if(sum(is.na(ages))) warning("Could not assign all samples to a layer.")
    	if(nseq!=length(ages)) warning("Number of ages does not match number of sequences!")
      haplotype(temp)->s
      haploNet(s)->net
      #layout_with_kk(make_graph(c(t(net[,1:2])),directed=F),weights=net[,3])
      
      if(is.na(layernm[1])) layernm=paste("Layer",1:length(unique(ages)))
    	if(is.na(color[1])) color=rainbow(length(unique(ages)))
    	if(invert) { ages=max(ages)-ages+1; layernm<-rev(layernm) }
          
     	lvls=1:max(ages)
      if(length(text.col)<length(lvls)) text.col=rep(text.col,length.out=length(lvls))
      matrix(0,nr=nrow(s),nc=length(lvls))->lvl
      for(i in 1:nrow(s)) for(j in lvls) lvl[i,j]=sum(ages[attr(s,"index")[[i]]]==j) #how many of each hapgrp on each lvl
      if(useprop) { #adjust proportions so 1=1/max(n)
        perlvl=apply(lvl,2,sum)
        for(i in lvls) lvl[,i]=max(perlvl)/perlvl[i]*lvl[,i]
      }
      colmat=ifelse(lvl,col(lvl),0)
      npie=0; piepar=list() #haplogroups that need to be pie graphs
      if(length(mergelayer)) { #multi-colored layers; set mc=T, create color matrix
          todrop=c()
          for(i in 1:length(mergelayer)) {
              for(j in 1:nrow(lvl)) {
                  if(sum(lvl[j,mergelayer[[i]]]!=0)<2) colmat[j,mergelayer[[i]][1]]=max(colmat[j,mergelayer[[i]]])
                  else {
                      x=rep(0,ncol(lvl)); x[mergelayer[[i]]]=lvl[j,mergelayer[[i]]]
                      piepar[[(npie=npie+1)]]=x
                      colmat[j,mergelayer[[i]][1]]=-npie
                  }
              } 
              lvl[,mergelayer[[i]][1]]=apply(lvl[,mergelayer[[i]]],1,sum)
              todrop=c(todrop,mergelayer[[i]][-1])
          }
          lvl=lvl[,-todrop]; colmat=colmat[,-todrop]; mc=T
      } else mc=F
          
      table(c(net[,1:2]))->links
      
      #reorder net so the longest chain is roughly horizontal
      which.max(apply((x=combn(nrow(s),2)),2,function(x) sum(s[x[1],]!=s[x[2],])))->long
      x[,long]->long
      zp=path(long,net)
      lk=paste(net[,1],net[,2]) #for quicker searching
      newix=c(); curang=.5;tochk=F
      for(i in 2:length(zp)) {
          newix=c(newix,which(lk==paste(zp[i-1],zp[i]) | lk==paste(zp[i],zp[i-1])))
          if(tochk) {
              if((x=which(newix==tail(newix,1)))[1]!=length(newix)) newix=newix[-x[1]]
              else newix=newix[-(length(newix)-1)]
          }
          
          skip=which.min(abs((curang+.5+1:links[zp[i]]/links[zp[i]])%%1-.5))-1
          curang=curang+(skip+1)/links[zp[i]]-.5
          if(skip>0) {
              keyr=which(net[,1]==zp[i] | net[,2]==zp[i])
              keyr=keyr[keyr!=tail(newix,1)]
              newix=c(newix,keyr[0:skip+1]) #what if the "skip" is the one we need?
              tochk=T
          } else tochk=F
      }
      net=net[c(newix,(1:nrow(net))[-newix]),] #reorder
    }
    fnm=sub("[.][^.]+$",".tnf", file)
    if(load & file.exists(fnm)) elip=read.table(fnm)
    else {
      elip=cbind(layout_with_kk(make_graph(c(t(net[,1:2])),directed=F),weights=net[,3]),0)
      elip[,1:2]=elip[,1:2]*mut.len
    }
    #draw dat crazy graph
    par()$mar->oldmar
    clk=1; lbl.coords=NA
  	if(.Platform$OS.type=="mac") windows<-function(...) quartz(...)
  	if(.Platform$OS.type=="unix") windows<-function(...) X11(...)
  	windows(title="TempNet (press ESC when done)"); par(mar=c(0,0,0,0))
    while(clk[1]!=0) {
		range(c(elip[,2]-sqrt(lvl/pi),elip[,2]+sqrt(lvl/pi)))->rng
		diff(rng)/1.75->h #z-distance between layers
		range(elip[,1])+c(-max(sqrt(lvl[which.min(elip[,1]),]/pi)),max(sqrt(lvl[which.max(elip[,1]),]/pi)))->xrng
    in3d(expand.grid(xrng*(1+planes/10),c(rng[1],rng[2]),h*(1:ncol(lvl)-1)),theta,"pl",phi,type="n")
		if(is.na(lbl.coords[1])) lbl.coords=c(xrng[1],rng[1])
		else {
			if(lbl.coords[1]<par('usr')[1]) lbl.coords[1]=par('usr')[1]
			if(lbl.coords[2]<par('usr')[3]) lbl.coords[2]=par('usr')[3]
		}
    seq(0,2*pi,len=31)->seg #can't draw real elipses, so draw a regular 30-gon instead
    for(i in 1:ncol(lvl)) {
			if(planes) in3d(cbind(c(xrng,rev(xrng)),rep(rng,e=2),h*(i-1)),theta,"poly",phi,col=rgb(.2,.2,.2,a=.4))
			for(j in rev(order(elip[,2]))) { #elipses
				sqrt(lvl[j,i]/pi)->r; if(r==0) r=sqrt(nohap_size)/pi
        if(colmat[j,i]>=0) {
				  in3d(cbind(r*cos(seg)+elip[j,1],r*sin(seg)+elip[j,2],h*(i-1)),theta,"poly",phi,col=ifelse(lvl[j,i]==0,"white",color[colmat[j,i]]))
        } else {
          zz=piepar[[-colmat[j,i]]]; th=0
          for(k in 1:length(zz)) if(zz[k]) {
            xseg=seq(0,zz[k]/sum(zz)*2*pi,len=31)+th; th=th+zz[k]/sum(zz)*2*pi
            in3d(cbind(r*c(0,cos(xseg),0)+elip[j,1],r*c(0,sin(xseg),0)+elip[j,2],h*(i-1)),theta,"poly",phi,col=color[k])
          }
        }
			}
			for(j in rev(order(elip[,2]))) #vertical lines
    		if(i<ncol(lvl)) if(lvl[j,i+1]>0 && lvl[j,i]>0) #haplogroup exists at adjacent time points
    		  in3d(cbind(elip[j,1]+c(-1,1)*sqrt(lvl[j,i]/pi),rep(elip[j,2],2),h*(i-1),elip[j,1]+c(-1,1)*sqrt(lvl[j,i+1]/pi),rep(elip[j,2],2),h*i), theta,"seg",phi,col=ifelse(vcol,color[i+1],"black"))
      for(j in rev(order(elip[,2]))) if(lvl[j,i]>0) { 
        if(useprop) in3d(cbind(elip[j,1],elip[j,2],h*(i-1)),theta,"text",phi,labels=switch(hapid+1,paste0(round(100*lvl[j,i]/max(perlvl),1),"%"),j,""),col=text.col[i]) #n in haplogroup
        else in3d(cbind(elip[j,1],elip[j,2],h*(i-1)),theta,"text",phi,labels=switch(hapid+1,lvl[j,i],j,""),col=text.col[i]) #n in haplogroup
      }
        
      for(j in 1:nrow(net)) { #lines; solid if both haplogroups exist at that time, dashed if not
      	sqrt(lvl[net[j,1:2],i]/pi)->r; r[r==0]=sqrt(nohap_size)/pi
      	th=atan2(diff(elip[net[j,1:2],2]),diff(elip[net[j,1:2],1]))
        x=elip[net[j,1:2],1]+c(1,-1)*r*cos(th)
        y=elip[net[j,1:2],2]+c(1,-1)*r*sin(th)
				in3d(cbind(x,y,h*(i-1)),theta,"line",phi,lwd=1+(sum(lvl[net[j,1:2],i]==0)==0),lty=c("dotted","solid")[1+(sum(lvl[net[j,1:2],i]==0)==0)])
        if(net[j,3]>1) in3d(cbind(seq(from=x[1],to=x[2],len=net[j,3]+1)[-c(1,net[j,3]+1)], seq(from=y[1],to=y[2],len=net[j,3]+1)[-c(1,net[j,3]+1)],h*(i-1)),theta,"pt",phi,pch=20,cex=mut_size)
      }
			in3d(matrix(c(lbl.coords,h*(i-1)),nrow=1),theta,"text",phi,label=layernm[i],col=color[i],pos=4,cex=1.5) #layer labels
		}
		if(confirm==T) {
			spot=z0plane(unlist(locator(1)))		
			apply(elip[,1:2]-rep(spot,e=nrow(elip)),1,function(x) sum(x^2))->dst
			if(min(dst)<1) {
				clk=which.min(dst)
				cat("\nMoving haplogroup",clk)
                spot=z0plane(unlist(locator(1)))
				apply(elip[,1:2]-rep(spot,e=nrow(elip)),1,function(x) sum(x^2))->dst
                if(min(dst)<.6) { 
                    cat(" swap with haplogroup",(clk[2]=which.min(dst)))
                    if(clk[1]==clk[2]) { #move whole chain
                        ix=c(net[net[,1]==clk[1],2],net[net[,2]==clk[1],1])
                        if(length(ix)>1) { #need at least one more haplogroup
                            ix=lapply(ix,children,net,not=clk[1]); ixx=c(clk[1])
                            for(ji in (1:length(ix))[-which.max(lapply(ix,length))])
                                ixx=c(ixx,ix[[ji]])
                            in3d(cbind(elip[ixx,1:2],0),theta,drawcmd="pt",cex=4) #mark grp to move
                            spot=z0plane(unlist(locator(1)))
                            dx=spot-elip[clk[1],1:2]
                            elip[ixx,1:2]=elip[ixx,1:2]+rep(dx,e=length(ixx))
                        } 
                    } else {                    
                        ix=lapply(clk,function(x) c(net[net[,1]==x,2],net[net[,2]==x,1]))
                        pivot=ix[[1]][ix[[1]]%in%ix[[2]]]
                        if(length(pivot)==1) {
                            ix=lapply(clk,children,net,not=pivot)
                            th=diff(unlist(lapply(clk,function(x) atan2(elip[pivot,2]-elip[x,2],elip[pivot,1]-elip[x,1]))))
                            rotm=matrix(c(cos(th),sin(th),-sin(th),cos(th)),nr=2) #re-rotate
                            dx=elip[clk[1],1:2]
                            elip[ix[[1]],1:2]=t(elip[clk[2],1:2]+rotm%*%(t(matrix(elip[ix[[1]],1:2],nc=2))-elip[clk[1],1:2]))
                            rotm=matrix(c(cos(-th),sin(-th),-sin(-th),cos(-th)),nr=2) #re-rotate
                            elip[ix[[2]],1:2]=t(dx+rotm%*%(t(matrix(elip[ix[[2]],1:2],nc=2))-elip[clk[2],1:2]))
                        } else cat("\nERROR: Cannot swap!")
                    }
                    
                } else elip[clk,1:2]=spot
			} else { 
				if(sum(spot-lbl.coords)<1) {
                    cat("\n Moving labels"); lbl.coords=z0plane(unlist(locator(1)))
				} else { 
				    cat("\nAre you sure you want to quit (y/n)? ")
                    clk=readline()
                    if(substr(clk,1,1)=="y") break; 
				}
			}
		} else break
	}
        
	par(mar=oldmar)
	rep(NA,nrow(temp))->x
	names(x)=rownames(temp)
	attr(s,"index")->y
	for(i in 1:length(y)) x[y[[i]]]=i
    if(save) write.table(elip,fnm)
	list(x,lvl,net)
}

in3d<-function(pts,theta,drawcmd,phi=pi/6,...) { #pts must be a matrix where columns are x,y,z coords; returns 2d coords
	rotm<-matrix(c(cos(theta),-sin(theta),0,0,-sin(phi),cos(phi)),nrow=3)
	if(drawcmd=="seg") {
		pts[,1:2]<-pts[,1:3]%*%rotm
		pts[,3:4]<-pts[,4:6]%*%rotm
		pts[,1:4]->pts
	} else as.matrix(pts)%*%rotm->pts
	if(drawcmd=="pl") plot(pts,...)
	if(drawcmd=="pt") points(pts,...)
	if(drawcmd=="poly") polygon(pts,...)
	if(drawcmd=="seg") segments(pts[,1],pts[,2],pts[,3],pts[,4],...)
	if(drawcmd=="line") lines(pts,...)
	if(drawcmd=="text") text(pts,...)
	pts
}

