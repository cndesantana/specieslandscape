#----------------------------------------------------------
#Neutral metacommunities on ice landscape dynamics 
#Calculating the components of the landscape (how many fragments and their sizes)
#Input ice cover (P/A matrices)
#Charles N. de Santana & Carlos J. Melian, EAWAG Sept 2013
#----------------------------------------------------------

library(spatstat)

files<-read.csv("./list_of_mat.dat",sep="\n",header=FALSE);
nfiles<-length(files$V1);

##Reading all matrix described in list_of_matrix file
for(i in 1:nfiles){
	inputfile<-as.character(files$V1)[i];
	date_sufix<-unlist(strsplit(unlist(strsplit(inputfile,'/'))[2],"[.]"))[1]
	basalmat<-read.csv(inputfile,sep=" ",header=FALSE);
	basalmat<-as.matrix(basalmat[,-length(basalmat[1,])]);#there is an extra column in each file

	for(threshold in 1:99){
	#	threshold<-15;#Percentual of sea-ice coverage to be considered as 'ice'
		mat<-basalmat;
		threshold_sufix<-as.character(sprintf("%02d",threshold));
		dir.create(paste("distancemat/thsld_",threshold_sufix,sep=""),showWarnings = FALSE);
		outputfile_compsizes<-paste("distancemat/thsld_",threshold_sufix,"/output_componentsizes_thsld_",threshold_sufix,".dat",sep="");
		if(!file.exists(outputfile_compsizes)) file.create(outputfile_compsizes,overwrite=FALSE);

	#identifying the grid-points with sea-ice-concentration higher than the threshold
		mat[which(mat>=(threshold/100.),arr.ind=TRUE)]<-1;
		mat[which(mat<(threshold/100.),arr.ind=TRUE)]<-0;
	#converting the matrix to 'spatstat' format
	#identifying the components (blocks) and their sizes
		conn<-connected(as.im(mat),background=0,method="C");
		nc<-length(levels(conn));
		components<-sort(as.numeric(summary(conn)$table),decreasing=TRUE);
	#Identifying the CENTROIDS of the components
		W<-tiles(tess(image=conn))
		centroid<-t(matrix(unlist(lapply(W,centroid.owin)),nrow=2));
	#Calculating the MATRIX OF DISTANCES among the components
		matdist<-as.matrix(dist(centroid));
	
	#OUTPUTS
	#Components size
		cat(paste(date_sufix," ",sep=""),file=outputfile_compsizes,append=TRUE);
		cat(components,file=outputfile_compsizes,append=TRUE);
		cat("\n",file=outputfile_compsizes,append=TRUE);
	
	#Centroids and matrix of distances
		outputfile_centroids<-paste("distancemat/thsld_",threshold_sufix,"/output_centroids_",date_sufix,"_thsld_",threshold_sufix,".dat",sep="");
		outputfile_matdist<-paste("distancemat/thsld_",threshold_sufix,"/output_matdist_",date_sufix,"_thsld_",threshold_sufix,".dat",sep="");
		outputfile_icemat<-paste("distancemat/thsld_",threshold_sufix,"/output_icemat_",date_sufix,"_thsld_",threshold_sufix,".mat",sep="");
		if(!file.exists(outputfile_centroids)) file.create(outputfile_centroids);
		if(!file.exists(outputfile_matdist)) file.create(outputfile_matdist);
		for(j in 1:nc){
			cat(centroid[j,],file=outputfile_centroids,append=TRUE);
			cat("\n",file=outputfile_centroids,append=TRUE);
			
			cat(matdist[j,],file=outputfile_matdist,append=TRUE);
			cat("\n",file=outputfile_matdist,append=TRUE);
		}

		if(!file.exists(outputfile_icemat)) file.create(outputfile_icemat);
		for(j in 1:length(basalmat[,1])){
			cat(mat[j,],file=outputfile_icemat,append=TRUE);
			cat("\n",file=outputfile_icemat,append=TRUE);
		}
	}
}
