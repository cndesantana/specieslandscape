#----------------------------------------------------------
#Neutral metacommunities on ice landscape dynamics 
#Calculating the components of the landscape (how many fragments and their sizes)
#Input ice cover (P/A matrices)
#Charles N. de Santana & Carlos J. Melian, EAWAG Sept 2013
#----------------------------------------------------------

library(spatstat)

#Defining a test random 5x5 matrix
#nrows=15;ncols=15;
#mat<-matrix(sample(c(0,1),1,size=nrows*ncols),ncol=ncols,nrow=nrows)
files<-read.csv("./list_of_mat.dat",sep="\n",header=FALSE);
nfiles<-length(files$V1);


#for(threshold in 1:99){
threshold<-15;#Percentual of sea-ice coverage to be considered as 'ice'
threshold_sufix<-as.character(sprintf("%02d",threshold));

outputfile_compsizes<-paste("distancemat/output_componentsizes_thsld_",threshold_sufix,".dat",sep="");
file.create(outputfile_compsizes);
##Reading all matrix described in list_of_matrix file
for(i in 1:nfiles){
	inputfile<-as.character(files$V1)[i];
	date_sufix<-substr(inputfile,17,38);#sufix considering a filename like 'data/arctic_sic_t_00193_date_01_10_1994.mat'
	mat<-read.csv(inputfile,sep=" ",header=FALSE);
	mat<-as.matrix(mat[,-305]);#there is an extra column in each file
#identifying the grid-points with sea-ice-concentration higher than the threshold
#	mat[which(mat>=(threshold/100.),arr.ind=TRUE)]<-1;
#	mat[which(mat<(threshold/100.),arr.ind=TRUE)]<-0;
#converting the matrix to 'spatstat' format
#identifying the components (blocks) and their sizes
	immat<-as.im(mat);
	conn<-connected(immat,background=0,method="C");
	nc<-length(levels(conn));
	components<-sort(as.numeric(summary(conn)$table),decreasing=TRUE);
#Identifying the CENTROIDS of the components
	W<-tiles(tess(image=conn))
	centroid<-t(matrix(unlist(lapply(W,centroid.owin)),nrow=2));
#Calculating the MATRIX OF DISTANCES among the components
	matdist<-as.matrix(dist(centroid));

#OUTPUTS
	cat(components,file=outputfile_compsizes,append=TRUE);
	cat("\n",file=outputfile_compsizes,append=TRUE);

	outputfile_centroids<-paste("distancemat/output_centroids_",date_sufix,"_thsld_",threshold_sufix,".dat",sep="");
	file.create(outputfile_centroids)
	for(j in 1:nc){
		cat(centroid[j,],file=outputfile_centroids,append=TRUE);
		cat("\n",file=outputfile_centroids,append=TRUE);
		
		outputfile_matdist<-paste("distancemat/output_matdist_",date_sufix,"_thsld_",threshold_sufix,".dat",sep="");
		cat(matdist[j,],file=outputfile_matdist,append=TRUE);
		cat("\n",file=outputfile_matdist,append=TRUE);
	}
}
#}
