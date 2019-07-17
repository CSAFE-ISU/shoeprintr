#' An example of an input shoeprint
#' @format a data frame
"input_example"

#' An example of a reference shoeprint
#' @format a data frame
"reference_example"

#' @title Implementation of a smart sampling algorithm
#'
#' @description Define a smart sampling function which will ensure uniform samples so that we can get better results even with a lower ncross value. This process will divide the input circle area into hexbins and assign probabilities to the bins which will be correlated with the number of points in the bin. This will result in sampling the points from high density areas while still ensuring a more uniform distribution.
#' @name smart_sample
#' @param pts_xy The set of points (for example, the input circle of a shoeprint)
#' @param sample_bins The number of hexbins to sample
#' @param sample_size The number of points within each bin to sample
#' @param xbins Defines the grid of hexbins to divide the overall image (by default, 20x20)
#' @param seed The random seed for reproducing results
#' @return A data frame consisting of the sampled points from the original input circle
#' @export
#' @import hexbin
#' @examples
#' \dontrun{
#' data(input_example)
#'
#' smart_sample(input_example, sample_bins = NULL, sample_size = NULL) # All points from all bins
#' smart_sample(input_example, sample_bins = 30 & sample_size = NULL) # All points from 30 bins
#' smart_sample(input_example, sample_bins = NULL, sample_size = 30) # 30 points from all bins
#' smart_sample(input_example, sample_bins = 30, sample_size = 1) # one point from 30 selected bins
#' }
##############################################################################################################
smart_sample <- function(pts_xy, sample_bins = NULL, sample_size = NULL, xbins = 20, seed = 1) {
    ## Generate a hexbin object
    pts_hxb <- hexbin(pts_xy,xbins=xbins,IDs=TRUE)

    ## Assign cell ID to data points
    pts_xy$cell_id <- pts_hxb@cID

    ## Get number of points in every bin and sample bins to choose points from
    hxb_stats <- data.frame(cell_id=pts_hxb@cell,count=pts_hxb@count)
    set.seed(seed)
    candidate_hxbs <- sample(hxb_stats$cell_id,ifelse(is.null(sample_bins),length(pts_hxb@cell),min(nrow(hxb_stats),sample_bins)),prob=hxb_stats$count/sum(hxb_stats$count))

    ## Candidate points and final sample of points
    pts_xy <- pts_xy[pts_xy$cell_id %in% candidate_hxbs,]
    if(is.null(sample_size)) return(pts_xy[,-3])
    set.seed(seed)
    return(pts_xy[unlist(tapply(1:nrow(pts_xy),pts_xy$cell_id,function(x) sample(x,min(length(x),sample_size)))),-3])
}


##############################################################################################################
## Defining a vectorized subfunction to get adjacency list for a single vertice in product graph
##
## To compare ot to earlier loops, In this function i1 and j1 are fixed and the adjacency listed for all
## i2's and j2's inside this single function.
## This function use 'expand.grid' function to generate all possible i2,j2 pairs and a vectorize substraction
## to get the edge matches from a single vertice in product graph.
## Parameter b_a is combination of b and a where
## b is single row from the distl_in matrix. Means i1 is fixed
## a is single row from the distl_ref matrix. Means j1 is fixed
## The a length is also cut accordingly to reduce computations by so as to aoid calculations of upper triangulation
## of adjacency matrix
##############################################################################################################
get_edge_vertice <- function(b_a,eps,la,lb)
{
    grid<- expand.grid(unlist(b_a[1]),unlist(b_a[2]))
    edge_match <- which(abs(grid[,1]-grid[,2])<eps)
    if(length(edge_match)==0) return(NULL)
    return(edge_match)
}


##############################################################################################################
## Function to return stats on largest found clique
##############################################################################################################
#' @import vec2dtransf
#' @import sp
#' @import graphics
#' @import grDevices
#' @importFrom stats dist
#' @importFrom stats median
get_clique_stats <- function(clique,c_in,circle_in,c_ref,circle_ref,la,lb,plot=TRUE)
{
    if (length(clique) < 3) {
      cat("Cannot calculate clique, too few points\n")

      return(NA)
    }

    ## Subset clique points from input and refrence data
    cq_idx <- as.data.frame(t(sapply(clique,function(z) c(i=ifelse(z %% lb==0,(z %/% lb),(z %/% lb)+1),j=ifelse(z %% lb==0,lb,z %% lb)))))
    c_in_cq <- c_in[cq_idx$i,]
    c_ref_cq <- c_ref[cq_idx$j,]

    ## Calculate affine transformation parameters i.e Rotation and Translation
    affine_tfr_mat <- AffineTransformation(as.data.frame(cbind(c_in_cq,c_ref_cq)))
    calculateParameters(affine_tfr_mat)

    ## Apply affine transformation to input circle and convert classes to spatial points
    circle_in_sps <- applyTransformation(affine_tfr_mat,SpatialPoints(circle_in))
    c_in_cq_sps <- applyTransformation(affine_tfr_mat,SpatialPoints(c_in_cq))
    circle_ref_sps <- SpatialPoints(circle_ref)

    ## Prepare the match return statistics
    circle_in_ref_dist <- apply(spDists(circle_in_sps,circle_ref_sps),1,min)
    circle_in_dist_idx <- circle_in_sps@coords[circle_in_ref_dist<2,]
    circle_in_ref_dist_idx <- (circle_ref_sps[apply(spDists(circle_in_sps,circle_ref_sps),1,which.min)]@coords)[circle_in_ref_dist<2,]
    center_new <- apply(circle_in_ref_dist_idx,2,function(x) ((min(x) + max(x))/2))
    radius_new <- max(sqrt((circle_in_ref_dist_idx[,1]-center_new[1])^2+(circle_in_ref_dist_idx[,2]-center_new[2])^2))

    ## Rotation Angle
    R= matrix(affine_tfr_mat@parameters[-c(3,6)],2,byrow=FALSE)
    out = svd(R)
    s1=sign(out$u[1,1]*out$u[2,2])
    s2=sign(out$v[1,1]*out$v[2,2])
    out$u = out$u%*%diag(c(s1,1))
    out$v = out$v%*%diag(c(s2,1))
    out$d = out$d%*%diag(c(s1*s2,1))
    theta = atan(out$u[2,1]/out$u[1,1])*180/pi
    phi = atan(out$v[1,2]/out$v[1,1])*180/pi
    rot_angle = min(abs(theta+phi),180-abs(theta+phi))


    if(plot)
    {
        ## 2*2 Plot layout
        par(mfrow=c(2,2))

        ## Plot distance
        plot(as.matrix(dist(c_in_cq)),as.matrix(dist(c_ref_cq)),xlab="",ylab="",main="Corresponding Distances Plot")

        ## Plot matching points
        plot(c_ref_cq,col="red",pch=19,main="Matching Points")
        points(c_in_cq_sps,col="blue")

        ## Plot all points
        plot(	circle_ref,col="red",
              xlim=c(min(circle_ref$X,circle_in_sps@coords[,1]),max(circle_ref$X,circle_in_sps@coords[,1])),
              ylim=c(min(circle_ref$Y,circle_in_sps@coords[,2]),max(circle_ref$Y,circle_in_sps@coords[,2])),
              pch=19,main="All Points"
        )
        points(circle_in_sps,col="blue",pch=19)

        ## Plot Close points
        plot(	circle_in_ref_dist_idx,col="red",pch=19,
              xlim=c(min(circle_in_dist_idx[,1],circle_in_ref_dist_idx[,1]),max(circle_in_dist_idx[,1],circle_in_ref_dist_idx[,1])),
              ylim=c(min(circle_in_dist_idx[,2],circle_in_ref_dist_idx[,2]),max(circle_in_dist_idx[,2],circle_in_ref_dist_idx[,2])),
              main="Close Points"
        )
        points(circle_in_dist_idx,col="blue",pch=19)
    }
    gc()
    ## Return Statistics
    mylist <- list(clique_stats = data.frame(	clique_size=length(clique),
                                              rotation_angle = rot_angle,
                                              reference_overlap=nrow(circle_in_ref_dist_idx)/nrow(circle_ref),
                                              input_overlap=nrow(circle_in_dist_idx)/nrow(circle_in),
                                              med_dist_euc = round(median(apply((circle_in_ref_dist_idx-circle_in_dist_idx)^2,1,sum)),5),
                                              new_center_x=center_new[1],new_center_y=center_new[2],
                                              new_radius=radius_new),
                   affine = affine_tfr_mat)

    return(mylist)
}

get_os <- function() {
    if (grepl("darwin", version$os)) {
        return("mac64")
    } else if (grepl("linux", version$os)) {
        return("lin64")
    } else {
        return("win64")
    }
}

#' Function to perform boosted clique on two input circles
#'
#' @name boosted_clique
#' @param circle_in The input circle
#' @param circle_ref The reference circle
#' @param ncross_in_bins Number of bins in the input circle (See smart_sample)
#' @param xbins_in Number of bins along each axis in the hexbin grid for the input circle
#' @param ncross_in_bin_size Number of points to sample from each bin in the input circle
#' @param ncross_ref_bins Number of bins in the reference circle
#' @param xbins_ref Number of bins along each axis in the hexbin grid for the reference circle
#' @param ncross_ref_bin_size Number of points to sample from each bin in the reference circle
#' @param eps Distance tolerance for declaring an edge match
#' @param seed The random seed for reproducing results
#' @param num_cores The number of processor cores for parallel processing
#' @param plot If TRUE, produce a plot of the clique results
#' @param verbose If TRUE, print out the timing results for each portion of the algorithm
#' @param cl Optionally, a parallel cluster that has already been created
#' @return The statistics for the matching between the two circles
#' @export
#' @import parallel
#' @importFrom stats dist
boosted_clique <- function(circle_in, circle_ref, ncross_in_bins = 30, xbins_in = 20, ncross_in_bin_size = 1, ncross_ref_bins = NULL, xbins_ref = 30, ncross_ref_bin_size = NULL, eps = .75, seed = 1, num_cores = 8, plot = TRUE, verbose = FALSE, cl = NULL) {
    start_time_i <- Sys.time()
    if(verbose) cat("Preparing circles data for edge matching.\n")
    ##############################################################################################################
    ## Sample data using function smart_sample (Not in use as to keep the seed and results same to soyoung's code)
    ##############################################################################################################
    ## Smart sample input circle points
    c_in <- smart_sample(circle_in,sample_bins=min(ncross_in_bins,nrow(circle_in)),xbins=xbins_in,sample_size=ncross_in_bin_size,seed=seed)
    ## Smart sample reference circle points or full points
    c_ref <- smart_sample(circle_ref,sample_bins=min(ncross_ref_bins,nrow(circle_ref)),xbins=xbins_ref,sample_size=ncross_ref_bin_size,seed=seed)

    # set.seed(seed)
    # c_in <- as.matrix(circle_in[sample(1:nrow(circle_in),ncross_in_bins),])
    # c_ref <- as.matrix(circle_ref[sample(1:nrow(circle_ref),ifelse(is.null(ncross_ref_bins),nrow(circle_ref),ncross_ref_bins)),])

    ## Data dimesion for later reference within function
    la <- nrow(c_in)					## Number of minutiae in image 1
    lb <- nrow(c_ref)					## Number of minutiae in image 2
    lg <- la * lb						## Number of vertices in the "product graph"
    ##############################################################################################################

    ##############################################################################################################
    ## Prepare distance matrix for adjancey lists calculations
    ##############################################################################################################
    ## Pairwise distance matrix between vertices in input and reference data
    dist_in <- as.matrix(dist(c_in))
    dist_ref <- as.matrix(dist(c_ref))

    ## Fill NA's to avoid loops
    dist_in[row(dist_in)==col(dist_in)] <- NA
    dist_ref[row(dist_ref)==col(dist_ref)] <- NA

    ## Split matrix to list of columns so as to vectorize calculation of adjacency list
    distl_in <- as.list(as.data.frame(dist_in))
    distl_ref <- as.list(as.data.frame(dist_ref))
    dist_l_grid <- expand.grid(distl_ref,distl_in)

    ## Adding k1 verice 1 name to list so as to get edge match for k2<k1 only and thus avoiding values
    ## in upper triangulation matrix (Visualising edge match matrix lg*lg)
    dist_l_grid$node <- 1:nrow(dist_l_grid)
    dist_l_grid$Var2_n <- mapply(function(x,y) unlist(x)[1:unlist(y)],dist_l_grid[,2],as.list(((1:nrow(dist_l_grid)) %/% lb)+1))
    dist_l_grid <- dist_l_grid[,c("Var1","Var2_n","node")]
    rownames(dist_l_grid) <- 1:nrow(dist_l_grid)
    gc()
    ##############################################################################################################
    end_time_i <- Sys.time()
    if(verbose) cat(paste("Prepared circles data for edge matching. Took",round(difftime(end_time_i,start_time_i,units="secs")),"Seconds.\n\n\n"))





    start_time_al <- Sys.time()
    if(verbose) cat("Calculating adjacency list.\n")
    ##############################################################################################################
    ## Adjacency list prepasration (Vectorized Code with parallel processing [*nix style])
    ##############################################################################################################
    ## Split data to distribute work to cores
    dist_l_grid_core_assign <- split(dist_l_grid,floor(seq(0,((num_cores*4)-.01),length.out=nrow(dist_l_grid))))

    ## Windows
    if (!is.null(cl)) {
        final_ks <- do.call(c,parLapply(cl, dist_l_grid_core_assign,function(x,eps,la,lb) {gc();apply(x,1,get_edge_vertice,eps=eps,la=la,lb=lb)},eps=eps,la=la,lb=lb))
    ## Mac/Linux
    } else {
        final_ks <- do.call(c,mclapply(dist_l_grid_core_assign,function(x,eps,la,lb) {gc();apply(x,1,get_edge_vertice,eps=eps,la=la,lb=lb)},eps=eps,la=la,lb=lb,mc.cores=num_cores))
    }

    gc()
    ##############################################################################################################
    end_time_al <- Sys.time()
    if(verbose) cat(paste("Calculated adjacency list. Took",round(difftime(end_time_al,start_time_al,units="secs")),"Seconds.\n\n\n"))

    start_time_mc <- Sys.time()
    if(verbose) cat("Finding largest clique.\n")
    ##############################################################################################################
    ## Generate graph and get largest cliques
    ##############################################################################################################
    ## Prune edge list and write it to disk for parallel clique calculation
    fro_node_name <- rep(as.numeric(sapply(names(final_ks),function(x) strsplit(x,"\\.")[[1]][2],USE.NAMES=FALSE)),sapply(final_ks,length))
    to_node_name <- unlist(final_ks,use.names=FALSE)
    edge_list <- paste(fro_node_name,to_node_name," ")[to_node_name<fro_node_name]
    edge_list[length(edge_list)] <- gsub("  "," ",edge_list[length(edge_list)])
    edge_list_pmc_format <- c(	"%%MatrixMarket matrix coordinate pattern symmetric  ",
                               paste(nrow(dist_l_grid),nrow(dist_l_grid),length(edge_list)," "),
                               edge_list
    )
    mytempdir <- tempdir()
    writeLines(edge_list_pmc_format, file.path(mytempdir, "edge.mtx"))
    ext <- ifelse(get_os() == "win64", ".exe", "")
    clique_max <- system(paste0(system.file(package = "shoeprintr", "bin", get_os(), paste0("pmc", ext)), " -f ", file.path(mytempdir, "edge.mtx"), " -a 0"), intern = TRUE)
    clique_max <- as.numeric(unlist(strsplit(gsub("Maximum clique: ","",clique_max[grepl("^Maximum clique: ",clique_max)])," ")))
    ## Cleanup
    gc()
    ##############################################################################################################
    end_time_mc <- Sys.time()
    if(verbose) cat(paste("Calculated largest clique. Took",round(difftime(end_time_mc,start_time_mc,units="secs")),"Seconds.\n\n\n"))

    ##############################################################################################################
    ## Get and return clique stats
    ##############################################################################################################
    return(get_clique_stats(clique_max,c_in,circle_in,c_ref,circle_ref,la=la,lb=lb,plot=plot))
    ##############################################################################################################
}

#' Adjust it to (0,0) in the left lower corner
#'
#' @name leftcorner_cent
#' @param data The input data
#' @export
leftcorner_cent <- function(data) {
    return(apply((-1*data),2,function(x) x-min(x)))
}

##############################################################################################################
## Generate circle on given x ratio,y ratio and radius percentage
##############################################################################################################
get_circle <- function(data,x_rt,y_rt,rad_pct)
{
    cir_center <- apply(data,2,function(x) max(x)-min(x))*c(x_rt,y_rt)
    cir_rad <-  min(apply(data,2,function(x) ((max(x)-min(x))/2)*rad_pct))
    return(as.data.frame(data[apply(data,1,function(x,cir_center,cir_rad) sqrt(sum((x-cir_center)^2))<cir_rad,cir_center,cir_rad),]))
}

##############################################################################################################
## Generate circle on given x,y and radius percentage
##############################################################################################################
get_circle_fix <- function(x,y,rad_pct,data)
{
    cir_center <- c(x,y)
    if (rad_pct <= 1) {
        cir_rad <-  min(apply(data,2,function(x) ((max(x)-min(x))/2)*rad_pct))
    } else {
        cir_rad <- rad_pct
    }
    return(as.data.frame(data[apply(data,1,function(x,cir_center,cir_rad) sqrt(sum((x-cir_center)^2))<cir_rad,cir_center,cir_rad),]))
}

##############################################################################################################
## Generate circle on given x,y and radius
##############################################################################################################
get_circle_rad_fix <- function(x,y,cir_rad,data)
{
    cir_center <- c(x,y)
    return(as.data.frame(data[apply(data,1,function(x,cir_center,cir_rad) sqrt(sum((x-cir_center)^2))<cir_rad,cir_center,cir_rad),]))
}

#' @title Matches an input and a reference shoeprint
#'
#' @description Function to perform matching of footprints by matching all the candidate circles with input circles. This function performs initial matching between 3 circles on the input print and 27 candidate circles on the reference shoeprint. After the initial cliques, it select 2 best circles for each input circle and performs reinforcement matching on them. In reinforcement matching, full points on reference circle are considered instead of user defined. The function returns the final circle parameters and the statistics of the triangle formed.
#'
#' @name match_print
#' @param print_in The input print
#' @param print_ref The reference print
#' @param circles_input The input circles, specified as a matrix with columns (x, y, rad), or NULL to automatically generate
#' @param circles_reference The reference circles, specified as a matrix with columns (x, y, rad), or NULL to automatically generate
#' @param ncross_in_bins Number of bins in the input circle (See smart_sample)
#' @param xbins_in Number of bins along each axis in the hexbin grid for the input circle
#' @param ncross_in_bin_size Number of points to sample from each bin in the input circle
#' @param ncross_ref_bins Number of bins in the reference circle
#' @param xbins_ref Number of bins along each axis in the hexbin grid for the reference circle
#' @param ncross_ref_bin_size Number of points to sample from each bin in the reference circle
#' @param max_rotation_angle The maximum rotation angle, in degrees, for inclusion in the best circle matches
#' @param eps Distance tolerance for declaring an edge match
#' @param seed The random seed for reproducing results
#' @param num_cores The number of processor cores for parallel processing
#' @param plot If TRUE, produce a plot of the clique results
#' @param verbose If TRUE, print out the timing results for each portion of the algorithm
#' @return The statistics for the matching between the two circles
#' @import dplyr
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom stats dist
#' @export
#' @examples
#' \dontrun{
#' data(input_example)
#' data(reference_example)
#'
#' ## Transform all the points to have (0,0) at lower left corner.
#' print_in <- leftcorner_cent(reference_example)
#' print_ref <- leftcorner_cent(reference_example)
#'
#' ## Perform Print Match
#' print_stats <- match_print(print_in, print_ref,
#'                           ncross_in_bins = 30, xbins_in = 20, ncross_in_bin_size = 1,
#'                           ncross_ref_bins = NULL, xbins_ref = 30, ncross_ref_bin_size = NULL,
#'                           eps = .75, seed = 1, num_cores = parallel::detectCores(),
#'                           plot = TRUE, verbose = FALSE)
#' print_stats
#' }
match_print <- function(print_in, print_ref, circles_input = NULL, circles_reference = NULL, ncross_in_bins = 30, xbins_in = 20, ncross_in_bin_size = 1, ncross_ref_bins = NULL, xbins_ref = 30, ncross_ref_bin_size = NULL, max_rotation_angle = 365, eps = .75, seed = 1, num_cores = 8, plot = FALSE, verbose = FALSE) {
    ##############################################################################################################
    ## Cut circles on input print and reference print to perform matching
    ## 3 on Input circle and 27 on Reference circle
    ##############################################################################################################
    ## Generate 1 circles on input shoeprint
    if (is.null(circles_input)) {
        circles_dims <- apply(print_in,2,function(x) max(x)-min(x))
        circle_centers <- matrix(c(	.25,.8,			## 1/4 on x, 4/5 on y
                                    .25,.3,			## 1/4 on x, 1/5 on y
                                    .7,.6			## 3/4 on x, 1/2 on y
        ),3,byrow=TRUE,dimnames=list(NULL,c("x","y"))) %*% diag(circles_dims)
        circle_centers <- cbind(circle_centers,c(.4,.4,.4))
    } else {
        circle_centers <- circles_input
    }

    circles_in <- apply(circle_centers,1,function(x,data) get_circle_fix(x[1],x[2],x[3],data),data=print_in)

    ## Generate candidate circles on reference shoeprint
    if (is.null(circles_reference)) {
        rad_pct <- .5
        x_radius <- min(apply(print_ref,2,function(x) ((max(x)-min(x))/2)*rad_pct))
        x_ratios <- c(.25,.5,.75)
        y_ratios <- x_radius*(1:round((max(print_ref)-min(print_ref))/x_radius))/(max(print_ref)-min(print_ref))
        xy_ratios <- expand.grid(x_ratios,y_ratios)
        circles_ref <- apply(xy_ratios,1,function(xy_rt,rad_pct,data) get_circle(data,xy_rt[1],xy_rt[2],rad_pct),rad_pct=rad_pct,data=print_ref)
    } else {
        circles_ref <- apply(circles_reference,1,function(x,data) get_circle_fix(x[1],x[2],x[3],data),data=print_ref)
    }
    ##############################################################################################################
    ##############################################################################################################


    ##############################################################################################################
    ## Perform initial matching on circles
    ## 3 on Input circle and ~27 on Reference circle
    ##############################################################################################################
    match_grid <- expand.grid(circles_ref,circles_in)

    ## Use a cluster for parallel computation if running on windows
    cl <- NULL
    if (get_os() == "win64") cl <- makeCluster(num_cores, renice = 0)

    ## Loop over all circle to circle comparisons
    match_result <- lapply(1:nrow(match_grid), function(match_idx) {
      cat(paste("Matching circle pair",match_idx,"out of", nrow(match_grid), "\n"))

      ## Get the circles
      circle_in <- match_grid[match_idx,2][[1]]
      circle_ref <- match_grid[match_idx,1][[1]]

      ## Affine transformation requires 3 control points
      if (nrow(circle_ref) > 2 && nrow(circle_in) > 2) {
        myboost <- boosted_clique(circle_in, circle_ref,
                       ncross_in_bins=ncross_in_bins,xbins_in=xbins_in,ncross_in_bin_size=ncross_in_bin_size,
                       ncross_ref_bins=ncross_ref_bins,xbins_ref=xbins_ref,ncross_ref_bin_size=ncross_ref_bin_size,
                       eps=eps,seed=seed,num_cores=num_cores,plot=plot,verbose=verbose,cl=cl
        )

        return(myboost$clique_stats)
      } else {
        cat("Skipping circle pair", match_idx, "due to no points contained\n")

        return(NA)
      }
    })

    ##############################################################################################################
    ##############################################################################################################

    ##############################################################################################################
    ## Select best 2 circles against each of input circle for reinforcement learning
    ## Best criterion is whichever have largest input circle overlap
    ##############################################################################################################
    match_result <- do.call(rbind,match_result)
    match_result$circle1 <- rep(1:length(circles_in), each = length(circles_ref))
    match_result$circle2 <- rep(1:length(circles_ref), times = length(circles_in))

    ## Take only ones with a low enough rotation angle
    match_result_best <- match_result %>%
      filter(rotation_angle <= max_rotation_angle) %>%
      group_by(circle1) %>%
      arrange(desc(input_overlap)) %>%
      slice(1:2)

    #best_2_match <- tapply(match_result$input_overlap,rep(1:length(circles_in), each=length(circles_ref)),function(x) order(x,decreasing=TRUE)[1:2])
    #best_2_match <- mapply(function(x,y,z) x+(y*z),best_2_match,(1:length(circles_in))-1,MoreArgs=list(z=length(circles_ref)),SIMPLIFY=FALSE)
    #best_2_match_params <- match_result[unlist(best_2_match),c("new_center_x","new_center_y")]
    best_2_match_params <- match_result_best %>%
      ungroup() %>%
      select(new_center_x, new_center_y) %>%
      as.data.frame()

    circles_ref_reinf <- apply(best_2_match_params,1,function(x,rad_pct,data) get_circle_fix(x[1],x[2],rad_pct,data),rad_pct=.6,data=print_ref )
    circles_in_reinf <- rep(circles_in, length.out = nrow(match_result_best))
    ##############################################################################################################
    ##############################################################################################################

    ##############################################################################################################
    ## Perform reinforcement learning
    ##############################################################################################################
    match_result_reinf <- lapply(1:length(circles_in_reinf), function(match_idx_reinf) {
        cat(paste("Reinforcement matching circle pair",match_idx_reinf,"out of", length(circles_in_reinf), "\n"))
        myboost <- boosted_clique(circle_in=circles_in_reinf[[match_idx_reinf]],circle_ref=circles_ref_reinf[[match_idx_reinf]],
                       ncross_in_bins=ncross_in_bins,xbins_in=xbins_in,ncross_in_bin_size=ncross_in_bin_size,ncross_ref_bins=NULL,xbins_ref=30,ncross_ref_bin_size=NULL,
                       eps=eps,seed=seed,num_cores=num_cores,plot=plot,verbose=verbose,cl=cl
        )

        return(myboost$clique_stats)
    })
    if (!is.null(cl)) stopCluster(cl)

    ##############################################################################################################
    ## Prepare plot and  return output
    ##############################################################################################################
    match_result_reinf <- do.call(rbind,match_result_reinf)
    match_result_reinf$circle1 <- match_result_best$circle1
    match_result_reinf$circle2 <- match_result_best$circle2

    #best_match <- tapply(match_result_reinf$input_overlap,rep(1:length(circles_in), each=2),function(x) order(x,decreasing=TRUE)[1],simplify=FALSE)
    #best_match <- mapply(function(x,y,z) x+(y*z),best_match,(1:length(circles_in))-1,MoreArgs=list(z=2),SIMPLIFY=FALSE)
    #circles_match <- match_result_reinf[unlist(best_match),]
    #best_match_params <- match_result_reinf[unlist(best_match),c("new_center_x","new_center_y","new_radius")]

    reinf_result_best <- match_result_reinf %>%
      #filter(rotation_angle <= max_rotation_angle) %>%
      group_by(circle1) %>%
      arrange(desc(input_overlap)) %>%
      slice(1)

    circles_match <- reinf_result_best %>%
      ungroup() %>%
      as.data.frame()

    best_match_params <- circles_match %>%
      select(new_center_x, new_center_y, new_radius)

    circles_ref_out <- apply(best_match_params,1,function(x,data) get_circle_rad_fix(x[1],x[2],x[3],data),data=print_ref)

    circles_match <- cbind(circle_centers[,1:2, drop = FALSE],circles_match[,c("new_center_x","new_center_y","new_radius","rotation_angle","input_overlap")])
    names(circles_match) <- c("Fixed_circle_X","Fixed_circle_Y","Match_circle_X","Match_circle_Y","Match_circle_Radius","Rotation_angle","Input_circle_overlap_pct")

    ## Congurent Triangle output
    Fixed_triangle <- as.matrix(dist(circles_match[,c("Fixed_circle_X","Fixed_circle_Y")]))
    Match_Triangle <- as.matrix(dist(circles_match[,c("Match_circle_X","Match_circle_Y")]))

    p1 <- ggplot(data = as.data.frame(print_in), aes(x = x, y = y)) +
      geom_point() +
      geom_point(data = circles_in[[1]], color = "red") +
      theme_bw()
      if (length(circles_in) > 1) p1 <- p1 + geom_point(data = circles_in[[2]], color = "yellow")
      if (length(circles_in) > 2) p1 <- p1 + geom_point(data = circles_in[[3]], color = "green")

    p2 <- ggplot(data = as.data.frame(print_ref), aes(x = x, y = y)) +
      geom_point() +
      geom_point(data = circles_ref_out[[1]], color = "red") +
      theme_bw()
      if (length(circles_ref_out) > 1) p2 <- p2 + geom_point(data = circles_ref_out[[2]], color = "yellow")
      if (length(circles_ref_out) > 2) p2 <- p2 + geom_point(data = circles_ref_out[[3]], color = "green")

    ## Plot
    #dev.off()
    #par(mfrow=c(1,2))
    #plot(print_in,pch=19)
    #points(circles_in[[1]],col="red")
    #points(circles_in[[2]],col="yellow")
    #points(circles_in[[3]],col="green")
    #plot(print_ref,pch=19)
    #points(circles_ref_out[[1]],col="red")
    #points(circles_ref_out[[2]],col="yellow")
    #points(circles_ref_out[[3]],col="green")

    try(grid.arrange(p1, p2, ncol = 2))

    fix1 <- Fixed_triangle[row(Fixed_triangle)<col(Fixed_triangle)]
    fix2 <- Match_Triangle[row(Match_Triangle)<col(Match_Triangle)]
    if (length(fix1) == 0) fix1 <- NA
    if (length(fix2) == 0) fix2 <- NA

    old_result <- cbind( circles_match,
                         data.frame(	Triangle_sides = c("A-B","A-C","B-C")[1:nrow(circles_match)],
                                     Fixed_circle_side_length = fix1,
                                     Match_circle_side_length = fix2
                         )
    )

    new_result<- sum_result(old_result)

    ## Return Stats
    return(list(old_result, new_result))
}



sum_result<-function(data){

  R1<-mean(data[,7])
  R2<-mean(abs(data[,9]-data[,10]))
  R3<-sd(data[,6])

  Result<-c(R1,R2,R3)
  return(Result)

}

#' @title Matches a center circle input and a reference shoeprint
#'
#' @description Function to perform matching of center footprints by matching all the candidate circles with input circles. This function performs initial matching between 3 circles on the input print and 27 candidate circles on the reference shoeprint. After the initial cliques, it select 2 best circles for each input circle and performs reinforcement matching on them. In reinforcement matching, full points on reference circle are considered instead of user defined. The function returns the final circle parameters and the statistics of the triangle formed.
#'
#' @name centercircle_match
#' @param input The input print
#' @param reference The reference print
#' @param output The output obtained from the match_print function
#'
#' @export
centercircle_match<-function(input, reference, output){

  if (nrow(output) < 3) stop("Must have three circles in the output")

  #center 1 - X
  cx.in.1<-mean(c(output[1,1],output[2,1]))
  #center 1 - Y
  cy.in.1<-mean(c(output[1,2],output[2,2]))

  #center 2 - X
  cx.in.2<-mean(c(output[2,1],output[3,1]))
  #center 2 - Y
  cy.in.2<-mean(c(output[2,2],output[3,2]))

  #center 3 - X
  cx.in.3<-mean(c(output[3,1],output[1,1]))
  #center 3 - Y
  cy.in.3<-mean(c(output[3,2],output[1,2]))

  #center 4 - X
  cx.in.4<-mean(c(output[1,1],output[2,1],output[3,1]))
  #center 4 - Y
  cy.in.4<-mean(c(output[1,2],output[2,2],output[3,2]))

  #estimated circle location
  #center 1 - X
  cx.re.1<-mean(c(output[1,3],output[2,3]))
  #center 1 - Y
  cy.re.1<-mean(c(output[1,4],output[2,4]))

  #center 2 - X
  cx.re.2<-mean(c(output[2,3],output[3,3]))
  #center 2 - Y
  cy.re.2<-mean(c(output[2,4],output[3,4]))

  #center 3 - X
  cx.re.3<-mean(c(output[3,3],output[1,3]))
  #center 3 - Y
  cy.re.3<-mean(c(output[3,4],output[1,4]))

  #center 4 - X
  cx.re.4<-mean(c(output[1,3],output[2,3],output[3,3]))
  #center 4 - Y
  cy.re.4<-mean(c(output[1,4],output[2,4],output[3,4]))

  # step 2 : 4 of center circles matching

  input_circles <- matrix(c(cx.in.1, cy.in.1, 40), nrow = 1, ncol = 3)
  reference_circles <- matrix(c(cx.re.1, cy.re.1, 55), nrow = 1, ncol = 3)

  cc1 <- match_print(input, reference, circles_input = input_circles, circles_reference = reference_circles,
                     ncross_in_bins = 30, xbins_in = 20, ncross_in_bin_size = 1,
                     ncross_ref_bins = NULL, xbins_ref = 30, ncross_ref_bin_size = NULL,
                     eps = .75, seed = 1, num_cores = parallel::detectCores(),
                     plot = FALSE, verbose = TRUE)


  input_circles <- matrix(c(cx.in.2, cy.in.2, 40), nrow = 1, ncol = 3)
  reference_circles <- matrix(c(cx.re.2, cy.re.2, 55), nrow = 1, ncol = 3)

  cc2 <- match_print(input, reference, circles_input = input_circles, circles_reference = reference_circles,
                     ncross_in_bins = 30, xbins_in = 20, ncross_in_bin_size = 1,
                     ncross_ref_bins = NULL, xbins_ref = 30, ncross_ref_bin_size = NULL,
                     eps = .75, seed = 1, num_cores = parallel::detectCores(),
                     plot = FALSE, verbose = TRUE)

  input_circles <- matrix(c(cx.in.3, cy.in.3, 40), nrow = 1, ncol = 3)
  reference_circles <- matrix(c(cx.re.3, cy.re.3, 55), nrow = 1, ncol = 3)

  cc3 <- match_print(input, reference, circles_input = input_circles, circles_reference = reference_circles,
                     ncross_in_bins = 30, xbins_in = 20, ncross_in_bin_size = 1,
                     ncross_ref_bins = NULL, xbins_ref = 30, ncross_ref_bin_size = NULL,
                     eps = .75, seed = 1, num_cores = parallel::detectCores(),
                     plot = FALSE, verbose = TRUE)


  input_circles <- matrix(c(cx.in.4, cy.in.4, 40), nrow = 1, ncol = 3)
  reference_circles <- matrix(c(cx.re.4, cy.re.4, 55), nrow = 1, ncol = 3)

  cc4 <- match_print(input, reference, circles_input = input_circles, circles_reference = reference_circles,
                     ncross_in_bins = 30, xbins_in = 20, ncross_in_bin_size = 1,
                     ncross_ref_bins = NULL, xbins_ref = 30, ncross_ref_bin_size = NULL,
                     eps = .75, seed = 1, num_cores = parallel::detectCores(),
                     plot = FALSE, verbose = TRUE)


  comp<-c('1-2','2-3','3-1','1-2-3')
  re<-rbind(cc1[[1]][1,1:7],cc2[[1]][1,1:7],cc3[[1]][1,1:7],cc4[[1]][1,1:7])
  Result<-data.frame(comp,re)

  return(Result)
}


#' @title Starting plot to position three local areas of shoe Q
#'
#' @description Function to plot shoe Q and shoe K together. In shoe Q, three circles are colored.
#'
#' @name start_plot
#' @param input The input print (shoe Q)
#' @param reference The reference print (shoe K)
#' @param input_circle Centers and radius for three circles that we fix in the input print (shoe Q)
#'
#' @export

start_plot<-function(input, reference, input_circle){

  cx1<-input_circle[1,1]
  cx2<-input_circle[2,1]
  cx3<-input_circle[3,1]
  cy1<-input_circle[1,2]
  cy2<-input_circle[2,2]
  cy3<-input_circle[3,2]
  r1<-input_circle[1,3]
  r2<-input_circle[2,3]
  r3<-input_circle[3,3]

  P1<-ggplot(data.frame(input), aes(x=x, y=y))+ geom_point(data=data.frame(input), aes(x=x, y=y), color='black',size=0.1) +
    geom_point(data=data.frame(int_inside_center(data.frame(input), r1, nseg, cx1,cy1)),color="red",size=0.1)+
    gg_circle(r1, xc=cx1, yc=cy1, color="red") +
    geom_point(data=data.frame(int_inside_center(data.frame(input), r2, nseg, cx2,cy2)),color="orange",size=0.1)+
    gg_circle(r2, xc=cx2, yc=cy2, color="orange") +
    geom_point(data=data.frame(int_inside_center(data.frame(input), r3, nseg, cx3, cy3)),color="green",size=0.1)+
    gg_circle(r3, xc=cx3, yc=cy3, color="green")


  P2<-ggplot(data.frame(reference), aes(x=x, y=y))+ geom_point(data=data.frame(reference), aes(x=x, y=y), color='black',size=0.1)



  return(P1+P2)

}


#' @title Detecting points within the circle
#'
#' @description Function to find points within the circle
#'
#' @name int_inside_center
#' @param Data The input print
#' @param r radius of the circle
#' @param c1 x value of the circle
#' @param c2 y value of the circle
#'
#' @export

int_inside_center<-function(Data, r, nseg, c1,c2){
  nseg=360
  x <- c1+r*cos(seq(0,2*pi, length.out=nseg))
  y <- c2+r*sin(seq(0,2*pi, length.out=nseg))
  circle<-data.frame(x,y)
  shoe_cc<-subset(Data,  (c1-r)<Data[,1] & Data[,1]<(c1+r) & (c2-r)<Data[,2] & Data[,2]<(c2+r))
  intersect_points<-NULL
  CIR<-(shoe_cc[,1]-c1)^2+(shoe_cc[,2]-c2)^2
  shoe_cc<-cbind(shoe_cc, CIR)
  inside_cir<-subset(shoe_cc, shoe_cc[,3]<r^2)
  x<-inside_cir[,1]
  y<-inside_cir[,2]
  int_pts_cir<-cbind(x,y)
  return(unique(int_pts_cir))
}


#' @title Adding the circlular area
#'
#' @description Function to find the circular area
#'
#' @name gg_circle
#' @param r radius of the circle
#' @param xc x value of the circle
#' @param yc y value of the circle
#'
#' @export
#'
gg_circle<-function(r, xc, yc, color="black", fill=NA, ...) {
  x <- xc + r*cos(seq(0, pi, length.out=100))
  ymax <- yc + r*sin(seq(0, pi, length.out=100))
  ymin <- yc + r*sin(seq(0, -pi, length.out=100))
  annotate("ribbon", x=x, ymin=ymin, ymax=ymax, color=color, fill=fill, ...)
}



#' @title Find the initial circle position by quantile of x and y ranges
#'
#' @description Function to find the initial circle position by quantile of x and y ranges. Centers are found as (30%, 80%), (20%, 40%), (70%, 70%) with radius of 50 for three circles.
#'
#' @name initial_circle
#' @param input input impression
#'
#' @export
#'
#'
initial_circle<-function(input){
  circles_dims <- apply(input, 2, function(x) max(x) -  min(x))
  circle_centers1 <- matrix(c(0.3, 0.8, 0.2, 0.4, 0.7, 0.7), 3,
                            byrow = TRUE, dimnames = list(NULL, c("x", "y"))) %*% diag(circles_dims)

  circle_centers2 <- circle_centers1 + matrix(c(min(input[,1]), min(input[,1]), min(input[,1]), min(input[,2]),min(input[,2]),min(input[,2])), 3, byrow=FALSE)
  circle_centers3 <- cbind(circle_centers2, c(50, 50, 50))
  input_circles<-circle_centers3
  return(input_circles)
}



#' @title Cut x and y coordinate points less than 1% quantile and larger than 99% quantile
#'
#' @description Function to cut x and y coordinate points less than 1% quantile and larger than 99% quantile
#'
#' @name focus_data2
#' @param data data with x and y coordinate values
#'
#' @export
#'
#'
focus_data2<-function(data){
  data.new<-subset(data, quantile(data[,1], 0.01)<data[,1] &
                     data[,1]<quantile(data[,1], 0.99) &
                     quantile(data[,2], 0.01)<data[,2] &
                     data[,2]<quantile(data[,2], 0.99))
  data.new<-data.frame(data.new)
  names(data.new)<-c("x","y")
  return(data.new)
}

#' @title Reverse right side shoe to look like left side shoe
#'
#' @description Function to transform coordinate values to look like from the left side shoe
#'
#' @name rev_R
#' @param Data data with x and y coordinate values from right side of shoe
#'
#' @export
#'
rev_R<-function(Data){
  newx<-max(Data[,1])-Data[,1]
  rev_Data<-data.frame(newx, Data[,2])
  names(rev_Data)<-c("x", "y")

  return(rev_Data)
}




#' @title  Draw multiple plots
#'
#' @description Draw multiple plots
#'
#' @name multiplot
#'
#' @export
#'

multiplot<-function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#' @title  First plot before starting matching to confirm the circle position.
#'
#' @description Draw edge coordiates with designated circles. Circle will be colored in both input and reference
#'
#' @name step0_plot
#' @param input Input image for fixing three circles
#' @param reference Input image for finding correspondence
#' @param input_circles Three circles which will be fixed in the input
#'
#' @export
#'

step0_plot<-function(input, reference, input_circles){

  cx1<-input_circles[1,1]
  cx2<-input_circles[2,1]
  cx3<-input_circles[3,1]
  cy1<-input_circles[1,2]
  cy2<-input_circles[2,2]
  cy3<-input_circles[3,2]
  r1<-input_circles[1,3]
  r2<-input_circles[2,3]
  r3<-input_circles[3,3]

  P1<-ggplot(data.frame(input), aes(x=x, y=y))+ geom_point(data=data.frame(input), aes(x=x, y=y), color='black',size=0.1) +
    geom_point(data=data.frame(int_inside_center(data.frame(input), r1, nseg, cx1,cy1)),color="red",size=0.1)+
    gg_circle(r1, xc=cx1, yc=cy1, color="red") +
    geom_point(data=data.frame(int_inside_center(data.frame(input), r2, nseg, cx2,cy2)),color="orange",size=0.1)+
    gg_circle(r2, xc=cx2, yc=cy2, color="orange") +
    geom_point(data=data.frame(int_inside_center(data.frame(input), r3, nseg, cx3, cy3)),color="green",size=0.1)+
    gg_circle(r3, xc=cx3, yc=cy3, color="green")


  P2<-ggplot(data.frame(reference), aes(x=x, y=y))+ geom_point(data=data.frame(reference), aes(x=x, y=y), color='black',size=0.1) +
    geom_point(data=data.frame(int_inside_center(data.frame(reference), r1, nseg, cx1,cy1)),color="red",size=0.1)+
    gg_circle(r1, xc=cx1, yc=cy1, color="red") +
    geom_point(data=data.frame(int_inside_center(data.frame(reference), r2, nseg, cx2,cy2)),color="orange",size=0.1)+
    gg_circle(r2, xc=cx2, yc=cy2, color="orange") +
    geom_point(data=data.frame(int_inside_center(data.frame(reference), r3, nseg, cx3,cy3)),color="green",size=0.1)+
    gg_circle(r3, xc=cx3, yc=cy3, color="green")

  return(multiplot(P1,P2,cols=2))

}




#' @title  Match input and reference with given circle information in input. The searching area is confined.
#'
#' @description Find corresponding areas in reference for fixed circles in input. To reduce the time, it confines the area for candidate circles in reference image.
#' For the first circle, it confines the area into left toe, for the second circle, it searches only left bottom, for the third circle, it searches right toe area of reference shoe.
#'
#' @name match_print_subarea
#' @param input Input image for fixing three circles
#' @param reference Input image for finding correspondence
#' @param input_circles Three circles which will be fixed in the input
#' @param max_rotation_angle maximum rotation angle that we allow for comparisons.
#'
#' @export
#'

match_print_subarea<-function(input, reference, input_circles, max_rotation_angle){

  #original : c(0.3, 0.85, 0.25, 0.3, 0.75, 0.75)


  if (is.null(input_circles)) {
    circles_dims <- apply(input, 2, function(x) max(x) -  min(x))
    circle_centers1 <- matrix( c(0.3, 0.85, 0.25, 0.3, 0.75, 0.75), 3, byrow = TRUE, dimnames = list(NULL, c("x", "y"))) %*% diag(circles_dims)
    circle_centers2 <- circle_centers1 + matrix(c(min(input[,1]), min(input[,1]), min(input[,1]), min(input[,2]),min(input[,2]),min(input[,2])), 3, byrow=FALSE)
    circle_centers3 <- cbind(circle_centers2, c(50, 50, 50))
    input_circles<-circle_centers3
  }   else {
    input_circles <- input_circles
  }


  ref<-data.frame(reference)
  ref_len_y<-(max(reference[,2])-min(reference[,2]))
  ref_len_x<-(max(reference[,1])-min(reference[,1]))

  ref_top<-subset(ref, y>0.5*ref_len_y+min(reference[,2]))
  ref_bottom<-subset(ref, y<=0.5*ref_len_y+min(reference[,2]))


  location_ref<-list()

  location_ref[[1]]<-subset(ref_top, x<(ref_len_x/2+min(reference[,1]))) # ref_top_left
  location_ref[[2]]<-subset(ref_bottom, x<(ref_len_x/2+min(reference[,1]))) #ref_bottom_left
  location_ref[[3]]<-subset(ref_top, x>=(ref_len_x/2+min(reference[,1]))) #ref_top_right
  location_ref[[4]]<-subset(ref_bottom, x>=(ref_len_x/2+min(reference[,1]))) #ref_bottom_right

  nseg=360
  in.cx<-NULL
  in.cy<-NULL
  in.r<-NULL
  final.cx<-NULL
  final.cy<-NULL
  final.r<-NULL
  MM<-NULL
  FM<-NULL
  ref_loc<-NULL
  rd_score<-NULL
  WM<-NULL
  ref_loc<-NULL
  rd_score<-NULL
  wrong_mat<-NULL
  WW<-NULL
  wrong_set<-NULL

  K<-matrix(c(2,3,4,1,3,4,1,2,4,1,2,3), nrow=4, byrow=T)

  for ( k in 1:3){

    print(paste("circle",k,"matching"))
    in.cx[k]<-input_circles[k,1]
    in.cy[k]<-input_circles[k,2]
    in.r[k]<-input_circles[k,3]

    circle_in<-data.frame(int_inside_center(data.frame(input), in.r[k], nseg, in.cx[k], in.cy[k]))


    R<-NULL
    r_ref<-(input_circles[1,3]+15)

    ref_loc<-location_ref[[k]]

    ref_loc_len_x<-(max(ref_loc$x)-min(ref_loc$x))
    ref_loc_len_y<-(max(ref_loc$y)-min(ref_loc$y))


    if(min(ref_loc$y)+r_ref<max(ref_loc$y)-r_ref){
      ref_cdd_y<-seq(min(ref_loc$y)+r_ref, max(ref_loc$y), by = r_ref)} else {
        ref_cdd_y<-seq(min(ref_loc$y)+r_ref, max(ref_loc$y)+r_ref, by = r_ref)
      }

    ref_cdd_x<-(min(ref_loc$x)+max(ref_loc$x))/2


    for ( j in 1:length(ref_cdd_y)){
      circle_cdd_ref<-data.frame(int_inside_center(data.frame(ref), r_ref, nseg, ref_cdd_x, ref_cdd_y[j]))

      if(nrow(circle_cdd_ref)>(nrow(circle_in)*0.2)){
        M<-try(boosted_clique(circle_in, circle_cdd_ref, ncross_in_bins = 30, xbins_in = 20,
                              ncross_in_bin_size = 1, ncross_ref_bins = NULL, xbins_ref = 30,
                              ncross_ref_bin_size = NULL, eps = 0.75, seed = 1, num_cores = parallel::detectCores()-1,
                              plot = FALSE, verbose = FALSE, cl = NULL))

        if (sum(is.na(M[[1]]))<1) {
          M2<-M$clique_stats
          Mat<-paste0('mat',j)
          R<-rbind(R,cbind(Mat,M2))}
      }
    }



    R1<-subset(R,rotation_angle<max_rotation_angle)

    if (nrow(R1)!=0){
      new.cx<-R1[which.max(R1$input_overlap),7]
      new.cy<-R1[which.max(R1$input_overlap),8]
      new.r<-R1[which.max(R1$input_overlap),9]+15} else{

        new.cx<-R[which.max(R$input_overlap),7]
        new.cy<-R[which.max(R$input_overlap),8]
        new.r<-R[which.max(R$input_overlap),9]+15
      }


    circle_ref<-data.frame(int_inside_center(data.frame(reference), new.r, nseg, new.cx, new.cy))
    step2_mat<-boosted_clique(circle_in, circle_ref, ncross_in_bins = 30, xbins_in = 20,
                              ncross_in_bin_size = 1, ncross_ref_bins = NULL, xbins_ref = 30,
                              ncross_ref_bin_size = NULL, eps = 0.75, seed = 1, num_cores = parallel::detectCores()-1,
                              plot = FALSE, verbose = FALSE, cl = NULL)$clique_stats

    final.cx[k]<-step2_mat[,6]
    final.cy[k]<-step2_mat[,7]
    final.r[k]<-step2_mat[,8]

    MM<-rbind(MM,step2_mat)


    #wrong_set<-rbind(location_ref[[K[k,1]]],location_ref[[K[k,2]]],location_ref[[K[k,3]]] )

    #count<-0
    #repeat{
    #  w_idx<-sample(1:nrow(wrong_set),1)
    #  rd.cx<-as.numeric(wrong_set[w_idx,][1])
    #  rd.cy<-as.numeric(wrong_set[w_idx,][2])
    #  new.r<-55
    #  rd_circle_ref<-data.frame(int_inside_center(data.frame(reference), new.r, nseg, rd.cx, rd.cy))
    #  count<-count+1
    #  if(nrow(rd_circle_ref)>0.3*nrow(circle_in) & nrow(rd_circle_ref)<1.3*nrow(circle_in) | count==20){
    #    break}
    #}



    #wrong_mat<-boosted_clique(circle_in, rd_circle_ref, ncross_in_bins = 30, xbins_in = 20,
    #ncross_in_bin_size = 1, ncross_ref_bins = NULL, xbins_ref = 30,
    #ncross_ref_bin_size = NULL, eps = 0.75, seed = 1, num_cores = parallel::detectCores(),
    #plot = FALSE, verbose = FALSE, cl = NULL)$clique_stats


    #WW<-rbind(WW,wrong_mat)

  }

  Input_X<-input_circles[,1]
  Input_Y<-input_circles[,2]
  Comp<-c('1-2','1-3','2-3')
  d_in_1<-sqrt((Input_X[1]-Input_X[2])^2+(Input_Y[1]-Input_Y[2])^2)
  d_in_2<-sqrt((Input_X[1]-Input_X[3])^2+(Input_Y[1]-Input_Y[3])^2)
  d_in_3<-sqrt((Input_X[2]-Input_X[3])^2+(Input_Y[2]-Input_Y[3])^2)
  Euc_input_dist<-c(d_in_1,d_in_2,d_in_3)
  Reference_X<-MM[,6]
  Reference_Y<-MM[,7]
  d_ref_1<-sqrt((Reference_X[1]-Reference_X[2])^2+(Reference_Y[1]-Reference_Y[2])^2)
  d_ref_2<-sqrt((Reference_X[1]-Reference_X[3])^2+(Reference_Y[1]-Reference_Y[3])^2)
  d_ref_3<-sqrt((Reference_X[2]-Reference_X[3])^2+(Reference_Y[2]-Reference_Y[3])^2)
  Euc_ref_dist<-c(d_ref_1,d_ref_2,d_ref_3)
  Reference_radius<-MM[,8]

  #FM<-data.frame(Input_X,Input_Y,Reference_X,Reference_Y,Reference_radius,MM[,c(1:5)],Comp,Euc_input_dist,Euc_ref_dist, WW)

  FM<-data.frame(Input_X,Input_Y,Reference_X,Reference_Y,Reference_radius,MM[,c(1:5)],Comp,Euc_input_dist,Euc_ref_dist)


  P1<-ggplot(data.frame(input), aes(x=x, y=y))+ geom_point(data=data.frame(input), aes(x=x, y=y), color='black',size=0.1) +
    geom_point(data=data.frame(int_inside_center(data.frame(input), in.r[1], nseg, in.cx[1], in.cy[1])),color="red",size=0.1)+
    gg_circle(in.r[1], xc=in.cx[1], yc=in.cy[1], color="red") +
    geom_point(data=data.frame(int_inside_center(data.frame(input), in.r[2], nseg, in.cx[2], in.cy[2])),color="yellow",size=0.1)+
    gg_circle(in.r[2], xc=in.cx[2], yc=in.cy[2], color="yellow") +
    geom_point(data=data.frame(int_inside_center(data.frame(input), in.r[3], nseg, in.cx[3], in.cy[3])),color="green",size=0.1)+
    gg_circle(in.r[3], xc=in.cx[3], yc=in.cy[3], color="green")


  P2<-ggplot(data.frame(reference), aes(x=x, y=y))+ geom_point(data=data.frame(reference), aes(x=x, y=y), color='black',size=0.1) +
    geom_point(data=data.frame(int_inside_center(data.frame(reference), final.r[1], nseg, final.cx[1],final.cy[1])),color="red",size=0.1)+
    gg_circle(final.r[1], xc=final.cx[1], yc=final.cy[1], color="red") +
    geom_point(data=data.frame(int_inside_center(data.frame(reference), final.r[2], nseg, final.cx[2],final.cy[2])),color="yellow",size=0.1)+
    gg_circle(final.r[2], xc=final.cx[2], yc=final.cy[2], color="yellow") +
    geom_point(data=data.frame(int_inside_center(data.frame(reference), final.r[3], nseg, final.cx[3],final.cy[3])),color="green",size=0.1)+
    gg_circle(final.r[3], xc=final.cx[3], yc=final.cy[3], color="green")

  try(multiplot(P1, P2, cols=2))


  return(FM)
}


