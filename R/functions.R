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
    return(data.frame(	clique_size=length(clique),
                       rotation_angle = rot_angle,
                       refrence_overlap=nrow(circle_in_ref_dist_idx)/nrow(circle_ref),
                       input_overlap=nrow(circle_in_dist_idx)/nrow(circle_in),
                       med_dist_euc = round(median(apply((circle_in_ref_dist_idx-circle_in_dist_idx)^2,1,sum)),5),
                       new_center_x=center_new[1],new_center_y=center_new[2],
                       new_radius=radius_new
    ))
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
match_print <- function(print_in, print_ref, circles_input = NULL, circles_reference = NULL, ncross_in_bins = 30, xbins_in = 20, ncross_in_bin_size = 1, ncross_ref_bins = NULL, xbins_ref = 30, ncross_ref_bin_size = NULL, max_rotation_angle = 60, eps = .75, seed = 1, num_cores = 8, plot = FALSE, verbose = FALSE) {
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
            boosted_clique(circle_in, circle_ref,
                           ncross_in_bins=ncross_in_bins,xbins_in=xbins_in,ncross_in_bin_size=ncross_in_bin_size,
                           ncross_ref_bins=ncross_ref_bins,xbins_ref=xbins_ref,ncross_ref_bin_size=ncross_ref_bin_size,
                           eps=eps,seed=seed,num_cores=num_cores,plot=plot,verbose=verbose,cl=cl
            )
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
    circles_in_reinf <- rep(circles_in,each=2)
    ##############################################################################################################
    ##############################################################################################################

    ##############################################################################################################
    ## Perform reinforcement learning
    ##############################################################################################################
    match_result_reinf <- lapply(1:length(circles_in_reinf), function(match_idx_reinf) {
        cat(paste("Reinforcement matching circle pair",match_idx_reinf,"out of 6...\n"))
        boosted_clique(circle_in=circles_in_reinf[[match_idx_reinf]],circle_ref=circles_ref_reinf[[match_idx_reinf]],
                       ncross_in_bins=ncross_in_bins,xbins_in=xbins_in,ncross_in_bin_size=ncross_in_bin_size,ncross_ref_bins=NULL,xbins_ref=30,ncross_ref_bin_size=NULL,
                       eps=eps,seed=seed,num_cores=num_cores,plot=plot,verbose=verbose,cl=cl
        )
    })
    if (!is.null(cl)) stopCluster(cl)

    ##############################################################################################################
    ## Prepare plot and  return output
    ##############################################################################################################
    match_result_reinf <- do.call(rbind,match_result_reinf)
    best_match <- tapply(match_result_reinf$input_overlap,rep(1:length(circles_in), each=2),function(x) order(x,decreasing=TRUE)[1],simplify=FALSE)
    best_match <- mapply(function(x,y,z) x+(y*z),best_match,(1:length(circles_in))-1,MoreArgs=list(z=2),SIMPLIFY=FALSE)
    circles_match <- match_result_reinf[unlist(best_match),]
    best_match_params <- match_result_reinf[unlist(best_match),c("new_center_x","new_center_y","new_radius")]
    circles_ref_out <- apply(best_match_params,1,function(x,data) get_circle_rad_fix(x[1],x[2],x[3],data),data=print_ref)

    circles_match <- cbind(circle_centers[,1:2],circles_match[,c("new_center_x","new_center_y","new_radius","rotation_angle","input_overlap")])
    names(circles_match) <- c("Fixed_circle_X","Fixed_circle_Y","Match_circle_X","Match_circle_Y","Match_circle_Radius","Rotation_angle","Input_circle_overlap_pct")

    ## Congurent Triangle output
    Fixed_triangle <- as.matrix(dist(circles_match[,c("Fixed_circle_X","Fixed_circle_Y")]))
    Match_Triangle <- as.matrix(dist(circles_match[,c("Match_circle_X","Match_circle_Y")]))

    p1 <- ggplot(data = as.data.frame(print_in), aes(x = x, y = y)) +
      geom_point() +
      geom_point(data = circles_in[[1]], color = "red") +
      geom_point(data = circles_in[[2]], color = "yellow") +
      geom_point(data = circles_in[[3]], color = "green") +
      theme_bw()

    p2 <- ggplot(data = as.data.frame(print_ref), aes(x = x, y = y)) +
      geom_point() +
      geom_point(data = circles_ref_out[[1]], color = "red") +
      geom_point(data = circles_ref_out[[2]], color = "yellow") +
      geom_point(data = circles_ref_out[[3]], color = "green") +
      theme_bw()

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

    old_result <- cbind(	circles_match,
                         data.frame(	Triangle_sides = c("A-B","A-C","B-C"),
                                     Fixed_circle_side_length = Fixed_triangle[row(Fixed_triangle)<col(Fixed_triangle)],
                                     Match_circle_side_length = Match_Triangle[row(Match_Triangle)<col(Match_Triangle)]
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

