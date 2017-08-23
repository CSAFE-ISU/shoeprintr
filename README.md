Shoeprint Boosted
=====================================

How-To
------
1. Make sure you have the `devtools` package installed
2. From R, run `devtools::install_github("CSAFE-ISU/shoeprintr")`
3. (If running on the CSAFE server or another machine that restricts executing binary files) From a Terminal, run `chmod +x /home/username/R/x86_64-pc-linux-gnu-library/3.4/shoeprintr/bin/lin64/pmc` replacing the path as necessary
4. From R, Load `shoeprintr` by typing `library(shoeprintr)`
5. Run the following sample code, modifying the input files and parameters as necessary:

```
    ## Load in the input and the reference images
    data(input_example)
    data(reference_example)
    
    ## Transform all the points to have (0,0) at lower left corner.
    print_in <- leftcorner_cent(input_example)
    print_ref <- leftcorner_cent(reference_example)

    ## Perform Print Match
    print_stats <- match_print(print_in, print_ref,
			       ncross_in_bins = 30, xbins_in = 20, ncross_in_bin_size = 1,
			       ncross_ref_bins = NULL, xbins_ref = 30, ncross_ref_bin_size = NULL,
			       eps = .75, seed = 1, num_cores = parallel::detectCores(), 
			       plot = TRUE, verbose = TRUE)
    print_stats
```
6. (Optional) We provide the `pmc` binary with the R package. If you would like to update this binary to use the latest development version, follow the instructions at the bottom of this README file.

Deliverables
--------
1. The cleaned up, commented, documented code as an R package. The R package will enable easy sourcing of the utility functions and seamless documentation
2. A 10-15 slide rmarkdown presentation walking through 1) the changes, 2) the speed improvement, 3) live demo
3. A README for the package which walks through its installation and use
4. A spreadsheet of profiling results showing the speedup with various parameters compared to the original script

PMC Instructions
--------
1.	Clone the pmc library
	+ 	git clone https://github.com/ryanrossi/pmc.git
2.	Go into pmc directory
	+	cd pmc
3.	Compile the library
	+	make
4.	Replace the binary corresponding to your platform in inst/bin/<platform> with the newly compiled pmc binary
