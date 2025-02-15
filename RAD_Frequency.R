## Load packages
list.of.packages <- c("JuliaCall", "reticulate", "optparse", "ggplot2",
                      "reshape2", "optparse", "dada2", "ShortRead")
lapply(list.of.packages,library,character.only = TRUE)

##Specifying Illumina vs. PacBio files, and what the sample name is.
option_list <- list(make_option(c("-s", "--script_path"), type="character", default=NULL, help="Directory where scripts are located.", 
                                metavar="character"),
                    make_option(c("-c", "--cpus"), type="character", default=FALSE, help="task.cpus", metavar="character"),
                    make_option(c("-d", "--directory"), type="character", default=NULL, help="Specify working directory", metavar="character"),
                    make_option(c("-m", "--metadata"), type="character", default=NULL, help="Specify metadata", metavar="character"),
                    make_option(c("-a", "--pacbio_sample"), type="character", default=NULL, help="Specify PacBio file to run."), 
                    make_option(c("-n", "--sample_name"), type="character", default=NULL, help="Specify sample name."), 
                    make_option(c("-p", "--pacbio"), type="character", default=FALSE, help="Specify if these files are only PacBio.", 
                                metavar="character", action="store_true"),
                    make_option(c("-i", "--illumina"), type="character", default=FALSE, help="Specify if these files are only Illumina.", 
                                metavar="character", action="store_true"));
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

path <- opt$directory
script.dir <- opt$script_path

#####

## To run this script manually in R, uncomment the following lines. You do not need to change the preceding lines of path and script.dir,
## but remember to recomment the lines if you want to run the script automatically in the pipeline.
## path refers to the folder your metadata.csv and sequencing files (.fastq) are. (This is the -directory option).
## script.dir refers to the folder where all the script files are located. This should point to where you saved the cloned GitHub.

#path <- "/Users/administrator/Desktop/Illumina_PacBio"
#script.dir <- "/Users/administrator/Desktop/Illumina_PacBio"

## This script can also be run from the command line.
## Usage: rscript \path\to\og_files_to_all_reads.R -s [script_path] -d [directory]

#####

PacBio_fns <- c(as.character(opt$pacbio_sample))

## Identify primers to go from ATG to stop in tprK
tprKF <- "GGAAAGAAAAGAACCATACATCC"
tprKR <- "CGCAGTTCCGGATTCTGA"
rc <- dada2:::rc
#noprimer_filenames <- paste(substr(basename(PacBio_fns),1,nchar(basename(PacBio_fns))-5),"noprimers.fastq",sep ='')
noprimer_filenames <- paste(gsub("fastq.gz$", "noprimers.fastq", PacBio_fns))
nop <- file.path(noprimer_filenames)

print(paste("length pacbio_sample", length(opt$pacbiosample)))
print(paste("pacbio_sample", opt$pacbiosample))
print(paste("PacBIOFNS", PacBio_fns[0]))
print(paste("length PacBIOFNS", length(PacBio_fns)))
print(paste("No Primer Filenames",noprimer_filenames))
print(paste("Nop",nop))
print(paste("PacBio Fns",PacBio_fns))

print("Ending")

#weffsdffsfsdffsdfsdfdfsfsdf




if(opt$illumina == FALSE) {
  ## Points to Julia install in docker "quay.io/greninger-lab/tprk"
  julia <- julia_setup(JULIA_HOME = "/usr/local/julia/bin")
  #julia <- julia_setup(JULIA_HOME = "/Applications/Julia-0.6.app/Contents/Resources/julia/bin/")
  ## Remove primers
  for (count in c(1:length(nop))) {
    if(file.exists(nop[count])) {
      print(paste(noprimer_filenames[count], " already exists. Skipping removing primers step...", sep=""))
    } else {
      print("Removing primers from PacBio...")
      print(opt$metadata)
      print(PacBio_fns)
      nop <- tempfile(fileext=".fastq.gz")

      PacBio_fns2 <- file.path('/Users/administrator/Desktop/Illumina_PacBio/',PacBio_fns)

      print(paste("PacBio_fns2", PacBio_fns2))

      prim <- removePrimers(PacBio_fns2, nop, primer.fwd=tprKF, primer.rev=rc(tprKR), orient=TRUE, verbose=TRUE)
      #prim <- removePrimers(PacBio_fns, nop, primer.fwd=tprKF, primer.rev=rc(tprKR), orient=TRUE, verbose=TRUE)
    }
  }

  print("Filtering PacBio reads...")
  ## Setting up file names to filter.
  
  #filter_filenames <- paste(substr(basename(PacBio_fns),1,nchar(basename(PacBio_fns))-5),"noprimers.filtered.fastq",sep ='')
  filter_filenames <- paste(gsub("fastq.gz$", "noprimers.filtered.fastq", PacBio_fns))
  #filterEE1_filenames <- paste(substr(basename(PacBio_fns),1,nchar(basename(PacBio_fns))-5),"noprimers.filtered.EE1.fastq",sep ='')
  filterEE1_filenames <- paste(gsub("fastq.gz$", "noprimers.filtered.fastq", PacBio_fns))
  
  #filt <- file.path(filter_filenames)
  filt <- file.path('/Users/administrator/Desktop/Illumina_PacBio', filter_filenames)
  #filtEE1 <- file.path(filterEE1_filenames)
  filtEE1 <- file.path('/Users/administrator/Desktop/Illumina_PacBio', filterEE1_filenames)
  
  ## Filter reads for tprK length and do not worry about expected errors.
  for (count in c(1:length(filt))) {
    if (file.exists(filt[count])) {
      print(paste(filter_filenames[count]," already exists. Skipping filtering step..."), sep="")
    } else {
      print(paste("Filtering ",nop[count],"...",sep=""))
      track <- fastqFilter(nop[count], filt[count], minLen=1400,maxLen=1800,
                           maxN=0,
                           compress=FALSE, multithread=TRUE)
    }
  }
  
  ##Consider: Filter reads for tprK length and allow only 1 expected error for the entire read.
  # for (count in c(1:length(filtEE1))) {
  #   track <- fastqFilter(nop[count], filtEE1[count], minLen=1400,maxLen=1800,
  #                        maxN=0, maxEE=1,
  #                        compress=FALSE, multithread=TRUE)
  # }
  
  RAD_filenames <- paste(gsub("fastq.gz$", "noprimers.filtered.RAD.fasta", PacBio_fns))
  #RAD_filenames <- paste(substr(basename(PacBio_fns),1,nchar(basename(PacBio_fns))-5),"noprimers.filtered.RAD.fasta",sep ='')
  
  print(paste("RAD_filenames", RAD_filenames))
  RAD_files <- file.path('/Users/administrator/Desktop/Illumina_PacBio/',RAD_filenames)
  
     
  ## Build RAD files for each PacBio sample. This step takes forever!!!
  for (count in c(1:length(filt))) {
    to_rad_name <- paste(RAD_filenames[count])
    # Skips RAD step if files already exist, because it takes forever.
    if(file.exists(to_rad_name)) {
      print(paste(to_rad_name, " already exists. Skipping RAD step...", sep=""))
    } else{
      # Only want to set up Julia once, takes forever
      if (count == 1) {
        print("Setting up Julia...")
        print("Constructing RAD files...")
        julia_command("Pkg.init(); Pkg.update(); Pkg.clone(\"https://github.com/MurrellGroup/NextGenSeqUtils.jl\"); using NextGenSeqUtils")
        julia_command("Pkg.clone(\"https://github.com/MurrellGroup/DPMeansClustering.jl.git\")")
        julia_command("Pkg.clone(\"https://github.com/MurrellGroup/RobustAmpliconDenoising.jl.git\"); using RobustAmpliconDenoising")
      }
      # julia_command("using Pkg")
      # julia_command("Pkg.build(\"SpecialFunctions\")")
      # julia_command("Pkg.add(PackageSpec(name=\"NextGenSeqUtils\", rev= \"1.0\", url = \"https://github.com/MurrellGroup/NextGenSeqUtils.jl.git\"))")
      # julia_command("Pkg.add(PackageSpec(name=\"DPMeansClustering\", rev=\"1.0\", url = \"https://github.com/MurrellGroup/DPMeansClustering.jl.git\"))")
      # julia_command("Pkg.add(PackageSpec(name=\"RobustAmpliconDenoising\", rev=\"1.0\", url = \"https://github.com/MurrellGroup/RobustAmpliconDenoising.jl.git\"))")

      # julia_command("using RobustAmpliconDenoising")


      julia_readfastq <- paste("seqs, QVs, seq_names = read_fastq(\"",filt[count],'")',sep="")
      julia_command(julia_readfastq)
      julia_command("templates,template_sizes,template_indices = denoise(seqs)")
      julia_writefasta <- paste("write_fasta(\"",RAD_files[count],'",templates,names = ["seqs$(j)_$(template_sizes[j])" for j in 1:length(template_sizes)])',sep="")
      julia_command(julia_writefasta)
    }
  }
  
  ## RAD denoised files are written.  Let's get some frequencies of different variable regions
  RAD_files_nolines <- paste(substr(RAD_files,1,nchar(RAD_files)-5),"nolines.fasta",sep ='')
  RAD_files_fix <- paste(substr(RAD_files,1,nchar(RAD_files)-5),"nolines.fix.fasta",sep ='')

  # Fixes up the fastas so they wrap and don't have awkward new lines.
  # TODO: fix this section so it works. For some reason the pipeline currently runs without it? But probably should fix this anyway.
  awk_command <- paste("awk '/^>/ {printf(\"\\n%s\\n\",$0);next; } { printf(\"%s\",$0);}  END {printf(\"\\n\");}' < ",RAD_files," > ",RAD_files_nolines," ;")
  fix_firstline <- paste("tail -n+2 ",RAD_files_nolines," > ",RAD_files_fix)
  for (count in c(1:length(awk_command))) {
    system(awk_command[count])
    system(fix_firstline[count])
  
  }
} else {
  print("Illumina option specified. Skipping making PacBio frequency files...")
}
