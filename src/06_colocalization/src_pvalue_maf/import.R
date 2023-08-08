#/opt/apps/R/3.6.2/bin/R

# import.R
# Caleb Matthew Radens
# cradens@mail.med.upenn.edu
# 2016_1_19

# A series of functions for importing eQTL and GWAS summary statistics

system("echo ===========================")
R_ver <- substr(version$version.string,1,15) # Get R version
system("echo inside import")
system(paste("echo",R_ver))
system("echo ===========================")

lp <- "~/R/x86_64-conda_cos6-linux-gnu-library/3.6"
.libPaths(lp)

require(hash,lib=lp)
require(data.table,lib=lp)

system("echo wrapper package and script dependencies loaded and checked",wait=FALSE)

read_eQTL <- function(File_path,
                      Columns,
                      Skip=1,
                      Sep="\t",
                      Gene,
                      Chr_col,
                      Rsid_col,
                      Pos_col,
                      N_col,
                      N_eQTL,
                      PV_col,
                      Beta_col,
                      Varbeta_col,
                      Var_is_SE = TRUE){
  # This function imports one eQTL file and returns
  #		a table with summary stats for coloc.abf()
  #
  # Arguments:
  #		file_path: points to a eQTL summary statistics file.
  #			typeof: character
  #			format: "/.../TissueInfo_genomewide_significance.txt"
  #			 Note: TissueInfo may have multiple '_' in it.
  #			 Ex: Adipose_subcutaneous, or Adipose_visceral
  #		Columns: Integer number of columns in the file
  #		Skip: If there is meta-data or a header, Skip them.
  #		Chr_col: Column index with chr info
  #		Rsid_col: Column index with Rsid info
  #		Pos_col: Column index with hg19 position info
  #		N_col: eQTL sample size column index
  #   N_eQTL: eQTL sample size for all SNPs
  #		PV_col: P-value column index
  #		Beta_col: beta column index
  #		Varbeta: variance of beta column index
  # NOTE:
  #		File_path, Columns, Chr, Rsid, Pos, and (N_eQTL or N_col) are required.
  #		Either PV *or* Beta and Varbeta required. NOT ALL 3!!!
  #   If a row has a NA in it, it will be omitted!
  #
  # Returns: table with columns (either PV or beta and varbeta):
  # gene     | chr       | chr_pos  |rsid      | position  | N_eQTL  | PV_eQTL | beta_eQTL | varbeta_eQTL
  # character| character | character|character | integer   | integer | double  | double    | double

  # Check for filepath validity
  if(missing(File_path)){stop("Please include file_path")}
  if(!(file.exists(File_path))){stop(paste("File not found:\n",File_path))}

  # Checking that Columns, Chr, Rsid, and Position columns were specified.
  if(missing(Columns)
     || missing(Gene)
     || missing(Chr_col)
     || missing(Rsid_col)
     || missing(Pos_col)){
    stop("Must include Columns, Chr_col, Rsid_col, and Pos_cols.")
  }
  # Check that N_col or N were specified
  if(missing(N_col) && missing(N_eQTL)){
    stop("Must include N_col or N_eQTL.")
  }

  # Check that Columns is an integer
  if(!(is.numeric(Columns))){stop("Columns needs to be an integer")}
  if(Columns%%1!=0){stop("Columns needs to be an integer")}

  # Check that Skip is an integer >=0
  if(!(is.numeric(Skip))){stop("Skip needs to be an integer")}
  if(Skip%%1!=0){stop("Skip needs to be an integer")}
  if(Skip<0){stop("Skip needs to be >= 0")}

  # Check that Sep is a character
  if(!(is.character(Sep))){
    stop("Skip needs to be a *character* specifying how data are separated.")}

  # Check that Gene is a non-empty character
  if(!(is.character(Gene)) || length(Gene) == 0){
    stop(paste("Gene needs to be a non-zero-length character, not: '",
               as.character(Gene),"'",sep=""))
  }

  # Check that Rsid_col is an integer
  if(!(is.numeric(Rsid_col))){stop("Rsid_col needs to be an integer")}
  if(Rsid_col%%1!=0){stop("Rsid_col needs to be an integer")}

  # Check that Pos_col is an integer
  if(!(is.numeric(Pos_col))){stop("Pos_col needs to be an integer")}
  if(Pos_col%%1!=0){stop("Pos_col needs to be an integer")}

  # Initializing colClasses:
  #  NULL means the column won't be imported.
  colClasses <- rep("NULL", Columns)
  # Initializing what the column headers will be
  table_names <- rep("",Columns)

  # If PV argument supplied:
  if(!(missing(PV_col))){
    # Make sure beta and varbeta were *not* also supplied
    if(!(missing(Beta_col)) || !(missing(Varbeta_col))){
      stop("Please only include PV_col *OR* Beta_col and Varbeta_col.")
    }
    else{
  		# Initializing PV as character because fread() doesn't like super
  		#  small numbers
			table_names[PV_col] <- "PV_eQTL"
			colClasses[PV_col] <- "double"
    }
  }
  # If beta and varbeta supplied:
  else if(!(missing(Beta_col)) && !(missing(Varbeta_col))){
  	table_names[Beta_col] <- "beta_eQTL"
  	colClasses[Beta_col] <- "double"

  	table_names[Varbeta_col] <- "varbeta_eQTL"
  	colClasses[Varbeta_col] <- "double"
  }
  # Else, use didn't supply the appropriate arguments:
  else{stop("Please supply PV_col *or* Beta_col and Varbeta_col")}

  # If N_col supplied:
  if(!missing(N_col)){
    table_names[N_col] <- "N_eQTL"
    colClasses[N_col] <- "integer"
  }
  # Else make sure N_eQTL was supplied:
  else if (missing(N_eQTL)){
    stop("Must specify N_col or N_eQTL.")
  }

  table_names[Chr_col] <- "chr"
  colClasses[Chr_col] <- "character"

  table_names[Pos_col] <- "position"
  colClasses[Pos_col] <- "integer"

  table_names[Rsid_col]<- "rsid"
  colClasses[Rsid_col] <- "character"

  # Remove all the empty columns from table_names
  table_names <- table_names[!table_names %in% ""]

  # na.omit() removes rows from a table if any of the columns have an 'NA'
  table <- na.omit(fread(         # Omit a row that has NA in it!
      input=File_path,
      skip=Skip,                  # Skip lines when importing
      sep=Sep,                    # Indicate what elements are separated by
      stringsAsFactors = FALSE,   # Strings should be strings, not factors
      colClasses = colClasses,    # Define column classes (see above)
      data.table = FALSE          # Output a data.frame, not table
      ))

  # Name the column headers
  names(table) <- table_names

  # P-values imported:
  if(!(missing(PV_col))){
    # Convert pvalues to doubles; any pvalue below R's minimum recognizable number
    #  is converted to 1e-300:
    table$PV_eQTL <- as.double(table$PV_eQTL)

    indeces <- which(table$PV_eQTL == 0)

    table$PV_eQTL[indeces] <- 1e-300
  }

  if (!(missing(Varbeta_col))){
    # If the variance column is actually SE,
    if (Var_is_SE){
      # This column is standard error, so square it to make it into variance
      table$varbeta_eQTL <- (table$varbeta_eQTL)^2
    }

  }

  # Add chr_pos column to table
  table$chr_pos <- paste(table$chr, as.character(table$position), sep=":")

  # Add gene name column
  table$gene <- Gene

  if (!missing(N_eQTL)){
    table$N_eQTL = N_eQTL
  }

  table <- table[table$PV_eQTL<0.01, ]

  return(table)
}

get_gene_names <- function(directory, Pattern=".txt"){
  # Given a directory containing eQTL .txt files,
  #   extract the gene names and their chromosmes. Pattern is
  #   defaulted to identify CLI snptest files.
  #
  # Arguments:
  #   directory:
  #		  typeof: character
  #   Pattern: grep pattern to identify relevant files
  #     typeof: character
  #
  # Return: dataframe with columns:
  #     gene  | filename

  # Extract all file names with matching grep pattern from directory
  base_file_name <- list.files(path=directory, pattern=Pattern)

  # Split file name at "."
  split <- strsplit(x=base_file_name, split=".txt")

  # Unlist the lists into matrices (the file has 2 segments)
  # Each row is a split file name
  split <- matrix(unlist(split),ncol=1,byrow=TRUE)

  # Build dataframe
  genes <- data.frame(
              gene=split[,1], 			#GeneName
              filename=base_file_name,
              filepath = paste(directory, base_file_name,sep=""),
              stringsAsFactors = FALSE	# Strings should be strings
              )

  return(genes)
}

import_genes <- function(Genes, Directory, Pattern=".txt"){
  # Get gene names with get_gene_names() and eQTL.txt files from
  #    directory using read_eQTL.
  #
  #	Argument:
  #   Genes: list of gene names to import
  #     typeof: character vector
  #     format: c("GENE1", "GENE2", ...)
  #   Directory:
  #		  typeof: character
  # 	  format: "/.../gene_name_eQTL.txt"
  #
  # Returns: a dictionary with gene names pointing to imported
  # 	eQTL summary stat tables.

  # Initialize gene dictionary
  gene_dict <- hash()

  # Import all gene names, chromosomes, and filepaths
  all_genes <- get_gene_names(Directory, Pattern)

  # Convert list of Genes to a data frame
  Genes <- data.frame(gene=Genes)

  # Merge dataframes, keeping only the genes to be imported
  #		all=FALSE means:
  #			-If one of the user-defined argument Genes isn't
  #       found in the Directory, it will be excluded from
  #       this dataframe.
  #			-Likewise, if a gene identified in the Directory
  #       isn't specified in argument Genes, it will be excluded.
  genes_to_import <- merge(Genes, all_genes, by="gene", all=FALSE)

  # Cycle through gene list
  for (row in seq(from=1, to=length(genes_to_import$gene))){
    g <- as.character(genes_to_import$gene[row])    #edit by zhec
    # Building filepath where gene eQTL file is located
    file_path <- paste(Directory,
                      genes_to_import$filename[row],
                      sep="")

    # Read in the imported eQTL summary stat table and
	  # 		add it to the gene dictionary


    gene_dict[g]<-read_eQTL(File_path=file_path,
                            Columns=5,
                            Skip=1,
                            Sep="\t",
                            Gene=g,
                            Chr_col=1,
                            Rsid_col=5,
                            Pos_col=2,
                            N_col=4 ,
                            PV_col=3,
                            Var_is_SE = FALSE)  #edit by zhec
  }
  index <- duplicated(gene_dict[,2]) #edit by HYM
  gene_dict <- gene_dict[!index,] #edit by HYM
  return(gene_dict)

}

read_GWAS <- function(File_path,
                      Columns,
                      Skip=1,
                      Sep="\t",
                      Chr_col,
                      Chr_pos_col,
                      Rsid_col,
                      Pos_col,
                      MAF_col,
                      N_col,
                      N_GWAS,
                      PV_col,
                      Beta_col,
                      Varbeta_col,
                      Var_is_SE = TRUE){
  # This function imports one GLGC .txt file and
  # 	returns a table with summary stats for coloc.abf()
  #
  # Argument:
  #   file_path:
  #		  typeof: character
  #		  format: "base_directory/.../_blah_GWAS_blah.txt"
  #
  # Arguments:
  #		File_path: points to a GWAS summary statistics file.
  #			typeof: character
  #			format: "/.../GLGC_blah.txt"
  #		Columns: Integer number of columns in the file
  #		Skip: If there is meta-data or a header, Skip them.
  #		Chr_col: Column index with chr info
  #   Chr_pos_col: Column index with Chr:hg19_pos###
  #		Rsid_col: Column index with Rsid info
  #		Pos_col: Column index with hg19 position info
  #		MAF_col: Column index withe the MAF info
  #		N_col: GWAS sample size column index
  #   N_GWAS: GWAS sample size for all SNPs
  #		PV_col: P-value column index
  #		Beta_col: beta column index
  #		Varbeta: variance of beta column index
  # NOTE:
  #		N_col or N_GWAS required!
  #   If a row has a NA in it, it will be omitted!
  #
  # Returns: table with columns:
  # chr       | chr_pos  | rsid      | position  | MAF    | N_GWAS  | PV_GWAS | beta_GWAS | varbeta_GWAS
  # character | character| character | integer   | double | integer | double  | double    | double
  #
  # Notes:
  #   - pvalues that are below R's minimum number are converted to 1e-300
  #   -  NOTE: values equal to 0 are also converted to 1e-300
  #   - If 'MAF' is >0.5, transform MAF into 1-MAF. If MAF is 1, remove the SNP.
  #   - If 'MAF' is <0.001, remove SNP because Giambartolomei et. al. did so...

  # Check for filepath validity
  if(missing(File_path)){stop("Please include file_path")}
  if(!(file.exists(File_path))){stop(paste("File not found:\n",File_path))}

  # Checking that all required arguments are specified.
  if(missing(Columns)
     || missing(Chr_col)
     || missing(Chr_pos_col)
     || missing(Rsid_col)
     || missing(Pos_col)
     || missing(MAF_col)
     || missing(PV_col)){
    stop("Must include Columns, Chr_col, Chr_pos_col, Rsid_col, Pos_cols., MAF_col, PV_col, and N_GWAS or N_Col")
  }
  # Check that N_col or N were specified
  if(missing(N_col) && missing(N_GWAS)){
    stop("Must include N_col or N_GWAS.")
  }

  # Check that Columns is an integer
  if(!(is.numeric(Columns))){stop("Columns needs to be an integer")}
  if(Columns%%1!=0){stop("Columns needs to be an integer")}

  # Check that Chr_col is an integer
  if(!(is.numeric(Chr_col))){stop("Chr_col needs to be an integer")}
  if(Chr_col%%1!=0){stop("Chr_col needs to be an integer")}

  # Check that Chr_pos_col is an integer
  if(!(is.numeric(Chr_pos_col))){stop("Chr_pos_col needs to be an integer")}
  if(Chr_pos_col%%1!=0){stop("Chr_pos_col needs to be an integer")}

  # Check that Skip is an integer >=0
  if(!(is.numeric(Skip))){stop("Skip needs to be an integer")}
  if(Skip%%1!=0){stop("Skip needs to be an integer")}
  if(Skip<0){stop("Skip needs to be >= 0")}

  # Check that Sep is a character
  if(!(is.character(Sep))){
    stop("Skip needs to be a *character* specifying how data are separated.")}

  # Check that Rsid_col is an integer
  if(!(is.numeric(Rsid_col))){stop("Rsid_col needs to be an integer")}
  if(Rsid_col%%1!=0){stop("Rsid_col needs to be an integer")}

  # Check that Pos_col is an integer
  if(!(is.numeric(Pos_col))){stop("Pos_col needs to be an integer")}
  if(Pos_col%%1!=0){stop("Pos_col needs to be an integer")}

  # Check that MAF_col is an integer
  if(!(is.numeric(MAF_col))){stop("MAF_col needs to be a double")}
  if(MAF_col%%1!=0){stop("MAF_col needs to be an integer")}

  # Initializing colClasses:
  #  NULL means the column won't be imported.
  colClasses <- rep("NULL", Columns)
  # Initializing what the column headers will be
  table_names <- rep("",Columns)

  if(!(missing(N_col))){
    if(!(missing(N_GWAS))){ # Make sure only N_col or N specified
      stop("Please only include N_col *or* N_GWAS.")
    }
    table_names[N_col] <- "N_GWAS"
    colClasses[N_col] <- "integer"
  }

  # If beta and varbeta supplied:
  if(!(missing(Beta_col)) && !(missing(Varbeta_col))){
  	table_names[Beta_col] <- "beta_GWAS"
  	colClasses[Beta_col] <- "double"

  	table_names[Varbeta_col] <- "varbeta_GWAS"
  	colClasses[Varbeta_col] <- "double"
  }

  table_names[Chr_col] <- "chr"
  colClasses[Chr_col] <- "character"

  table_names[Chr_pos_col] <- "chr_pos"
  colClasses[Chr_pos_col] <- "character"

  table_names[Pos_col] <- "position"
  colClasses[Pos_col] <- "integer"

  table_names[Rsid_col]<- "rsid"
  colClasses[Rsid_col] <- "character"

  table_names[MAF_col]<- "MAF"
  colClasses[MAF_col] <- "double"

  # PV is initialized as a character because fread() doesn't like
  #   super small numbers.
  table_names[PV_col] <- "PV_GWAS"
  colClasses[PV_col] <- "character"

  # Remove all the empty columns from the table names
  table_names <- table_names[!table_names %in% ""]

  # na.omit() removes rows from a table if any of the columns have an 'NA'
  table <- na.omit(fread(         # Omits rows that have any NAs in them!
      input=File_path,
      skip=Skip,                  # Skip lines when importing
      sep=Sep,                    # Indicate what the elements are separated by
      stringsAsFactors = FALSE,   # Strings should be strings, not factors
      colClasses = colClasses,    # Define column classes (NULL is skipped)
      data.table = FALSE          # Output a data.frame, not table
      ))

  # Name the column headers
  names(table) <- table_names

  # P-values imported:
  if(!(missing(PV_col))){
  	# Convert pvalues to doubles; any pvalue below R's minimum recognizable number
  	#  is converted to 1e-300:
  	table$PV_GWAS <- as.double(table$PV_GWAS)

  	indeces <- which(table$PV_GWAS == 0)

  	table$PV_GWAS[indeces] <- 1e-300
  }
  if (!(missing(Varbeta_col))){
    # If the variance column is actually SE,
    if (Var_is_SE){
      # This column is standard error, so square it to make it into variance
      table$varbeta_GWAS <- as.numeric((table$varbeta_GWAS))^2
    }

  }

  # MAF > 0.5? Turn it into 1-MAF
  table$MAF <-as.numeric(table$MAF)
  table$MAF[which(table$MAF>0.5)] <- 1 - table$MAF[which(table$MAF>0.5)]

  # Remove SNPs with MAF < 0.001. Giambartolomei et. al. do so...
  # Also, MAF of 0 causes errors in coloc.abf()
  table <- table[table$MAF>0.01, ] ### MAF>0.01 edit by hchen 20220119

  if (!(missing(N_GWAS))){
    table$N_GWAS = N_GWAS
  }

  return(table)
}


get_GWAS_traits <- function(directory, Pattern="GLGC"){
  # Given a directory of GWAS files, extract the traits
  #
  # Arguments:
  #   directory
  #		  typeof: character
  #   Pattern (optional)
  #     typeof: character
  #
  # 	Expected format of GLGC files: "GLGC_blah.txt"
  #
  # Return dataframe:
  #     trait		| 	filename
  #		character	|	character

  # Grep extract all file names with Pattern
  base_file_name <- list.files(path=directory, pattern=Pattern)

  # Split file names at "_"
  split_underscore <- strsplit(x=base_file_name, split="_")
  # Unlist list of split file names (three segments) into a matrix
  # Each row is a filename with columns:
  #   [1] == "GLGC"
  #   [2] == "trait"
  #   [3] == "results.txt"
  split_underscore <- matrix(unlist(split_underscore),ncol=3,byrow=TRUE)
  # Build traits dataframe
  traits<- data.frame(
    trait=split_underscore[,2],		# The trait from the file name
    filename=base_file_name,
    filepath=paste(directory,base_file_name,sep=""),
    stringsAsFactors = FALSE	    # Strings should be strings, not factors
    )

  return(traits)
}

import_GWAS_tables <- function(traits, GWAS_dir){
  # Given a GWAS trait or list of traits (character vector) and directory,
  # 	cross-check that traits matches those in the directory, then call
  #		read_GWAS to import the appropriate data. If the cross-check
  #		fails for a trait, the trait won't be imported, but the function
  #		will still import all the other trait-associated files.
  #
  # Arguments:
  #		traits: list of GWAS traits to import GWAS stats from
  #			typeof: character (vector)
  #			length: at least 1
  #			format: c("TRAIT1", "TRAIT2", etc)
  #		GWAS_dir: directory containing GWAS .txt file(s)
  #			typeof: character
  #
  # Returns:
  #		a dictionary with trait(s) pointing to imported GWAS table(s)

  # Initialize the trait dictionary
  trait_dict <- hash()

  # Import all traits from directory into a data frame
  all_traits <- get_GWAS_traits(GWAS_dir)

  # Convert list of traits to a data frame
  traits <- data.frame(trait=traits)

  # Merge dataframes, keeping only the traits to be imported
  #		all=FALSE means:
  #			-If one of the user-defined argument traits isn't
  #       found in the Directory, it will be excluded from
  #       this dataframe.
  #			-Likewise, if a trait identified in the Directory
  #       isn't specified in argument traits, it will be excluded.
  traits_to_import <- merge(traits, all_traits, by="trait", all=FALSE)

  for (trait in traits_to_import$trait) {
    # Building filepath where trait GWAS file is located
    file_path <- paste(GWAS_dir,
                       traits_to_import$filename[which(traits_to_import$trait==trait)],
                       sep="")

    # Assign the imported GWAS table to the list of traits
    trait_dict[trait]<-read_GWAS(File_path=file_path,    #edit by zhec
                                  Columns=9,
                                  Skip=1,
                                  Sep="\t",
                                  Chr_col=1,
                                  Chr_pos_col=3,
                                  Rsid_col=4,
                                  Pos_col=2,
                                  MAF_col=9,
                                  N_col=7,
                                  PV_col=8,
                                  Beta_col=5,
                                  Varbeta_col=6,
                                  Var_is_SE = TRUE)      #edit by zhec
  }
  return(trait_dict)
}

not_too_tiny <- function(numbers, new_tiny = 1e-300){
  #   R converts numbers that are too small to 0. To prevent this, use mpfr to identify
  #     these very small numbers, then convert them to a new arbitrarily small number.
  #
  #   By default, 'too small' is less than 1e-300 and 1e-300 is the new smallest number.
  #
  #   Arguments:
  #     characters: a double vector of numbers to be checked for smallness
  #       typeof: double
  #       format: #.# [# may be any representation of a number that R recognizes]
  #
  #     new_tiny: the cutoff minimum number
  #       typeof: double
  #
  #   Returns:
  #     A double vector with all numbers less than new_tiny and greater than 0 converted
  #       to new_tiny

  # Identify the numbers that are too small (excpet for 0 or negative numbers)
  too_tiny <- which(numbers<new_tiny)

  # Overwrite the too small numbers with the cutoff value
  numbers[too_tiny]<-new_tiny

  return(numbers)
}


not_too_big <- function(numbers, new_big = 1e+300){
  #   R converts numbers that are too big to Inf. To prevent this,  identify
  #     these very big numbers, then convert them to a new max number.
  #
  #   By default, 'too big' is greater than 1e+300 and 1e+300 is the new biggest number.
  #
  #   Arguments:
  #     characters: a double vector of numbers to be checked for bigness
  #       typeof: double
  #       format: #.# [# may be any representation of a number that R recognizes]
  #
  #     new_tiny: the cutoff maximum number
  #       typeof: double
  #
  #   Returns:
  #     A double vector with all numbers greater than new_big converted
  #       to new_big

  # Identify the numbers that are too small (excpet for 0 or negative numbers)
  too_big <- which(numbers>new_big)

  # Overwrite the too small numbers with the cutoff value
  numbers[too_big]<-new_big

  return(numbers)
}

bedify <- function(Table, Value, Chr="chr", Start="position", End="position",Loci="position", remove_0s=TRUE) {
  # Given a table, return a bedfile style table. User needs to specify which column to get
  #   values from, but the start and end positions as well as the loci of interest are defaulted
  #   to the position column.
  bedtable <- data.frame(c(
      Table[Chr],
      Table[Start],
      Table[End],
      Table[Loci],
      Table[Value]
      ),
      stringsAsFactors = FALSE)

  colnames(bedtable) <- c("chr", "start", "end", "loci", "value")

  # If specified, remove rows with value of 0
  if(remove_0s){
    bedtable<-bedtable[bedtable$value>0,]
  }

  return(bedtable)
}

get_chrom_lengths <- function(Chromosomes=as.character(seq(from=1,to=22))){
  # Extracts hg19 chromsomes lengths from ucsc.
  # Defualts to extracting chromsomes 1->22 (excludes X and Y)

  # Build table with chromsome start and end info
  # https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
  # backup table based on hg19 as of 11/13/2015:
  #   V1=c("chr1", "chr2",
  #        "chr3", "chr4",
  #        "chr5", "chr6",
  #        "chr7", "chr8",
  #        "chr9", "chr10",
  #        "chr11", "chr12",
  #        "chr13", "chr14",
  #        "chr15", "chr16",
  #        "chr17", "chr18",
  #        "chr19", "chr20",
  #        "chr21", "chr22")
  #   V2=c(249250621, 243199373,
  #        198022430, 191154276,
  #        180915260, 171115067,
  #        159138663, 146364022,
  #        141213431, 135534747,
  #        135006516, 133851895,
  #        115169878, 107349540,
  #        102531392, 90354753,
  #        81195210, 78077248,
  #        59128983, 63025520,
  #        48129895, 51304566)
  #   genome<-data.frame(V1, V2, stringsAsFactors = FALSE)

  # Website with chromosome length info:
  #link <- '~/2021-10-31-cancer-GWAS/00data/hg19.chrom.sizes'

  # Import table
  table <- read.table("~/2021-10-31-cancer-GWAS/00data/hg19.chrom.sizes", stringsAsFactors = FALSE)

  # Initialize chromsomes and lengths vectors
  chromsomes<- c()
  lengths <- c()

  # for each chromosme in Chromosme, add "chr" before it
  for (num in Chromosomes) {
    chr <- paste("chr",num,sep="")

    # If chr# is in the table:
    if (chr %in% table[,1]){

      # Extract the chromosome length from the table
      lengths <- append(lengths, table[table[,1]==chr,2])
      chromsomes <- append(chromsomes, chr)

      # Else chr# isn't in the table, throw an error:
    } else{
      stop(paste(chr, " was not found in the ucsc hg19 chrom sizes table."))
    }
  }
  return(data.frame(chromsomes, lengths, stringsAsFactors = FALSE))
}

get_command_args <- function(args){
  # If a script is instantiated in PMACS as such:
  # >bsub Rscript my_fav_script.R --argument1_handle argument1 --argument2_handle argument2
  #
  # And then the script does:
  # args <- commandArgs(trailingOnly = TRUE)
  #
  # and passes args into this function,
  # return the arguments as a list.
  #
  # The arguments may be extracted from the list by:
  # list[[argument1_handle]]

  temp <- paste(unlist(args),collapse=' ')
  listoptions <- unlist(strsplit(temp,'--'))[-1]
  options_args <- sapply(listoptions,function(x){
    unlist(strsplit(x, ' '))[-1]
  })
  options_names <- sapply(listoptions,function(x){
    option <-  unlist(strsplit(x, ' '))[1]
  })
  names(options_args) <- unlist(options_names)
  return(options_args)
}


import_microarray_table <- function(){
  microArray_alnTable_file <- "~/2021-10-31-cancer-GWAS/aQTL_coloc/output/microArray/final/allTissue_hg19_microArray.txt"    #edit by zhec: change the path of microArray_alnTable
  microArray_alnTable <- na.omit(fread(input=microArray_alnTable_file,
                                       stringsAsFactors = FALSE,
                                       data.table = FALSE,
                                       colClasses = c("character",      #edit by zhec: change the colClasses
                                                      "character",
                                                      rep("integer",2))))
  names(microArray_alnTable) <- c("gene","chr","TSS", "TES")
  g <- microArray_alnTable$gene
  unique_gene_indeces <- order(g)[!duplicated(sort(g))]
  microArray_alnTable <- microArray_alnTable[unique_gene_indeces,]
  t <- microArray_alnTable$TES
  unique_TES_indeces <- order(t)[!duplicated(sort(t))]
  microArray_alnTable <- microArray_alnTable[unique_TES_indeces,]
  TES_dictionary <- hash()
  for (row in seq(from=1,to=length(microArray_alnTable$gene))){
    TES_dictionary[[microArray_alnTable$gene[row]]]<-list(chr=microArray_alnTable$chr[row],
                                                          TSS=microArray_alnTable$TSS[row],
                                                          TES=microArray_alnTable$TES[row])
  }
  return(TES_dictionary)
}

get_gene_TssTes <- function(Genes, only_TssTes = TRUE,
                            transcripts_file="~/2021-10-31-cancer-GWAS/aQTL_coloc/output/bed/hg19/allTissue_hg19.bed"){
  # This function imports a human genome transcripts file, and returns a dictionary of
  # gene TSSs and TESs given a vector of genes to look up.
  #
  # Arguments:
  #     transcripts_file: filepath to the transcript file. By default, it is from gencode V19
  #       and looks like:
  #   V1   | V2  | V3  | V4               | V5        | V6
  #   chr# | TSS | TES | common_gene_name | ensemblID | strand
  #
  #     Genes: c("vector", "of", "common", "genes", "like", "SORT1",...)

  #       data in the file be returned, too?
  #
  # Returns:
  #   gene dictionary
  #     COMMON_GENE_NAME -> list(TSS=c(vector of TSSs), TES=c(vector of TESs))
  #
  # Note: if a gene to lookup isn't in the table, the dictionary points to a string
  #   instead of a list

  # Get rid of potential duplicates
  Genes <- unique(Genes)

  col_classes <- c("character",
                   "integer",
                   "integer",
                   rep("character",3))
  table <- fread(transcripts_file,
                 stringsAsFactors = FALSE,
                 colClasses = col_classes)

  gene_dict <- hash()
  for (gene in Genes){
    indeces<-which(table$V4==gene)
    if(length(indeces)==0){
      gene_dict[[gene]] <- paste(gene,"not in transcript file")
    } else if (only_TssTes){
      gene_dict[[gene]] <- list(TSS=table$V2[indeces],
                                TES=table$V3[indeces])
    } else {
      gene_dict[[gene]] <- list(chr=unique(table$V1[indeces]),
                                TSS=table$V2[indeces],
                                TES=table$V3[indeces],
                                ENS=table$V5[indeces],
                                Strand=table$V6[indeces])
      if(length(unique(table$V1[indeces])) > 1){stop(paste("There is more than 1 chr for:",gene))}
    }
  }
  return(gene_dict)
}




get_sentinal_snp_regions <- function(GWAS_table,Range=1e6,Cutoff_PV=5e-8){
  # Given a GWAS table, return sentinal SNPs:
  #   Sentinal SNPs are those SNPs whose pvalues are the lowest within each region in the genome
  #   Each sentinal SNP is at least ~1MB away from other sentinal SNPs (thus definining the 'regions')
  #
  # Arguments:
  #   GWAS_table: table with a least the following columns (with indicated names):
  #     chr_pos | PV_GWAS
  #   Range: how far apart sentinal SNPs need to be
  #   Cutoff_PV: what is considered to be genome-wide significant (minimum allowed sentinal SNP PV)
  #
  # Returns:
  #   GWAS_table[at loci determined to be the senitnal SNPs,]
  #
  significant_PV_indeces <- which(GWAS_table$PV_GWAS<=Cutoff_PV)
  print(length(significant_PV_indeces))
  if (length(significant_PV_indeces) <=0){
    write.table("hello world",file=paste("sentinalSNP","_was_empty.txt",sep=""))
    stop(paste("Table was empty:"," sentinalSNP"))
  }
  #print(significant_PV_indeces)
  significant_PVs <- GWAS_table$PV_GWAS[significant_PV_indeces]
  chr_loci <- matrix(unlist(strsplit(GWAS_table$chr_pos[significant_PV_indeces],split = ":")),ncol=2,byrow=TRUE)
  df <- data.frame(pvalues = significant_PVs,
                   index = significant_PV_indeces,
                   chr =chr_loci[,1],
                   position = as.integer(chr_loci[,2]),
                   stringsAsFactors = FALSE)
  # Sort significant PV data frame by chrome, then by PV (maintaining trait GWAS_table indeces and loci)
  df <- df[order(df$chr,df$pvalues),]
  # Starting with the top (smallest) PV index in each chrome, add PV indeces to keepers
  #   if they are more than Range away on genome from the PVs already in keepers
  keepers <- data.frame()
  for (chrome in unique(as.character(df$chr))){
    df_temp <- df[which(df$chr==chrome),]
    keepers <- rbind(keepers,df_temp[1,])
    # Keep the loci of the SNP with the smallest PV in the current chromsome
    while (length(df_temp$position>0)){
      last_keeper <- keepers[length(keepers$position),]
      bottom_range <- last_keeper$position - Range
      top_range <- last_keeper$position + Range
      not_too_close <- c(which(df_temp$position<=bottom_range),which(df_temp$position>=top_range))
      df_temp <- df_temp[not_too_close,]
      if(length(df_temp$pvalues)==0){
        break
      }
      keepers<-rbind(keepers, df_temp[1,])
    }
  }
  sentinal_indeces <- keepers$index
  return(GWAS_table[sentinal_indeces,])
}

get_microArray_table <- function(filepath){
  microArray_alnTable <- na.omit(fread(input=filepath,
                                       stringsAsFactors = FALSE,
                                       data.table = FALSE,
                                       sep="\t",
                                       colClasses = c("character",
                                                      "character",
                                                      rep("integer",2))))
  names(microArray_alnTable) <- c("gene","chr","TSS", "TES")
  return(microArray_alnTable)
}

read_1000genome_table <- function(file_path){
  # Imports a 1000 genome table and returns a table.
  #   Only imports specified columns:
  col_classes <- c(rep("character",5),rep("NULL",9))

  # Filepath minus the .gz, for when the file is unzipped
  non_gz_file_path<-substr(file_path,start=1,stop=nchar(file_path)-3)

  # fread cannot hangle .gz files...
  # unzip it using linux's gunzip (wait until it is unzipped before proceeding)
  # -k to keep the zipped file
  system(paste("gunzip -c",file_path,">",non_gz_file_path),wait = TRUE)

  # na.omit() removes rows from a table if any of the columns have an 'NA'
  table <- fread(
    input=non_gz_file_path,
    sep=" ",                 # elements are separated by a space
    stringsAsFactors = FALSE,# Strings should be strings, not factors
    header = TRUE,
    # skip=1,
    colClasses = col_classes,  # Only import specified columns (see above)
    data.table=FALSE
  )

  # Delete unzipped file
  unlink(non_gz_file_path)


  return(table)
}
