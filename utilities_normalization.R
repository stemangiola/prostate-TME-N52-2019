#' Eliminate lowly tanscribed genes for normalization
#'
#' @param df A matrix
#' @param cpm_threshold A number
#' @param prop A number
#' @return A tibble filtered
rnaseq_norm.get_cpm = function(df,cpm_threshold = 0.5){

	cpm_threshold =
		cpm_threshold /
		(
			df %>%
				dplyr::group_by(sample) %>%
				dplyr::summarise(s = sum(value)) %>%
				dplyr::ungroup() %>%
				dplyr::summarise(m =median(s)) %>%
				dplyr::pull(m) /
				1e6
		)

	df %>%
	{ if("ct"%in%names(.)) dplyr::select(-ct)	else .}	%>%
		dplyr::select(gene, sample, value) %>%
		tidyr::spread(sample, value) %>%
		tidyr::drop_na() %>%
		dplyr::do(
			dplyr::bind_cols(
				gene = (.)$gene,
				tibble::as_tibble( edgeR::cpm(	(.) %>% dplyr::select(-gene) ) )
			)
		) %>%
		tidyr::gather(sample, cpm, -gene) %>%
		dplyr::mutate(cpm_threshold = cpm_threshold)

}


#' Eliminate lowly tanscribed genes for normalization
#'
#' @param df A matrix
#' @param cpm_threshold A number
#' @param prop A number
#' @return A tibble filtered
rnaseq_norm.get_low_expressed = function(df,cpm_threshold = 0.5, prop = 3/4){

	df %>%
		tidyr::spread(sample, value) %>%
		tidyr::drop_na() %>%
		
		# Filter low transcribed
		dplyr::filter(
			rowSums(
				edgeR::cpm((.) %>% dplyr::select(-gene)	) > 
					
					# Adjust threshold
					cpm_threshold /
					(
						df %>%
							dplyr::group_by(sample) %>%
							dplyr::summarise(s = sum(value)) %>%
							dplyr::ungroup() %>%
							dplyr::summarise(m =median(s)) %>%
							dplyr::pull(m) /
							1e6
					)
			) <
				
			# Filter based on how many samples are lowly transcribed
			floor(ncol((.) %>% dplyr::select(-gene)	) * prop)
		) %>%
		
		# Pull information
		dplyr::pull(gene) %>%
		as.character()
}

error_if_log_transformed = function(x){
	if(length(x$value)>0) if(max(x$value, na.rm=T)<50)
		stop("ARMET: The input was log transformed in: check_if_sd_zero_and_correct")
}

#' Calculate the norm factor with calcNormFactor from limma
#'
#' @param df A matrix
#' @param reference A reference matrix
#' @param cpm_threshold A number
#' @param prop A number
#' @return A list including the filtered data frame and the normalization factors
rnaseq_norm.calcNormFactor = function(
	df, 
	reference = NULL, 
	cpm_threshold = 0.5, 
	prop = 3/4
){

	error_if_log_transformed(df)

	# Get list of low transcribed genes
	gene_to_exclude =
		rnaseq_norm.get_low_expressed(
			df %>% dplyr::filter(sample!="reference") ,
			cpm_threshold = cpm_threshold,
			prop = prop
		)

	if(length(gene_to_exclude) == df %>% dplyr::distinct(gene) %>% nrow()) stop("ARMET: The gene expression matrix has been filtered completely for lowly expressed genes")

	writeLines(sprintf("ARMET: %s genes excluded for normalization", length(gene_to_exclude)))

	df.filt =
		df %>%
		dplyr::filter(!gene %in% gene_to_exclude) %>%
		droplevels()

	list(
		gene_to_exclude = gene_to_exclude,
		nf =
			tibble::tibble(
				sample = factor(levels(df.filt$sample)),
				nf = edgeR::calcNormFactors(
					df.filt %>%
						tidyr::spread(sample, value) %>%
						dplyr::select(-gene),
					refColumn=which(reference == factor(levels(df.filt$sample))),
					method="TMM"
				)
			) %>%
			dplyr::left_join(
				df.filt %>%
					dplyr::group_by(sample) %>%
					dplyr::summarise(tot_filt = sum(value, na.rm = T)) %>%
					dplyr::mutate(sample = as.factor(as.character(sample))),
				by = "sample"
			)
	)
}

#' Normalize a RNA seq data set using rnaseq_norm.calcNormFactor
#'
#' @param df A matrix
#' @param reference A reference matrix
#' @param cpm_threshold A number
#' @param prop A number
#' @return A list including the filtered data frame and the normalization factors
norm_RNAseq = function(
	input.df, 
	reference = NULL, 
	cpm_threshold = 0.5, 
	prop = 3/4, 
	genes_to_keep = c(), 
	method = "TMM",
	sample_column = "sample",
	gene_column = "gene",
	value_column = "value"
){

	# Reformat input data set
	df = input.df %>% 
		
		# Rename
		dplyr::select(!!sample_column, !!gene_column, !!value_column) %>%
		setNames(c("sample", "gene", "value")) %>%
	
		# Set samples and genes as factors
		mutate(
			sample = factor(sample), 
			gene = factor(gene)
		) 
	
	# Check if columns are available
	if(length(intersect( c("sample", "gene", "value"), colnames(df))) < 3 ) 
		stop("ARMET: input table is not correctly formatted as gene, sample")

	# Get norm factor object
	reference = 
		df %>% 
		group_by(sample) %>% 
		summarise(sum = sum(value)) %>% 
		mutate(med = median(sum)) %>% 
		mutate(diff = abs(sum-med)) %>%
		arrange(diff) %>%
		head(n=1) %>% 
		pull(sample) %>%
		as.character()
	
	nf_obj = rnaseq_norm.calcNormFactor(df, reference, cpm_threshold, prop)

	# Calculate normalization factors
	nf = nf_obj$nf %>%
		dplyr::left_join(
			df %>%
				dplyr::group_by(sample) %>%
				dplyr::summarise(tot = sum(value, na.rm = T)) %>%
				dplyr::mutate(sample = as.factor(as.character(sample))),
			by = "sample"
		) %>%
		dplyr::mutate(
			multiplier = 
				1 / 
				(tot_filt * nf) * 
				((.) %>% filter(sample == reference) %>% pull(tot)) 
		) %>%
		{
			# I have correct the strange behaviour of edgeR of reference
			# sample not being 1
			if("reference" %in% ( (.) %>% dplyr::pull(sample) ))
				(.) %>%
				dplyr::mutate(
					multiplier =
						multiplier /
						(.) %>%
						dplyr::filter(sample == "reference") %>%
						dplyr::pull(multiplier)
				)
			else
				(.)
		}

	# Output
	input.df %>% 
		arrange(!!as.symbol(sample_column), !!as.symbol(gene_column)) %>%
		
		# Add normalised data set
		bind_cols(
			df %>%
				dplyr::mutate(sample = as.factor(as.character(sample))) %>%
				dplyr::left_join( nf, by = "sample") %>%
				
				# Calculate normalised values
				dplyr::mutate(value_normalised = value * multiplier) %>%
				
				# Format df for join
				dplyr::select(sample, gene, value_normalised, everything()) %>%
				dplyr::mutate(filt_for_calc = gene %in% nf_obj$gene_to_exclude) %>%
				dplyr::select(-value, -tot, -tot_filt) %>%
				dplyr::rename(TMM = nf) %>%
				setNames(c("sample", "gene", sprintf("%s normalised" ,value_column), colnames(.)[4:ncol(.)])) %>%
				arrange(sample, gene) %>%
				dplyr::select(-sample, -gene)
		)
}

