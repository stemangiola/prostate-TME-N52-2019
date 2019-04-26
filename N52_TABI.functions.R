print_summary_return = function(., other_print = ""){ (.) %>% summary() %>% print(); (.) }
print_return = function(.){ (.) %>% print(); (.) }

# Get annotation
age <- function(dob, age.day = lubridate::today(), units = "years", floor = TRUE) {
	calc.age = lubridate::interval(dob, age.day) / lubridate::duration(num = 1, units = units)
	if (floor) return(as.integer(floor(calc.age)))
	return(calc.age)
}

make_stats_and_plots = function(){
	stats_STAR = do.call("rbind", lapply(dir(path = sprintf("%s/alignment_hg38", my_data_dir), recursive=T, pattern="Log\\.final\\.out$", full.names = T), function(i){
		df <- suppressMessages(suppressWarnings(
			read_delim(i,  "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
		))
		
		df$X1 = sapply(df$X1, function(x) gsub(" |", "", x, fixed=T ))
		colnames(df) = c("info", "value")
		df$sample = strsplit(i, "/|_")[[1]][5]
		df$`N adapter` = strsplit(i, "_|-")[[1]][4]
		df$`S adapter` = strsplit(i, "_|-")[[1]][5]
		df$reference_genome = strsplit(i, "/")[[1]][1]
		
		df
	}))
	stats_STAR = stats_STAR %>% filter(info%in%c("Number of input reads", "Uniquely mapped reads number", "Number of splices: Total", "Number of reads mapped to multiple loci"))
	stats_STAR$value = as.numeric(stats_STAR$value)
	
	ggplot(stats_STAR %>% filter(info%in%c("Number of input reads", "Uniquely mapped reads number", "Number of splices: Total", "Number of reads mapped to multiple loci")), aes(x = factor(sample), y = value, fill=info)) +
		geom_bar(stat="identity",position ="identity") +
		scale_fill_brewer(palette="Set1") +
		facet_grid(reference_genome~.,scales = "free_x") +
		theme_bw() +
		theme(
			panel.border = element_blank(), 
			axis.line = element_line(),
			panel.grid.major = element_line(size = 0.2),
			panel.grid.minor = element_line(size = 0.1),
			text = element_text(size=12),
			legend.position="bottom",
			aspect.ratio=1,
			axis.text.x = element_text(angle = 90, hjust = 1),
			strip.background = element_blank(),
			axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
			axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
			
			
		) 
	ggsave(
		filename =	"STAR_stats.pdf",
		device = "pdf",
		useDingbats=FALSE,
		units = c("mm"),
		width = 183 ,
		height = 183 / 3
	)
	
	stats_RNAseqc = do.call("rbind", lapply(dir(path = sprintf("%s/alignment_hg38", my_data_dir), recursive=T, pattern="metrics\\.tmp\\.txt$", full.names = T), function(i){
		df <- suppressMessages(suppressWarnings(
			read_delim(i, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
		))
		df.formatted = tibble(
			info=unlist(df[c(1,3,5,7,9),]), 
			value=as.numeric(unlist(df[c(2,4,6,8,10),])), 
			sample=strsplit(basename(i), "_|-")[[1]][7],
			`N adapter` = strsplit(basename(i), "_|-")[[1]][5],
			`S adapter` = strsplit(basename(i), "_|-")[[1]][6],
			reference_genome  = strsplit(i, "/")[[1]][1]
		)
		df.formatted
	}))
	
	manipulate = function(x, name_orig, name_alt){
		mapped_count = (x %>% filter(info=="Mapped Pairs"))$value
		my_row = x %>% filter(info==name_orig) %>% mutate(info=name_alt) %>% mutate(value = mapped_count * value)
		x %>% rbind(my_row)
	}
	stats_RNAseqc = stats_RNAseqc %>% group_by(sample, reference_genome) %>% do(manipulate(., "Exonic Rate", "Exonic Count")) %>% ungroup()
	stats_RNAseqc = stats_RNAseqc %>% group_by(sample, reference_genome) %>% do(manipulate(., "Intragenic Rate", "Intragenic Count")) %>% ungroup()
	
	
	subread_stats = sam.counts$stat
	colnames(subread_stats)[1] = "info"
	#colnames(subread_stats)[-1]  = sapply(colnames(subread_stats)[-1], function(cn) strsplit(cn, "\\.|_")[[1]][7])
	subread_stats = as_tibble(reshape::melt(subread_stats))
	subread_stats$`N adapter` = sapply(as.character(subread_stats$variable), function(n) strsplit(n, "_|\\.")[[1]][5])
	subread_stats$`S adapter` = sapply(as.character(subread_stats$variable), function(n) strsplit(n, "_|\\.")[[1]][6])
	subread_stats$reference_genome = "alignment_hg38"
	colnames(subread_stats)[2] = "sample"
	
	ggplot(subread_stats %>% filter(info%in%c("Assigned", "Unassigned_NoFeatures", "Unassigned_FragmentLength")), aes(x = factor(sample), y = value, color=info)) +
		geom_bar(stat="identity",position ="identity", alpha=0) + 
		scale_color_brewer(palette="Set1") +
		theme_bw() +
		theme(
			panel.border = element_blank(), 
			axis.line = element_line(),
			panel.grid.major = element_line(size = 0.2),
			panel.grid.minor = element_line(size = 0.1),
			text = element_text(size=12),
			legend.position="bottom",
			aspect.ratio=1,
			axis.text.x = element_text(angle = 90, hjust = 1),
			axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
			axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
			
			
		)
	
	ggsave(
		filename =	"assignment_stats.pdf",
		device = "pdf",
		useDingbats=FALSE,
		units = c("mm"),
		width = 183 ,
		height = 183 / 3
	)
	
	
	exp_annot = tibble(
		file = colnames(d),
		`N adapter` = sapply(colnames(d), function(dd) strsplit(dd, "\\.|_")[[1]][5]),
		`S adapter` = sapply(colnames(d), function(dd) strsplit(dd, "\\.|_")[[1]][6])
	)
	library(gsheet)
	
	facs_annot = gsheet2tbl("https://docs.google.com/spreadsheets/d/19N9JS3AZAeBMeUrGQZo76GkERI_wKEvDWTMbr2mVgUU/edit#gid=0", sheetid = "diary")
	facs_annot$DOS = as.Date(facs_annot$date,format="%d/%m/%Y") 
	facs_annot$cell_counts = as.numeric(facs_annot$cell_counts)
	facs_annot$Batch = as.factor(as.numeric(facs_annot$Batch))
	seq_annot = gsheet2tbl("https://docs.google.com/spreadsheets/d/19N9JS3AZAeBMeUrGQZo76GkERI_wKEvDWTMbr2mVgUU/edit#gid=1898907242", sheetid = "sequencing")
	seq_annot$DOS = as.Date(seq_annot$DOS,format="%d/%m/%Y") 
	exp_seq_annot = inner_join(exp_annot, seq_annot, by = c("N adapter", "S adapter"))
	clinical_annot <- read_csv("/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/PhD/TME/N52_plus_12_run/fastq/run_all/clinical_info.csv")
	clinical_annot$DOS = as.Date(clinical_annot$DOS,format="%d/%m/%Y") 
	clinical_annot$gleason_total = paste(clinical_annot$PathGG1, clinical_annot$PathGG2, sep="")
	capra_annot <- read_csv("/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/PhD/TME/N52_plus_12_run/fastq/run_all/capra_N52_NMC.csv")
	capra_annot$DOS = as.Date(capra_annot$DOS,format="%d/%m/%y") 
	annot = clinical_annot
	annot = left_join(annot, exp_seq_annot,  by = c("DOS"))
	
	annot = left_join(annot, stats_STAR %>% spread(info, value), by = c("N adapter", "S adapter")) 
	annot = left_join(annot, stats_RNAseqc  %>% spread(info, value), by = c("N adapter", "S adapter")) 
	annot = left_join(annot, subread_stats %>% spread(info, value), by = c("N adapter", "S adapter")) 
	annot = left_join(annot, facs_annot,  by = c("DOS", "cell type"))
	annot = left_join(annot, capra_annot,  by = c("DOS", "RARP", "UBR"))
	
	colnames(annot)[colnames(annot)=="sample.x"] = "sample"
	
	#subsample only 
	annot = annot[match(colnames(d), annot$file),]
	
	annot$counts_sum = colSums(d$counts)
	annot$density_pathological = 0
	annot[sapply(sprintf("_%s.", c("S28", "S31", "S30", "S29", "S14", "S16", "S17", "S8", "S32", "S13", "S18", "S15", "S27")), function(s) grep(s, annot$file)),"density_pathological"] = 1
	annot$`cell_type_formatted` = sapply(annot$`cell type`, substr, 1, 1)
	annot$CAPRA_TOTAL_orig = annot$CAPRA_TOTAL
	annot[annot$label=="benign", "CAPRA_TOTAL"] = 0
	annot$gleason_total = as.numeric(annot$gleason_total)
	
	annot$is_mono_neutro = annot$`cell type` == "M"
	
	annot$CAPRA_groups = "low"
	annot$CAPRA_groups[annot$CAPRA_TOTAL > 2] = "high"
	
	annot$CAPRA_strips = 0
	annot$CAPRA_strips[annot$CAPRA_TOTAL%in%3:5] = 1
	annot$CAPRA_strips[annot$CAPRA_TOTAL>=6] = 2
	
	annot$batch = as.factor(annot$batch)
	
	# Set age
	annot = 
		annot %>%
		mutate(age = age(format(as.Date(DOB.x, "%d/%m/%Y"), "%Y-%m-%d"))) %>%
		mutate(age_rel = age - min(age))
	
	# # Plots
	library(RColorBrewer)
	n <- 10
	qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
	col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
	set.seed(1273)
	col_vector = sample(col_vector)
	
	p =
		annot %>% dplyr:::select(
			c(
				"cell_type_formatted" ,
				"sample", 
				"Number of input reads",
				"Mapped Pairs",
				"Uniquely mapped reads number", 
				#"Number of splices: Total", 
				#"Number of reads mapped to multiple loci",
				#"Assigned",
				"Exonic Count",
				#"Intragenic Count",
				#"Alternative Aligments",
				"Unassigned_NoFeatures"
			)) %>% 
		setNames(
			c(
				"Cell type",
				"Sample",
				"Input reads",
				"Mapped reads",
				"Uniquely mapped",
				"Exonic reads",
				"Unasigned reads"
			)
		) %>%
		gather(info, value, -Sample, -`Cell type` ) %>%
		mutate(info = as.factor(info)) %>%
		mutate(info = factor(info, levels=c(
			"Cell type",
			"Sample",
			"Input reads",
			"Mapped reads",
			"Uniquely mapped",
			"Exonic reads",
			"Unasigned reads"
		))) %>%
		ggplot(aes(x = factor(Sample), y = value, color=info, group=info)) +
		geom_line() +
		facet_grid(~`Cell type`, scale="free") +
		#scale_color_manual(values = col_vector) + 
		scale_color_brewer(palette="Set1", direction = 1) +
		theme_bw() +
		theme(
			panel.border = element_blank(), 
			axis.line = element_line(),
			panel.grid.major = element_line(size = 0.2),
			panel.grid.minor = element_line(size = 0.1),
			text = element_text(size=12),
			legend.position="bottom",
			legend.direction = "vertical",
			aspect.ratio=1,
			axis.text.x = element_text(angle = 90, hjust = 1),
			strip.background = element_blank(),
			axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
			axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
			
		) + xlab("Sample") + ylab("Read count")
	
	ggsave(
		plot=p, 
		filename =	"mapping_stats.pdf",
		device = "pdf",
		useDingbats=FALSE,
		units = c("mm"),
		width = 183 ,
		height = 183 
	)
	
	list(plot=p, annot=annot) 
	
}


make_pairs = function(df){
	ggpairs(
		df , 
		lower = list(continuous = wrap("points", alpha = 0.5), 
								 combo = wrap("box"), 
								 discrete = wrap("facetbar") ), 
		diag = list(continuous = wrap("densityDiag",  color = "blue", alpha = 0.5) )
	) +
		theme_bw() +
		theme(
			panel.border = element_blank(), 
			axis.line = element_line(),
			panel.grid.major = element_line(size = 0.2),
			panel.grid.minor = element_line(size = 0.1),
			text = element_text(size=12),
			legend.position="bottom",
			aspect.ratio=1,
			strip.text.x = element_text(angle = 90), strip.text.y = element_text(angle = 0),
			axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
			axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
			
		)
	
}

make_pair_plot_features = function(){
	
	features_to_observe = c(
		"Number of input reads", 
		#"label", 
		#"cell_counts" ,
		#"Batch", 
		"CAPRA_TOTAL", 
		"cell_type_formatted" 
	)
	features_to_observe_scale1 = c("Mapped")
	features_to_observe_scale2 = c(
		"Uniquely mapped reads number"
		#"Exonic Count",
		#"Unassigned_Ambiguity", 
		#"Unassigned_NoFeatures"
		#	"counts_sum", 
		#"Alternative Aligments",
		
	)
	
	 
	
	categories_ordered = 
		c(features_to_observe, 
			features_to_observe_scale1, 
			features_to_observe_scale2,
			"Bacterial reads proportion")
	
	df_4_pairs = annot %>% dplyr::select(categories_ordered)
	df_4_pairs[,colnames(df_4_pairs)%in%features_to_observe_scale1] = as_tibble(df_4_pairs[,colnames(df_4_pairs)%in%features_to_observe_scale1]/annot$`Number of input reads`/2)
	df_4_pairs[,colnames(df_4_pairs)%in%features_to_observe_scale2] = as_tibble(df_4_pairs[,colnames(df_4_pairs)%in%features_to_observe_scale2]/annot$Mapped)
	colnames(df_4_pairs) = sapply(colnames(df_4_pairs) , function(x) gsub(" ", "_", x))
	
	
	# Format colanames
	ggplot <- function(...) ggplot2::ggplot(...) 
	firstup <- function(x) {
		substr(x, 1, 1) <- toupper(substr(x, 1, 1))
		x
	}
	
	
	ggplot =	function(...) ggplot2::ggplot() +	scale_color_brewer(palette="Dark2") +		scale_fill_brewer(palette="Dark2")
	ggplot <- function(...) ggplot2::ggplot(...) + scale_fill_brewer(palette="Set2")
	
	ggally_mysmooth <- function(data, mapping, ...){
		ggplot(data = data, mapping=mapping) +
			geom_density(mapping = aes(color=`Cell type`), fill=NA)
	}
	
	p = df_4_pairs %>%
		setNames(c("Input reads", "CAPRA", "Cell type", "Mapped reads", "Uniquely mapped", "Bacterial reads proportion")) %>%
		select(c("Cell type", "CAPRA", "Input reads",  "Mapped reads", "Uniquely mapped", "Bacterial reads proportion")) %>%
		#setNames(firstup(gsub("_", " ", colnames(.)))) %>%
		ggpairs( 
			#my_ggplot,
			mapping = aes(color=`Cell type`),
			lower = list(
				continuous = wrap("points"), 
				combo = wrap("box"), 
				discrete = wrap("facetbar") 
			), 
			diag = list(continuous = ggally_mysmooth )
		) +
		theme_bw() + 
		theme(
			panel.border = element_blank(), 
			axis.line = element_line(),
			panel.grid.major = element_line(size = 0.2),
			panel.grid.minor = element_line(size = 0.1),
			text = element_text(size=12),
			legend.position="bottom",
			#aspect.ratio=1,
			strip.text.x = element_text(angle = 90), 
			strip.text.y = element_text(angle = 0),
			strip.background = element_blank(), axis.text.x = element_text(angle = 90),
			axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
			axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
		)
	
	p[1, 1] <- NULL  # set plot at position 1,1 to not appear
	p[1, 2] <- NULL  # set plot at position 1,1 to not appear
	p[1, 3] <- NULL  # set plot at position 1,1 to not appear
	p[1, 4] <- NULL  # set plot at position 1,1 to not appear
	p[1, 5] <- NULL  # set plot at position 1,1 to not appear
	p[1, 6] <- NULL  # set plot at position 1,1 to not appear

	
	scales <- scale_colour_brewer(palette = 'Dark2')  
	scales2 = scale_fill_brewer(palette = 'Dark2') 
	# 
	for (row in seq_len(p$nrow))
		for (col in seq_len(p$ncol)){
			p[row, col] <- p[row, col] + scales
			p[row, col] <- p[row, col] + scales2
		}
	p
	
	ggsave(
		plot=p, 
		filename =	"pairs_stats.pdf",
		device = "pdf",
		useDingbats=FALSE,
		units = c("mm"),
		width = 183 ,
		height = 183 
	)
	
}

# Average duplicated
average_duplicated_genes_tibble_spreaded = function(tbl){
	
	dup_genes =
		tbl %>%
		dplyr::group_by(symbol) %>%
		dplyr::summarise(tot = n()) %>%
		dplyr::filter(tot > 1) %>%
		dplyr::pull(symbol) %>%
		as.character()
	
	tbl %>%
		dplyr::mutate(symbol = as.character(symbol)) %>%
		dplyr::filter(!symbol %in% dup_genes) %>%
		dplyr::bind_rows(
			tbl %>%
				dplyr::mutate(symbol = as.character(symbol)) %>%
				dplyr::filter(symbol %in% dup_genes) %>%
				dplyr::group_by(symbol) %>%
				dplyr::summarise_if(is.numeric, median) %>%
				dplyr::ungroup()
		)
	
}

bezierCurve <- function(x, y, n=10)
{
	# http://rosettacode.org/wiki/Bitmap/B%C3%A9zier_curves/Cubic
	outx <- NULL
	outy <- NULL
	
	bez <- function(x, y, t)
	{
		outx <- 0
		outy <- 0
		n <- length(x)-1
		for (i in 0:n)
		{
			outx <- outx + choose(n, i)*((1-t)^(n-i))*t^i*x[i+1]
			outy <- outy + choose(n, i)*((1-t)^(n-i))*t^i*y[i+1]
		}
		
		return (list(x=outx, y=outy))
	}
	
	i <- 1
	for (t in seq(0, 1, length.out=n))
	{
		b <- bez(x, y, t)
		outx[i] <- b$x
		outy[i] <- b$y
		
		i <- i+1
	}
	
	return (list(x=outx, y=outy))
}

collect_inflections = function(ct, run_directory){
	foreach( fi = dir(path=run_directory, pattern = sprintf("tabi_res_%s_", ct), full.names = T), .combine = bind_rows) %dopar% {
		load(fi)
		
		# Add inflection 
		tabi_res$fit %>% spread_samples(inflection[gene_idx], log_y_cross[gene_idx], beta[cov_idx, gene_idx]) %>%
			ungroup() %>%
			filter(cov_idx == 2) %>%
			dplyr::select(-cov_idx) %>%
			dplyr::rename(slope = beta) %>%
			mutate(gene_idx = as.integer(gene_idx)) %>%
			
			# Add gene symbol back
			left_join(
				tabi_res$posterior_df %>%
					ungroup() %>% 
					distinct(gene_idx, gene)
			) 
	}
}

collect_res = function(ct, run_directory){
	foreach( fi = dir(path=run_directory, pattern = sprintf("tabi_res_%s_", ct), full.names = T), .combine = bind_rows) %dopar% {
		load(fi)
		
		tabi_res$posterior_df %>% 
			
			# Do stats
			dplyr::select(-gene)  %>% 
			mean_qi() %>%
			
			# Add gene symbol back
			left_join(tabi_res$posterior_df %>% distinct(gene_idx, gene)) %>%
			
			# Add inflection point estimate
			left_join(
				rstan::summary(tabi_res$fit, c("inflection", "log_y_cross"))$summary %>% 
					as_tibble(rownames="parameter") %>% 
					separate(parameter, c("term", "gene_idx"), sep = "[\\[\\],]") %>%
					dplyr::select(term ,  gene_idx, mean, `2.5%`, `97.5%`) %>% 
					gather(temp, score, -gene_idx, -term) %>%
					unite(temp1, term, temp, sep = ".") %>%
					spread(temp1, score) %>%
					setNames( gsub(".mean", "", (.) %>% colnames(), fixed=T) ) %>%
					mutate(gene_idx = as.integer(gene_idx))
			) %>%
			
			# Add factors of interest
			left_join(
				rstan::summary(tabi_res$fit, c("beta"))$summary %>% 
					as_tibble(rownames="parameter") %>% 
					separate(parameter, c("term", "covariate_idx", "gene_idx"), sep = "[\\[\\],]") %>%
					mutate(covariate = paste("b", as.integer(covariate_idx)-1, sep="")) %>%
					dplyr::select(covariate,  gene_idx, mean, `2.5%`, `97.5%`) %>% 
					gather(temp, score, -covariate, -gene_idx) %>%
					unite(temp1, covariate, temp, sep = ".") %>%
					spread(temp1, score) %>%
					setNames( gsub(".mean", "", (.) %>% colnames(), fixed=T) ) %>%
					mutate(gene_idx = as.integer(gene_idx))
			)
		
		
	}
}

collect_generated_quantities = function(ct, run_directory){
	foreach( fi = dir(path=run_directory, pattern =  sprintf("tabi_res_%s_", ct), full.names = T), .combine = bind_rows) %dopar% {
		load(fi)
		tabi_res$generated_quantities %>% 
			
			# Add gene symbol back
			left_join(tabi_res$posterior_df %>% ungroup() %>% distinct(gene_idx, gene)) %>%
			
			# Add count info
			left_join(
				tabi_res$input.y %>% mutate(sample_idx=1:n()) %>% gather(gene, `read count`, -sample_idx) ,
				by=c("sample_idx", "gene")
			) %>%
			
			# Add covariates
			left_join(
				tabi_res$input.X %>% mutate(sample_idx=1:n()) %>% dplyr::select(-`(Intercept)`),
				by="sample_idx"
			)
	}
}


get_log_inflection = function(y_cross, inflection, slope){
	#log1p_exp = function(w)	log(1+exp(w))
	#print(sprintf("k = %s", y_cross + log1p_exp(-( -inflection * slope ))))
	# y_cross + log1p_exp(-( -inflection * slope )) - log1p_exp(- ( -inflection * slope + slope * x ) ) 
	
	#http://www.wolframalpha.com/input/?i=solve+c+%2B+log(1%2Be%5E-(-n*b))+-+log(1%2Be%5E-(-n*b%2Bb*x))+%3D+(c+%2B+log(1%2Be%5E-(-n*b)))%2F2+for+x
	(slope * inflection - log(exp(y_cross/2) * sqrt(exp(slope * inflection) + 1) - 1))/slope
}

rescale_inflection = function(inflection){
	CAPRA_score = c(1, 0, 3, 3, 6, 0, 1, 3, 2, 0, 1, 7, 4)
	(inflection * sd(CAPRA_score)) + mean(CAPRA_score)
}

annotate_res = function(res_df){
	
	####################################### 
	# Offset
	# Put the zero on 3.5 CAPRA score
	# x = annot %>%
	# 	filter(cell_type_formatted==ct) %>%
	# 	arrange(file) %>%
	# 	pull(CAPRA_TOTAL)
	CAPRA_score = c(1, 0, 3, 3, 6, 0, 1, 3, 2, 0, 1, 7, 4)
	offset_x = (2 - mean(CAPRA_score)) / sd(CAPRA_score)
	
	#######################################
	
	
	# biomaRt::getBM(
	# 	filters=c("hgnc_symbol"),
	# 	values = "VIT", # (.)$symbol , 
	# 	attributes = c(
	# 		'COSMIC'
	# 	),
	# 	mart = biomaRt::useEnsembl(biomart="snp", dataset = "hsapiens_snp_som")
	# )
	
	res_df %>% 
		ungroup() %>%
		filter(covariate_idx == 2) %>% 
		filter(conf.low * conf.high > 0) %>% 
		
		# Correct inflection
		# Some inflection 2.5% are NaN that's why the transform produces NaNs
		
		mutate(
			inflection = get_log_inflection(log_y_cross, inflection, estimate),
			`inflection.2.5%` = get_log_inflection(log_y_cross, `inflection.2.5%`, estimate),
			`inflection.97.5%` = get_log_inflection(log_y_cross, `inflection.97.5%`, estimate)
		) %>%
		
		#mutate(inflection = inflection - !!offset_x) %>%
		mutate(
			inflection = rescale_inflection(inflection),
			`inflection.2.5%` = rescale_inflection(`inflection.2.5%`),
			`inflection.97.5%` = rescale_inflection(`inflection.97.5%`)
		) %>%
		
		# Annotate
		mutate(when = ifelse(inflection<2, "early", "late"), how = ifelse(estimate > 0, "up", "down")) %>%
		
		# Join secretome information
		dplyr::rename(symbol=gene) %>%
		left_join(
			bind_rows(
				read_delim("secretome/NOT.tsv", "\t",  escape_double = FALSE, trim_ws = TRUE),
				read_delim("secretome/protein_class_Predicted.tsv", "\t",  escape_double = FALSE, trim_ws = TRUE),
				read_delim("secretome/protein_class_Predicted (1).tsv", "\t",  escape_double = FALSE, trim_ws = TRUE),
				read_delim("secretome/protein_class_Predicted (2).tsv", "\t",  escape_double = FALSE, trim_ws = TRUE),
				read_delim("secretome/protein_class_Predicted (3).tsv", "\t",  escape_double = FALSE, trim_ws = TRUE)
			) %>% distinct() %>% dplyr::rename(symbol = Gene) %>%
				mutate(`Protein class` = strsplit(`Protein class`, ",")) %>%
				unnest(`Protein class`) %>%
				mutate(`Protein class` = gsub("^ ", "", `Protein class`)), 
			by = "symbol"
		) %>%
		
		
		# Join GO information
		#load("BM.RData")
		left_join(
			biomaRt::getBM(
				filters=c("hgnc_symbol"),
				values = (.)$symbol , 
				attributes = c(
					'ensembl_gene_id', 
					'entrezgene','hgnc_symbol',
					'description',
					'name_1006',    'namespace_1003', 'definition_1006'
				),
				mart = biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
			) %>%
				as_tibble() %>%
				dplyr::rename(symbol=hgnc_symbol),
			by = "symbol"
		)
}

# Filter too many 0s
filter_too_sparse = function(tbl, gen_df, prop = 13/3){
	tbl %>% left_join(
		gen_df %>% 
			mutate(`non zero` = as.integer(`read count`>1)) %>%
			group_by(gene) %>%
			summarise(`non zero` = sum(`non zero`)) %>% 
			dplyr::rename(symbol=gene)
	) %>%
		do({
			
			# Print
			# print(
			# 	(.) %>% 
			# 		distinct(symbol, conf.low, conf.high, `non zero`) %>% 
			# 		filter(conf.low * conf.high > 0) %>%
			# 		filter(`non zero` <= ceiling(prop)) %>% nrow
			# )
			
			# Return
			(.) %>% filter(`non zero` >= ceiling(prop))
		})
	
}

# Filter too manh outside CI
filter_out_CI = function(tbl, gen_df, prop = 13/3){
	tbl %>% left_join(
		gen_df %>% 
			group_by(gene) %>% 
			summarise(n_out_CI = sum( (conf.low-`read count`) * (conf.high-`read count`) > 0 )) %>%
			dplyr::rename(symbol=gene)
	) %>%
		do({
			
			# Print
			# print(
			# 	(.) %>% 
			# 		distinct(symbol, conf.low, conf.high, n_out_CI) %>% 
			# 		filter(conf.low * conf.high > 0) %>%
			# 		filter(n_out_CI >= floor(prop)) %>%
			# 		nrow
			# )
			
			# Return
			(.) %>% filter(n_out_CI < floor(prop))
		})
	
}

print_summary = function(res_df, gen_df, annot_df){
	
	# Summary statistics
	writeLines("Table of cell-cell interface genes")
	annot_df %>%
		filter_too_sparse(gen_df) %>%
		filter_out_CI(gen_df) %>%
		#filter(grepl("secreted|membrane", `Protein class`)) %>%
		distinct(when, how, symbol) %>%
		group_by(when, how) %>% 
		tally() %>%
		spread(how, n) %>%
		print()
	
	# Plot drug targets
	gen_df %>%
		left_join(
			# Select drug targets
			annot_df %>%
				# Filter genes that have too many 0s
				filter_too_sparse() %>%
				filter_out_CI() %>%
				filter(grepl("secreted|membrane", `Protein class`)) %>%
				mutate(`Safe estimate` = ifelse(how=="up", conf.low, conf.high)) %>%
				distinct(when, how, symbol, `Safe estimate`) %>%
				group_by(when, how) %>% 
				arrange(`Safe estimate` %>% abs %>% desc) %>%
				
				# Select top genes
				filter(row_number()%in%1:5 ) %>%
				mutate(`target` = 1) %>%
				dplyr::rename(gene=symbol) %>% 
				dplyr::select(gene, when, how, target)
		) %>%
		filter(`target` == 1) %>% 
		distinct() %>%
		group_by(when, how) %>%
		do(
			plots = 
				(.) %>%
				ggplot(aes(x=CAPRA_TOTAL, y=`read count`, color=gene)) +
				geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width = 0) +
				geom_jitter() +
				geom_vline(xintercept = 0, linetype="dashed") +
				facet_wrap(~gene, scale="free") +
				scale_y_log10() +
				theme_bw() +
				theme(
					strip.background = element_blank(), 
					axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
					axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
					
				) +
				ggtitle((.) %>% distinct(when, how) %>% as.data.frame() %>% paste(collapse=" ")  )
		) %>%
		
		# Compose the plots
		pull(plots) %>%
		gridExtra::grid.arrange(grobs=., top="Top cell-cell interface genes" ) 
	
	# Print first classes
	foreach(go_class = c("biological_process") )  %:% foreach(h = c("up", "down")) %do% {
		annot_df %>% 
			filter_too_sparse(gen_df) %>%
			filter_out_CI(gen_df) %>%
			#filter(grepl("secreted|membrane", `Protein class`)) %>%
			filter(namespace_1003 == go_class) %>% 
			filter(how==h) %>% 
			# Take last word
			rowwise() %>%
			mutate(name_1006_last = 
						 	# Temporary function TAKE LAST ELEMENT
						 	(function(...) {
						 		rev(strsplit(..., " ")[[1]])[1] %>% 
						 		{
						 			ifelse(
						 				grepl("^negative", ..., ignore.case = T), 
						 				paste(strsplit(..., " ")[[1]], .),
						 				.
						 			)
						 		}
						 	})(name_1006)
			) %>%
			ungroup() %>%
			# Annotate with how many non 0s
			left_join(
				gen_df %>% 
					mutate(`non zero` = as.integer(`read count`>0)) %>%
					group_by(gene) %>%
					summarise(`non zero` = sum(`non zero`)) %>% 
					dplyr::rename(symbol=gene)
			) %>%
			filter(`non zero` > ceiling(13/3)) %>%
			# Organise data set for grouping
			distinct(namespace_1003, how,symbol, name_1006, name_1006_last) %>% 
			# group_by(namespace_1003, when, how, name_1006_last) %>% 
			# summarise(n = n(), genes = paste(symbol, collapse = ','), name_1006 = paste(name_1006, collapse = ',')) %>% 
			group_by(namespace_1003, how, name_1006) %>%
			summarise(n = n(), genes = paste(symbol, collapse = ',')) %>%
			rowwise() %>%
			mutate(n = strsplit(genes, ",")[[1]] %>% unique() %>% length() ) %>%
			arrange(desc(n)) 
	}
}

print_gene_info = function(annot_df, gen_df, inflection_df, prop_filter_CI = 13/3, do_filter = T){
	browser()
	# Summary statistics
	annot_df %>%
	{
		if(do_filter)
			(.) %>% 
			filter_too_sparse(gen_df) %>%
			filter_out_CI(gen_df, prop_filter_CI) 
		else
			(.)
	}	%>%
		filter(namespace_1003 == "biological_process") %>%
		
		# Safe estimate
		mutate(`Safe estimate` = ifelse(how=="up", conf.low, conf.high)) %>%
		
		# Arrange
		arrange(desc(abs(`Safe estimate`)), symbol) %>%
		mutate(symbol = factor(symbol, levels = (.) %>% pull(symbol) %>% unique())) %>%
		
		# Loop genes
		group_by(symbol) %>%
		distinct() %>%
		do({
			
			# View
			(.) %>% dplyr::select(symbol, estimate, `Safe estimate`, name_1006, definition_1006) %>% View()
			
			# Print
			(.) %>% dplyr::select(symbol, estimate, `Safe estimate`, name_1006, definition_1006) %>% print(n=30)
			
			temp = (.)
			plot(
				gen_df %>% filter(gene== (temp %>% pull(symbol) %>% unique() %>% as.character())) %>%
					
					ggplot(aes(x=CAPRA_TOTAL, y=`read count`, color=gene)) +
					geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width = 0) +
					geom_jitter() +
					geom_vline(xintercept = 0, linetype="dashed") +
					facet_wrap(~gene, scale="free") +
					scale_y_log10() +
					theme_bw() +
					theme(
						aspect.ratio=1, 
						strip.background = element_blank(),
						axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
						axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
						
						
					) +
					ggtitle((.) %>% distinct(when, how) %>% as.data.frame() %>% paste(collapse=" ")  )
			)
			browser()
			
			
		})
	
}

core_gene_regression_plot = function(tbl, inflection_df, gen_df, multiple_ct = F, rescale_infl = F){

	if(multiple_ct)
		m_theme = theme(
			aspect.ratio=1, 				
			axis.title.x  = element_blank(),
			axis.title.y  = element_blank(),
			plot.title = element_text(size = 3),
			plot.margin = unit(c(0, 0, 0, 0), "cm")
		)
	else
		m_theme = theme(
			aspect.ratio=1, 				
			axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
			axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
			plot.title = element_text(size = 3),
			plot.margin = unit(c(0, 0, 0, 0), "cm")
		)
	
	
	pmap(
		inflection_df %>% 
			mutate(b0 = - inflection * slope) %>% 
			#rename(symbol=gene) %>% 
			filter(symbol == tbl %>% pull(symbol) %>% as.character %>% unique) %>% 
			sample_frac(0.5) %>%
			select(slope, b0, log_y_cross),
		function(x, slope, b0, log_y_cross) {
			stat_function(
				data = data.frame(x=c(-1,2)),
				fun = function(x, slope, b0, log_y_cross) {
					exp(log_y_cross) * (1 + exp(-b0)) / ( 1 + exp(- ( b0 + x * slope ) ) ) + 1;
				},
				args = list(
					slope = slope, 
					b0 = b0, 
					log_y_cross = log_y_cross
				) , 
				geom="line",
				alpha=0.05,
				color = "gray45"
			)
		}
	) %>% 
		reduce(
			.init = (
				gen_df %>% 
					filter(	symbol== (tbl %>% pull(symbol) %>% unique() %>% as.character())	) %>% 
					mutate(`read count` = `read count`+1) %>%
					ggplot() 
			),
			`+`
		) +
		geom_jitter(
			aes(x=CAPRA_TOTAL, y=`read count`), 
			color = ifelse(
				multiple_ct,
					switch(
						tbl %>% separate(symbol, c("Cell type", "g")) %>% distinct(`Cell type`) %>% pull(`Cell type`),
						"E" = "#199e77",
						"F" = "#d86012",
						"T" = "#e42e89",
						"M" = "#7570b2"
					),
				"#db2423"
				), 
			size=2
		) +
		geom_errorbar(aes(x=CAPRA_TOTAL,ymin=conf.low, ymax=conf.high), width = 0.2, alpha=ifelse(multiple_ct, 0.5, 0.8)) +
		scale_y_log10() +
		scale_color_brewer(palette="Dark2") +
		theme_bw() +
		m_theme +
		ggtitle(tbl %>% pull(symbol) %>% as.character %>% unique)
	
}

print_gene_regression = function(
	
		annot_df, 
		gen_df, 
		inflection_df, 
		prop_filter_CI = 13/3, 
		do_filter = T, 
		multiple_ct = F, 
		num_cols = NULL,
		gene_order = NULL
	){
		
		gen_inv_logit = function(x, slope, b0, log_y_cross) {
			exp(log_y_cross) * (1 + exp(-b0)) / ( 1 + exp(- ( b0 + x * slope ) ) );
		}
		
		# Summary statistics
		annot_df %>%
		{
			if(do_filter)
				(.) %>% 
				filter_too_sparse(gen_df) %>%
				filter_out_CI(gen_df, prop_filter_CI) 
			else
				(.)
		}	%>%
			
			distinct(symbol, conf.low, conf.high) %>%
			
			# Add read counts
			left_join(	gen_df %>% distinct(symbol, `read count`, CAPRA_TOTAL) 	) %>%
			
			# Order if necessarily
			{
				if(gene_order %>% is.null %>% `!`)
					(.) %>% mutate(symbol = factor(symbol, levels = gene_order)) 
				else 
					(.)
			} %>%
			# Loop genes
			group_by(symbol) %>%
			#distinct() %>%
			#multidplyr::partition(symbol) %>%
			
			do({
				
				`%>%` = magrittr::`%>%`
				library(tidyverse)
				library(magrittr)
				source("N52_TABI.functions.R")
					
				tibble(
					symbol = (.) %>% distinct(symbol) %>% pull(symbol),
					plots = core_gene_regression_plot(
						(.) %>% distinct() , 
						inflection_df, 
						gen_df, 
						multiple_ct
					) %>% list()
				)
			}) %>%
			#collect() %>%
			pull(plots) %>%
			#gridExtra::arrangeGrob( grobs=., ncol = 1 ) %>%
			cowplot::plot_grid(plotlist = ., align = "v",  axis="b", rel_widths = 1 , ncol = ifelse(num_cols %>% is.na %>% `!`, num_cols, annot_df %>% pull(symbol) %>% unique %>% length %>% sqrt %>% ceiling))
		
	}

plot_inflections = function(
	gen_df_E, annot_df_E, inflec_df_E, 
	gen_df_F, annot_df_F, inflec_df_F, 
	gen_df_T, annot_df_T, inflec_df_T, 
	gen_df_M, annot_df_M, inflec_df_M
){
	
	join_my_df = function(gen_df, annot_df, inflec_df){
		annot_df %>%
			filter_too_sparse(gen_df) %>%
			filter_out_CI(gen_df) %>%
			distinct(symbol, estimate, log_y_cross, b1) %>%
			
			# Add whole posterior
			inner_join(	
				inflec_df %>% 
					sample_frac(0.2) %>%
					distinct(gene, inflection) %>%
					dplyr::rename(symbol = gene)
			) %>%
			
			# Rescale
			rowwise() %>%
			mutate(inflection = get_log_inflection(log_y_cross, inflection, b1)) %>%	 
			mutate(inflection = rescale_inflection(inflection)) %>%	 
			ungroup()
	}
	
	# Plot distrbution of inflections
	foreach(
		my_list = 
			list( 
				E = list(gen_df_E, annot_df_E, inflec_df_E),
				F = list(gen_df_F, annot_df_F, inflec_df_F),
				T = list(gen_df_T, annot_df_T, inflec_df_T),
				M = list(gen_df_M, annot_df_M, inflec_df_M)
			),
		.combine = bind_rows,
		ct = c("E", "F", "T", "M")
	) %dopar% {
		join_my_df(my_list[[1]], my_list[[2]], my_list[[3]]) %>% 
			mutate(`Cell type` = ct)
	} %>%
		ggplot(aes(inflection, color=`Cell type`, linetype = estimate>0)) + 
		annotate("rect", xmin = 0, xmax = 7, ymin = 0, ymax = Inf,alpha = .2) +
		geom_density() +
		xlim(-25, 25) +
		facet_grid(`Cell type`~., scales = "free") +
		scale_color_brewer(palette="Dark2") +
		theme_bw() +
		theme(
			panel.border = element_blank(), 
			axis.line = element_line(),
			panel.grid.major = element_line(size = 0.2),
			panel.grid.minor = element_line(size = 0.1),
			text = element_text(size=8),
			legend.position="bottom",
			strip.background = element_blank(),
			axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
			axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
			
			
		) +
		ylab("density") + xlab("CAPRA score")
	
}

# Set scales for plot

signed_log = scales::trans_new(
	"signed_log",
	transform=function(x) sign(x)*log(abs(x)),
	inverse=function(x) sign(x)*exp(abs(x))
)

signed_log_breaks = scales::trans_breaks(function(x) sign(x)*log(abs(x)), function(x) sign(x)*exp(abs(x)),6)

my_trans <- function(x, i_low = interval_low, i_high = interval_high, r = reducer) {
	sapply(x, function(x) {
		if(x >= i_low & x <= i_high ) x
		else if(x < i_low) x / r + i_low
		else (x - i_high) / r + i_high
	})
}

magnify_trans <- function(interval_low = 0, interval_high = 7,  reducer = 20) {
	
	trans <- function(x, i_low = interval_low, i_high = interval_high, r = reducer) {
		sapply(x, function(x) {
			#print(sprintf("- %s",x))
			if(x >= i_low & x <= i_high ) x
			else if(x < i_low) x / r + i_low
			else (x - i_high) / r + i_high
		})
	}
	
	inv <- function(x, i_low = interval_low, i_high = interval_high, r = reducer) {
		sapply(x, function(x) {
			#print(sprintf("-- %s",x))
			if(!is.na(x)) {
				if(x >= i_low & x <= i_high ) x
				else if(x < i_low) (x - i_low) * r
				else ((x - i_high) * r) + i_high
			}
			else print("NA in transformation of y axis in ggplot")
		})
	}
	
	trans_new(name = 'custom',	transform = trans,inverse = inv	)
}


# Set functions for plot

log_gen_inv_logit = function( y_log,  b0,  log_y_cross) {
	log1p_exp = function(x) log( 1 + exp(x) )
	log_y_cross + log1p_exp(-b0) - log1p_exp(- y_log  ) ;
}

effect_size = function(b0, b1, log_y_cross, cap_to_one = F) {
	CAPRA0_standatdised = -1.06
	CAPRA7_standardised = 2.05
	delta_array = c(
		log_gen_inv_logit( b0 + b1 * CAPRA0_standatdised,  b0,  log_y_cross),
		log_gen_inv_logit( b0 + b1 * CAPRA7_standardised,  b0,  log_y_cross)
	) %>% exp()
	
	max(delta_array) /
		ifelse(	min(delta_array) < 1 & cap_to_one, 1, min(delta_array))  * 
		ifelse( delta_array[1]<delta_array[2], 1, -1)
}

cusotm_root_trans = function() scales::trans_new("cusotm_root",function(x) x^(1/4), function(x) x^(4))


# Calculate CAPRA-S score
capra_s = function(
	psa, 
	sur_mar, 
	semin_ves_invas, 
	gleason, 
	extracap_ext,
	lymph_inv, 
	sample
) {
	psa_score = function(psa){
		psa = as.numeric(psa)
		if(psa <= 6) 0
		else if(psa >6 & psa <= 10) 1
		else if(psa >10 & psa <=20) 2
		else if(psa > 20) 3
	}
	
	gleason_score = function(gleason){
		gleason = as.numeric(gleason)
		if(gleason %>% gsub("\\+", "", .) <= 33) 0
		else if(gleason %>% gsub("\\+", "", .) == 34) 1
		else if(gleason %>% gsub("\\+", "", .) == 43) 2
		else if(gleason %>% gsub("\\+", "", .) > 43) 3
	}
	
	sum(
		psa_score(psa) ,
		switch(sur_mar %>% as.character, "0" = 0, "1" = 2) ,
		switch(semin_ves_invas %>% as.character, "0" = 0, "1" = 2) ,
		gleason_score(gleason),
		switch(extracap_ext %>% as.character, "0" = 0, "1" = 2) ,
		switch(lymph_inv %>% as.character, "0" = 0, "1" = 2), na.rm = T
	)
}

get_redundant_genes_kmeans = function(
	input.df,
	replicates_column = "replicates",
	cluster_by_column = "cluster_by",
	value_column = "value",
	components_list = list(c(1,2)),
	log_transform = F,
	number_of_clusters,
	add_MDS_components
	
){

	# If number of genes is less than nuber of clusters bypass everything
	if(input.df %>% distinct(symbol) %>% nrow %>% `<=` (number_of_clusters))	
		return( input.df %>%distinct(symbol) %>%mutate(is_redundant = F))
	
	input.df %>%
		
		# Derive components
		add_MDS_components(
			replicates_column = replicates_column,
			cluster_by_column = cluster_by_column,
			value_column = value_column, 
			log_transform = T,
			components_list = components_list
		) %>%
		
		# Calculate clusters
		add_kmeans_clusters(
			replicates_column = "Component",
			cluster_by_column = "symbol",
			value_column = "Component value", 
			number_of_clusters = number_of_clusters
		) %>%
		
		# # Plot
		# distinct(symbol, cluster, Component, `Component value`) %>%
		# spread(Component, `Component value`) %>%
		# {
		# 	(.) %>%
		# 	ggplot(aes(x = `1`, y = `2`, color = factor(cluster), label=symbol)) + 
		# 	geom_point() 
		# } %>% 
		# ggplotly()
		
	distinct(symbol, `Cell type`, cluster, `Unique score`) %>%
		
		# Add scores
		group_by(cluster) %>%
		arrange(`Unique score` %>% desc) %>%
		mutate(is_redundant =  row_number() > 1) %>%
		ungroup() %>%
		distinct(symbol, is_redundant) 
	
}

get_redundant_genes_cor = function(
	input.df,
	replicates_column = "file",
	cluster_by_column = "symbol",
	value_column = "read count",
	log_transform = T,
	cor_threshold = 0.9
) 

{
	
	input.df %>%
		select(!!replicates_column, !!cluster_by_column, !!value_column) %>%
		distinct %>%
		{
			if(log_transform) 
				(.) %>% mutate(	!!value_column := !!as.name(value_column) %>% `+` (1) %>% log)
			else 
				(.)
		} %>%
		spread(!!replicates_column, !!value_column) %>%
		as_matrix(rownames = cluster_by_column) %>%
		t %>%
		cor %>%
		as_tibble(rownames = "symbol_1") %>%
		gather(symbol_2, correlation, -symbol_1) %>%
		filter(symbol_1 != symbol_2) %>%
		
		# Add scores
		left_join(
			input.df %>%
				distinct(symbol, `Unique score`) %>%
				rename(symbol_1 = symbol, `Unique score 1` = `Unique score`)
		) %>%
		left_join(
			input.df %>%
				distinct(symbol, `Unique score`) %>%
				rename(symbol_2 = symbol, `Unique score 2` = `Unique score`)
		) %>%
		mutate(symbol = ifelse(`Unique score 1` > `Unique score 2`, symbol_2, symbol_1)) %>%
		mutate(is_redundant = ifelse(correlation > cor_threshold, T, F)) %>%
		filter(is_redundant) %>%
		distinct(symbol, is_redundant) %>%
		mutate(`Cell type` = input.df %>% head(n=1) %>% pull(`Cell type`)) %>%
		
		# Add back the non redundant
		bind_rows({
			temp_df = (.)
			input.df %>% 
				distinct(symbol, `Cell type`) %>%
				unite(symbol_ct, c("symbol", "Cell type"), remove = F) %>%	
				filter(! symbol_ct %in% ( temp_df %>% unite(symbol_ct, c("symbol", "Cell type")) %>% pull(symbol_ct))) %>%
				select(-symbol_ct)
		}) %>%
		mutate(is_redundant = ifelse(is_redundant %>% is.na, FALSE, is_redundant))
		
}

get_TCGA_prostate_clinical_annotaton = function(){
	foreach(
		fi = 
			mycgds %>% 
			cgdsr::getCancerStudies() %>% 
			filter(grepl("TCGA", name)) %>% 
			pull(cancer_study_id) %>% 
			grep("prad", ., value = T), 
		.combine = bind_rows
	) %dopar% {
		print(fi)
		mycaselist = cgdsr::getCaseLists(mycgds,fi) %>% as_tibble() %>% filter(grepl("RNA Seq", case_list_name ) )
		print(ncol(mycaselist))
		cgdsr::getClinicalData(mycgds,caseList = mycaselist[1,] %>% pull(case_list_id)) %>%
			as_tibble(rownames="sample") %>%
			mutate_all(funs(as.character)) %>%
			mutate(cancer_study_id = fi)
	} %>%
		# Clean from absent data
		filter(! (DFS_MONTHS == "" | is.na(DFS_MONTHS))) %>%
		# Select duplicate more up to date
		group_by(sample) %>%
		arrange(desc(DFS_MONTHS)) %>%
		filter(row_number()==1) %>%
		ungroup() %>%
		separate(sample, c("t1", "t2", "t3"), sep = "\\.") %>%
		unite(sample, c("t1", "t2", "t3"), sep = "-")  %>%
		
		#3 Calculate CAPRA-S score
		unite(GLEASON_SCORE_FORMATTED, c("GLEASON_PATTERN_PRIMARY", "GLEASON_PATTERN_SECONDARY"), sep="", remove=F) %>%
		mutate(	sur_mar = ifelse(RESIDUAL_TUMOR == "R0", 0, 1)) %>%
		mutate(LYMPH_NODES_EXAMINED_HE_COUNT = ifelse(is.na(LYMPH_NODES_EXAMINED_HE_COUNT), 0, LYMPH_NODES_EXAMINED_HE_COUNT)) %>%
		filter(!is.na(PSA_MOST_RECENT_RESULTS) & !is.na(GLEASON_SCORE_FORMATTED)) %>%
		rowwise() %>%
		mutate(
			`CAPRA-S` = capra_s(
				psa = PSA_MOST_RECENT_RESULTS, 
				sur_mar = sur_mar, 
				semin_ves_invas = 0, 
				gleason = GLEASON_SCORE_FORMATTED, 
				extracap_ext = 0,
				lymph_inv = LYMPH_NODES_EXAMINED_HE_COUNT, sample=sample
			)
		) %>%
		ungroup() %>%
		
		# Format further
		distinct(
			sample, 
			DFS_MONTHS,
			DFS_STATUS,
			`CAPRA-S`
		) %>%
		mutate(is_recurred = ifelse(DFS_STATUS=="Recurred/Progressed", T, F)) %>%
		mutate(DFS_MONTHS = as.numeric(DFS_MONTHS))
}
