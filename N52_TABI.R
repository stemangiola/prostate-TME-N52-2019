
library(tidyverse)
library(magrittr)
library(foreach)
#library(plotly)
library(rstan)
library(tidybayes)
library(doParallel)
registerDoParallel()
source("https://raw.githubusercontent.com/betanalpha/knitr_case_studies/master/qr_regression/stan_utility.R")
# For plotting
library(scales)
library(googledrive)
library(googlesheets)
my_url = "https://docs.google.com/spreadsheets/d/1IJPo_ZAW2cxawrGOoWU9JaU77ipxiFFpF69ytxPJE_U/edit?usp=sharing"
source("N52_TABI.functions.R")
library(edgeR)
# Plot graphics
library(ggrepel)
library(grid)
source("utilities_normalization.R")
library(sva)
library(GGally)
# as_matrix function
source("https://gist.githubusercontent.com/stemangiola/dd3573be22492fc03856cd2c53a755a9/raw/cd5621d6eaf78c6706929f0abccbc8bc124eb9f8/tidy_extensions.R")
source("https://gist.githubusercontent.com/stemangiola/90a528038b8c52b21f9cfa6bb186d583/raw/4d5fa8dfd52c7687d803608c3e71ace1b71a2038/transcription_tool_kit.R")

my_theme = 	
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

#############################################################
# Setup data directory from 
my_data_dir = '..'
#############################################################


#############################################################
# Load data #################################################

# Get counts
load("counts_paired_end.RData")
d <- edgeR::DGEList(counts=sam.counts$counts, genes=sam.counts$annotation[,c("GeneID","Length")])
d$genes$symbol <- AnnotationDbi:::mapIds(org.Hs.eg.db::org.Hs.eg.db,
																				 keys=as.character(d$genes$GeneID),
																				 column="SYMBOL",
																				 keytype="ENTREZID",
																				 multiVals="first")

# Plots
stat = make_stats_and_plots()
annot = stat$annot %>% mutate_if(is.character, as.factor)




ex = as_tibble(d$counts, rownames="GeneID") %>%
	mutate(GeneID = as.integer(GeneID)) %>%
	left_join(
		as_tibble(d$genes) %>% dplyr::select(-Length),
		by = "GeneID"
	) %>%
	dplyr::select(-GeneID)

# Save stats
save(list=c("ex", "annot"), file="input_parallel_TABI.RData")

#############################################################
# MDS plots #################################################

foreach(ct = c("E", "F", "M", "T"), .combine = bind_rows) %do% 
{
	
	# collect expression data
	d$counts[,which(annot$cell_type_formatted == ct)] %>% 
		as_tibble %>%
		mutate(symbol = d$genes$symbol) %>%
		gather(sample, value, -symbol) %>%
		drop_na() %>%
		
		# Normalise
		norm_RNAseq( 
			sample_column = "sample", 
			gene_column = "symbol", 
			value_column = "value"
		) %>%
		filter(!filt_for_calc) %>%
		dplyr::select(symbol, sample, `value normalised`) %>%
		mutate(`value normalised log` = log(`value normalised` + 1)) %>%
		
		# Annotate
		left_join(
			annot %>% 
				distinct(file, ClinStageT, CAPRA_groups, CAPRA_TOTAL, batch, sample) %>% 
				dplyr::rename(sample_label = sample) %>% 
				dplyr::rename(sample = file) %>% 
				mutate(ClinStageT = as.character(ClinStageT)) %>%
				mutate(ClinStageT = ifelse(is.na(ClinStageT), "T1c", ClinStageT))
		) %>%
		mutate_if(is.character, as.factor) %>%
		
		# Batch correction
		{
			my_df = (.)
			
			my_df %>% 
				dplyr::select(symbol, sample, `value normalised log`) %>%
				spread(sample, `value normalised log`) %>%
				{
					mat = (.) %>% dplyr::select(-symbol) %>% as.matrix
					rownames(mat) = (.) %>% pull(symbol)
					mat
				} %>%
				ComBat(
					batch=my_df %>% distinct(sample, batch) %>% pull(batch), 
					mod=model.matrix(~ my_df %>% distinct(sample, CAPRA_groups) %>% pull(CAPRA_groups))
				) %>%
				{
					
					if(ct=="F")
						(.) %>%
						ComBat(
							batch=my_df %>% distinct(sample, ClinStageT) %>% pull(ClinStageT), 
							mod=model.matrix(~ my_df %>% distinct(sample, CAPRA_groups) %>% pull(CAPRA_groups))
						)
					else (.)
				} %>%
				
				# MDS calculation
				{
					
					my_df_mds = (.)
					
					foreach(
						components = list(c(1, 2), c(3, 4), c(5, 6), c(7, 8)), 
						#my_df_mds = (.),
						.combine = bind_rows
					) %do% {
						my_df_mds %>%
							plotMDS(dim.plot = components) %>% 
							{
								tibble(sample = names((.)$x), x = (.)$x, y = (.)$y, PCx = components[1], PCy = components[2])
							}
					} %>%
						left_join(my_df %>% distinct(sample, ClinStageT, batch, CAPRA_groups, CAPRA_TOTAL, sample_label)) %>%
						mutate(ct = ct)
				}
		}
} %>%
	
	# Plot and save
{
	ggplot(data=(.), aes(x = x, y = y, label = sample_label, PCx = PCx, PCy = PCy)) + 
		geom_point(aes(fill = CAPRA_TOTAL), size=3, shape=21, color="grey20") +
		ggrepel::geom_text_repel(
			size = 1.6, 
			point.padding = 0.3, 
			#fontface = 'bold', 
			# label.padding = 0.1, 
			# label.size = 0,
			segment.size = 0.2,
			seed = 654321
		) +
		#geom_text(color = "grey20", size = 2 ) +
		scale_fill_distiller(
			palette = "Spectral",
			na.value = 'white'
		) +
		facet_grid(interaction(PCx, PCy)~ct) +
		theme_bw() +
		theme(
			panel.border = element_blank(), 
			axis.line = element_line(),
			panel.grid.major = element_line(size = 0.2),
			panel.grid.minor = element_line(size = 0.1),
			text = element_text(size=12),
			legend.position="bottom",
			aspect.ratio=1,
			strip.background = element_blank(),
			axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
			axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
			
			
		) +
		xlab("Principal component x") + ylab("Principal component y") 
} %>%
	ggsave(plot = .,
				 sprintf("mds_plot.pdf"),
				 useDingbats=FALSE,
				 units = c("mm"),
				 width = 183 ,
				 height = 183 
	)


#############################################################
# Sample composititon #######################################


library(ARMET)

# ARMET_wrap = function(ct){
# 	ARMET_tc(
# 		mix = 
# 			d$counts[,which(annot$cell_type_formatted == ct)] %>% 
# 			as_tibble %>%
# 			mutate(gene = d$genes$symbol) %>%
# 			dplyr::select(gene, everything()) %>%
# 			filter(!is.na(gene)),
# 		
# 		my_design =
# 			model.matrix(
# 				~
# 					annot %>% filter(cell_type_formatted == ct) %>% pull(CAPRA_TOTAL) %>% scale + 
# 					annot %>% filter(cell_type_formatted == ct) %>% pull(batch) 
# 			) %>%
# 			as.data.frame() %>%
# 			as_tibble()	 %>%
# 			mutate(sample = annot %>% filter(cell_type_formatted == ct) %>% pull(file) ) %>%
# 			setNames(c("(Intercept)", "CAPRA_TOTAL", "batch", "sample")),
# 		cov_to_test = "CAPRA_TOTAL", do_debug = F, save_fit = T, verbose = T
# 	)
# }

# Save figures
foreach(ct = c("E", "F", "M", "T")) %do% {
	load(sprintf("%s/dtc_%s.RData", my_data_dir, ct))
	dtc %>%
		ARMET_plotPolar(
			size_geom_text = 1.8, 
			my_breaks=c(0, 0.01, 0.1,0.5,1),
			prop_filter = 0.01,
			barwidth = 0.5, barheight = 3,
			legend_justification = 0.78
		) +
		ggtitle(ct)
} %>%
	gridExtra::grid.arrange(grobs=.) %>%	
	ggsave(plot = .,
				 "dtc_polar.pdf",
				 useDingbats=FALSE,
				 units = c("mm"),
				 width = 183 ,
				 height = 183 
	)


# Cibersort
source("CIBERSORT_annotated.R")

foreach(ct = c("E", "F", "M", "T"), .combine = bind_rows) %dopar% {
	
	d$counts[,which(annot$cell_type_formatted == ct)] %>% 
		as_tibble %>%
		mutate(symbol = d$genes$symbol) %>%
		gather(sample, value, -symbol) %>%
		drop_na() %>%
		norm_RNAseq( 
			sample_column = "sample", 
			gene_column = "symbol", 
			value_column = "value"
		) %>%
		select(symbol, sample, `value normalised`) %>%
		spread(sample, `value normalised`) %>%
		{
			file_name = sprintf("expression_for_cibersort_%s.tab", ct)
			write_delim((.), path = file_name, delim = "\t")
			file_name
		} %>%
		{
			CIBERSORT(
				"LM22.txt",
				(.)
			)$proportions %>%
				as_tibble %>% 
				select(-`P-value`, -Correlation, -RMSE) %>% 
				gather(`Cell type`, value) %>% 
				group_by(`Cell type`) %>% 
				summarise(m = mean(value)) %>%
				mutate(`Group` = ct)
		}
}

#############################################################
# Run DE ####################################################

#############################################################
#  unix311 524 % parallel --ungroup --linebuffer 'Rscript run_parallel_TABI.R E {}' ::: 1 2 3 4 5
#  unix311 524 % parallel --ungroup --linebuffer 'Rscript run_parallel_TABI.R F {}' ::: 1 2 3 4 5
#  unix311 524 % parallel --ungroup --linebuffer 'Rscript run_parallel_TABI.R M {}' ::: 1 2 3 4 5
#  unix311 524 % parallel --ungroup --linebuffer 'Rscript run_parallel_TABI.R T {}' ::: 1 2 3 4 5
#############################################################

# tabi_res_E = collect_res("E", sprintf("%s/tabi_res_11_07_2018", my_data_dir))
# tabi_inflection_E = collect_inflections("E", sprintf("%s/tabi_res_11_07_2018", my_data_dir))
# tabi_gen_E = collect_generated_quantities("E", sprintf("%s/tabi_res_11_07_2018", my_data_dir))
# tabi_res_E.annot = annotate_res(tabi_res_E)
# 
# tabi_res_T = collect_res("T", sprintf("%s/tabi_res_11_07_2018", my_data_dir))
# tabi_inflection_T = collect_inflections("T", sprintf("%s/tabi_res_11_07_2018", my_data_dir))
# tabi_gen_T = collect_generated_quantities("T", sprintf("%s/tabi_res_11_07_2018", my_data_dir))
# tabi_res_T.annot = annotate_res(tabi_res_T)
# 
# tabi_res_M = collect_res("M", sprintf("%s/tabi_res_11_07_2018", my_data_dir))
# tabi_inflection_M = collect_inflections("M", sprintf("%s/tabi_res_11_07_2018", my_data_dir))
# tabi_gen_M = collect_generated_quantities("M", sprintf("%s/tabi_res_11_07_2018", my_data_dir))
# tabi_res_M.annot = annotate_res(tabi_res_M)
# 
# tabi_res_F = collect_res("F", sprintf("%s/tabi_res_11_07_2018", my_data_dir))
# tabi_inflection_F = collect_inflections("F", sprintf("%s/tabi_res_11_07_2018", my_data_dir))
# tabi_gen_F = collect_generated_quantities("F", sprintf("%s/tabi_res_11_07_2018", my_data_dir))
# tabi_res_F.annot = annotate_res(tabi_res_F)
# 
# save(
# 	list = c(
# 		"tabi_res_E", "tabi_inflection_E", "tabi_gen_E", "tabi_res_E.annot",
# 		"tabi_res_F", "tabi_inflection_F", "tabi_gen_F", "tabi_res_F.annot",
# 		"tabi_res_M", "tabi_inflection_M", "tabi_gen_M", "tabi_res_M.annot",
# 		"tabi_res_T", "tabi_inflection_T", "tabi_gen_T", "tabi_res_T.annot"
# 	), file="TABI_formatted_results.RData"
# )

load("TABI_formatted_results.RData")

#############################################################
# Plot all genes ############################################

tabi_res_E %>% ungroup %>% rename(symbol=gene) %>% filter_too_sparse(tabi_gen_E) %>% filter_out_CI(tabi_gen_E) %>% mutate(`Cell type` = "E")  %>% 
bind_rows(	tabi_res_F %>% rename(symbol=gene) %>% ungroup %>% filter_too_sparse(tabi_gen_F) %>% filter_out_CI(tabi_gen_F) %>% mutate(`Cell type` = "F")) %>%
bind_rows(	tabi_res_T %>% rename(symbol=gene) %>% ungroup %>% filter_too_sparse(tabi_gen_T) %>% filter_out_CI(tabi_gen_T) %>% mutate(`Cell type` = "T")) %>%
bind_rows(	tabi_res_M %>% rename(symbol=gene) %>% ungroup %>% filter_too_sparse(tabi_gen_M) %>% filter_out_CI(tabi_gen_M) %>% mutate(`Cell type` = "M")) %>%
	mutate(DE = conf.low * conf.high > 0) %>%
	#sample_frac(0.1) %>%
	arrange(estimate) %>%
	group_by(`Cell type`) %>%
	mutate(symbol = 1:n()) %>%
	ungroup() %>%
	{
		ggplot((.), aes(x=symbol, color=DE)) +
			geom_errorbar(aes(ymin=conf.low, ymax=conf.high, alpha = DE), width = 0) +
			facet_grid(`Cell type`~., scales = "free_x") +
			theme_bw() +
			theme(
				panel.border = element_blank(), 
				axis.line = element_line(),
				panel.grid.major = element_line(size = 0.2),
				panel.grid.minor = element_line(size = 0.1),
				text = element_text(size=12),
				legend.position="bottom",
				axis.text.x = element_blank(),
				axis.ticks.x = element_blank(),
				strip.background = element_blank(),
				axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
				axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
				
				
			) +
			scale_colour_manual(values = c("grey", "#db2323")) +
			scale_alpha_manual(values=c(0.1, 1)) + 
			xlab("Genes") + 
			ylab("Estimate")
	} 

ggsave(
	plot = (.),
	filename =	"DE_summary.pdf",
	device = "pdf",
	useDingbats=FALSE,
	units = c("mm"),
	width = 183 ,
	height = 183 
)

#############################################################
# Print filtered genes for erticle ##########################

dummy=tabi_res_E %>% ungroup() %>% dplyr::rename(symbol=gene) %>% filter_out_CI(tabi_gen_E)
dummy=tabi_res_T %>% ungroup() %>% dplyr::rename(symbol=gene) %>% filter_out_CI(tabi_gen_T)
dummy=tabi_res_F %>% ungroup() %>% dplyr::rename(symbol=gene) %>% filter_out_CI(tabi_gen_F)
dummy=tabi_res_M %>% ungroup() %>% dplyr::rename(symbol=gene) %>% filter_out_CI(tabi_gen_M)

#############################################################
# Print DE genes ############################################

tabi_res_E.annot %>% filter_out_CI(tabi_gen_E) %>% filter_too_sparse(tabi_gen_E) %>% distinct(symbol, estimate ) %>% group_by(estimate > 0) %>% summarise(n = n())
tabi_res_T.annot %>% filter_out_CI(tabi_gen_T) %>% filter_too_sparse(tabi_gen_T) %>% distinct(symbol, estimate ) %>% group_by(estimate > 0) %>% summarise(n = n())
tabi_res_F.annot %>% filter_out_CI(tabi_gen_F) %>% filter_too_sparse(tabi_gen_F) %>% distinct(symbol, estimate ) %>% group_by(estimate > 0) %>% summarise(n = n())
tabi_res_M.annot %>% filter_out_CI(tabi_gen_M) %>% filter_too_sparse(tabi_gen_M) %>% distinct(symbol, estimate ) %>% group_by(estimate > 0) %>% summarise(n = n())

#############################################################
# Print prostate cancer genes ###############################

read_sheet(ss = my_url, sheet = "E") %>% 
	bind_rows(
		read_sheet(ss = my_url, sheet = "E_intracellular") %>%
			filter(is.na(`bad fit`))
	) %>%
	dplyr::select(gene, `Prostate related`) %>%
	dplyr::rename(symbol=gene) %>%
	left_join(
		tabi_res_E.annot %>% dplyr::distinct(symbol, `Gene description`, estimate, `Protein class`)
	) %>%
	group_by(symbol, `Gene description`, estimate, `Prostate related`) %>%
	summarise(`Protein class` = paste(`Protein class`, collapse=", ")) %>%
	ungroup() %>%
	arrange(is.na(`Prostate related`)) %>%
	{
		write_csv((.), "tabi_prostate_cancer_related_genes.csv")
		
		(.) %>% 
			filter(!is.na( `Prostate related`)) %>% 
			group_by(estimate * `Prostate related` > 0) %>% 
			summarise(n = n()) %>% 
			{ (.) %>% print(); (.) } %>%
			
			# Enrichment test
			{ 
				prop.test(
					table(
						c ( 
							rep("pro", (.)[2,] %>% pull(n)),
							rep("anti", (.)[1,] %>% pull(n)) 
						)
					),
					alternative = "less"
				)$p.value %>% 
					scientific()
			}
		
	}

#############################################################
# Plot inflections ##########################################

plot_inflections(
	tabi_gen_E, tabi_res_E.annot, tabi_inflection_E, 
	tabi_gen_F, tabi_res_F.annot, tabi_inflection_F,
	tabi_gen_T, tabi_res_T.annot, tabi_inflection_T,
	tabi_gen_M, tabi_res_M.annot, tabi_inflection_M
) %>% 
	ggsave(filename  = "inflection_densities.pdf", useDingbats=FALSE, units = c("mm"), width = 89 , height = 120)


#############################################################
# Create summary data frame #################################

summary_df = foreach(
	ct_list = 
		list( 
			E = list(tabi_res_E.annot, tabi_inflection_E),
			F = list(tabi_res_F.annot, tabi_inflection_F),
			T = list(tabi_res_T.annot, tabi_inflection_T),
			M = list(tabi_res_M.annot, tabi_inflection_M)
		),
	.combine = bind_rows, .verbose = T,
	ct = c("E", "F", "T", "M")
) %dopar% {
	ct_list[[1]] %>%
		distinct(
			symbol, 
			log_y_cross, `log_y_cross.2.5%`, `log_y_cross.97.5%`, 
			b0, `b0.2.5%`, `b0.97.5%`, 
			b1, `b1.2.5%`, `b1.97.5%`
		) %>%
		inner_join(
			read_sheet(ss = my_url, sheet = ct) %>%
				gather(Category, value, 1:(which((.) %>% names() == "gene") - 1)) %>%
				filter(!is.na(value)) %>%
				dplyr::rename(symbol = gene) %>%
				mutate(`Cell type` = ct) %>%
				
				# Convert 0 to both +1 and -1
				{
					(.) %>% 
						mutate(value = ifelse(value == 0, 1, value)) %>%
						bind_rows(
							(.) %>% 
								filter(value == 0) %>% 
								mutate(value = -1)
						)
				}
			
		) %>%
		
		# Add inflection posterior
		inner_join(	
			ct_list[[2]] %>% 
				distinct(gene, inflection) %>%
				dplyr::rename(symbol = gene), 
			by = "symbol"	
		)
	
} %>%
	
	# Calculate Effect size
	rowwise() %>%
	mutate(`Effect size` = effect_size(b0, b1, log_y_cross)) %>%
	ungroup() %>%
	
	# Flip effect size if inhibitor
	#mutate(`Effect size` = `Effect size` * value)	%>%
	
	# Correct inflection
	rowwise() %>%
	mutate(inflection = get_log_inflection(log_y_cross, inflection, b1)) %>%	 
	mutate(inflection = rescale_inflection(inflection)) %>%	 
	ungroup() %>%
	filter(inflection != "NaN") %>%
	
	# Factorise
	mutate(`Cell type` = factor(`Cell type`)) %>%
	
	# Distinguish same symbol of different cell types
	mutate(symbol_ct = interaction(symbol, `Cell type`) ) %>%
	
	# Arrange on inflection
	do({
		my_order = (.) %>% 
			distinct(symbol_ct, inflection) %>%
			group_by(symbol_ct) %>%
			summarise(mean_inflection = median(inflection)) %>%
			arrange(mean_inflection) %>% 
			pull(symbol_ct)
		(.) %>% mutate(	symbol_ct = factor(symbol_ct, levels = my_order))
	}) %>%
	
	#Add adjustment fot violin size
	left_join( read_sheet(ss = my_url, sheet = "SCALING FACTOR") ) %>%
	
	# See if pro or anti tumor
	mutate(positively_associated = `Effect size` / abs(`Effect size`) * value ) %>%
	mutate(positively_associated = ifelse(abs(positively_associated)>1, 2, positively_associated))

#############################################################
# Heatmap ###################################################

summary_df %>%
	distinct(Category, symbol, `Cell type`, inflection) %>%
	group_by(Category, symbol, `Cell type`) %>%
	summarise(inflection = inflection %>% median) %>%
	mutate(inflection = rescale_inflection(inflection)) %>%	 
	ungroup() %>%
	#rename(`Cell type original` = `Cell type`) %>%
	left_join(
		ex %>%
			filter(symbol %>% is.na %>% `!`)  %>%
			gather(sample, `read count`, -symbol) %>%
			norm_RNAseq( 
				sample_column = "sample", 
				gene_column = "symbol", 
				value_column = "read count"
			) %>%
			left_join(
				annot %>% 
					select(file, cell_type_formatted, CAPRA_TOTAL, UBR) %>%
					rename(`Cell type` = cell_type_formatted, sample = file) 
			)
	) %>%
	mutate(`read count normalised log` =  `read count normalised` %>% `+` (1) %>% log)%>%
	mutate(sample = gsub("alignment_hg38.merged_", "", sample)) %>%
	mutate(sample = gsub("_trimmedBBduk.Aligned.out.sam", "", sample)) %>%
	unite(sample_ct, `Cell type`, sample, remove = F) %>%
	unite(symbol_ct, `Cell type`, symbol, remove = F) %>%
	group_by(Category) %>%
	do({
		
		tbl = (.) 
		
		mat = tbl %>%
			select(UBR, `read count normalised log`, symbol_ct) %>%
			distinct() %>%
			spread( UBR, `read count normalised log`) %>%
			drop_na() %>%
			as_matrix(rownames = "symbol_ct") %>%
			t() %>%
			apply(2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))) %>%
			#scale() %>%
			t() 
		
		robust_dist = function(x, y) {
			qx = quantile(x, c(0.1, 0.9))
			qy = quantile(y, c(0.1, 0.9))
			l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
			x = x[l]
			y = y[l]
			sqrt(sum((x - y)^2))
		}
		
		ct_colors = function(ct) foreach(cc = ct, .combine = c) %do% switch(
			cc,
			"E" = "#199E78",
			"F" = "#D96013",
			"M" = "#7571B3",
			"T" = "#E52E89"
		)
		
		ct_names = function(ct) foreach(cc = ct, .combine = c) %do% switch(
			cc,
			"E" = "Epithelial",
			"F" = "Fibroblast",
			"M" = "Myeloid",
			"T" = "T cell"
		)
		
		pdf(
			sprintf("heatmap_pathway_interface_%s.pdf", tbl %>% pull(Category) %>% unique %>% gsub(" |/", "_", .)), 
			height = 0.5 * tbl %>% distinct(symbol) %>% nrow %>% `sum` (3), 
			useDingbats=F
		)
		
		mat %>% 
			Heatmap(
				column_title = tbl %>% head(n=1) %>% pull(Category),
				row_title = NULL,
				name="mat", 
				width = unit(0.5 * 13, "cm"),
				height = unit(0.5 * tbl %>% distinct(symbol) %>% nrow, "cm"),
				col = circlize::colorRamp2(c(-2, 0, 2), viridis::viridis(3)),
				column_dend_reorder = tbl %>% distinct(UBR, CAPRA_TOTAL) %>% arrange(UBR) %>% pull(CAPRA_TOTAL) %>% `+` (1),
				row_split = 
					tbl %>% distinct(symbol_ct, `Cell type`) %>% 
					arrange(symbol_ct) %>% pull(`Cell type`), 
				cluster_row_slices = FALSE,
				#	clustering_distance_columns = robust_dist,
				left_annotation =
					rowAnnotation(
						ct = anno_block(
							gp = gpar(fill = ct_colors(
								tbl %>% distinct(`Cell type`) %>% arrange(`Cell type`) %>% pull(`Cell type`) %>% as.character
							)), 
							labels = ct_names(
								tbl %>% distinct(`Cell type`) %>% arrange(`Cell type`) %>% pull(`Cell type`) %>% as.character
							), 
							labels_gp = gpar(col = "white")
						),
						
						inflection =  anno_points(
							tbl %>% distinct(symbol_ct, inflection) %>% 
								arrange(symbol_ct) %>% pull(inflection)
						)
						
					),
				top_annotation  =
					HeatmapAnnotation(
						`CAPRA-S` =  tbl %>% distinct(UBR, CAPRA_TOTAL) %>% 
							arrange(UBR) %>% pull(CAPRA_TOTAL) ,
						col = list(	`CAPRA-S`  = circlize::colorRamp2(0:7, colorRampPalette(RColorBrewer ::brewer.pal(11,"Spectral") %>% rev)(8)))
					)
			) %>% print()
		
		dev.off()
		
		tbl %>% distinct(Category)
	}	)



#############################################################
# Plot summary ##############################################

# Plot for each category
summary_df %>%
	group_by(Category) %>%
	do({
		
		# Arrange symbol_ct
		p = (.) %>% 
			
			# Plot
			ggplot( aes(x =symbol_ct , y=inflection, fill=`Effect size`)) +
			geom_violin(
				scale = "width", adjust = 1, width = 1
				
			) +
			scale_fill_distiller(
				palette = "Spectral",
				na.value = 'white',
				direction = 1,
				trans = signed_log,
				breaks = signed_log_breaks,
				labels = scales::scientific,
				limits=c(
					-max(abs((.) %>% pull(`Effect size`))),
					max(abs((.) %>% pull(`Effect size`)))
				)
			) +
			geom_point( 
				data = (.) %>% 
					group_by(symbol_ct, Category, positively_associated) %>% 
					summarise(
						i = quantile(inflection, 0.5), 
						`Effect size` = unique(`Effect size`)
					),
				aes(y = i, x = symbol_ct),
				size = 0.5
			) +
			geom_boxplot(
				data = (.) %>% 
					distinct(
						symbol_ct, 
						`Effect size`,
						`Cell type`, 
						inflection, 
						Category,
						positively_associated
					) %>% 
					mutate(lower = -20, upper = -19),
				aes(x=symbol_ct, ymin = lower, lower = lower, middle = lower, upper = upper, ymax = upper, color=`Cell type`),
				stat = "identity"
			) +
			scale_color_brewer(palette="Dark2") +
			coord_flip(
				ylim=c(-25,17),
				xlim=c(1, (.) %>% distinct(symbol_ct, positively_associated) %>% group_by(positively_associated) %>% summarise(n = n()) %>% pull(n) %>% max )
			)+
			facet_wrap(
				~ positively_associated, 
				scales = "free", 
				#space = "free",
				drop = T
			) +
			scale_y_continuous(
				trans = magnify_trans(interval_low = 0, interval_high = 7,  reducer = 20),
				#trans = magnify_trans(intercept =  7,  reducer = 20),
				breaks = c(-10, 0, 1:7, 17)
			) +
			theme_bw() +
			theme(
				panel.border = element_blank(), 
				axis.line = element_line(),
				panel.grid.major = element_line(size = 0.2),
				panel.grid.minor = element_line(size = 0.1),
				text = element_text(size=12),
				legend.position="bottom",
				legend.text = element_text(angle = 60, vjust = 0.5),
				axis.text.x = element_text(angle = 60, vjust = 0.5),
				strip.background = element_blank(),
				axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
				axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
				
				
			) + 
			ylab("CAPRA score") + xlab("Genes") +
			guides(
				color = guide_legend(keywidth = 0.5, keyheight = 0.5, direction = "vertical"),
				fill = guide_colorbar(barwidth = 5, barheight = 0.3, title.vjust = 1)
			) + 
			ggtitle( unique((.)$Category) %>% as.character())
		
		ggsave(
			sprintf("violin_%s.pdf", (.) %>% pull(Category) %>% unique %>% gsub(" |/", "_", .)),
			plot = p,
			useDingbats=FALSE,
			units = c("mm"),
			width = 183 ,
			height = 50 + ( 4 * (.) %>% distinct(symbol_ct, positively_associated) %>% group_by(positively_associated) %>% summarise(n = n()) %>% pull(n) %>% max )
		)
		
		tibble(Category = (.) %>% pull(Category) %>% unique, done = "yes")
		
	}) 

#############################################################
# Illustration plot cartoons ################################

summary_df %>% 
	
	# Add back membrane association
	left_join(
		tabi_res_E.annot %>% 
			bind_rows(tabi_res_F.annot) %>% 
			bind_rows(tabi_res_T.annot) %>% 
			bind_rows(tabi_res_M.annot) %>%
			filter(grepl("membrane", `Protein class`)) %>%
			distinct(symbol,  `Protein class`)
	) %>%
	mutate(is_membrane = !is.na(`Protein class`)) %>%
	
	# Correct automatic annotation - in the future use COMPARTMENTS
	mutate(is_membrane = ifelse(symbol %in% c("CSF1"), FALSE, is_membrane)) %>%
	
	# Eliminate for illustration sake the EMT > 2
	filter(abs(positively_associated) < 2) %>%	
	
	# Analyse each category
	group_by(Category) %>%
	do({
		browser()
		tibble(
			Category = (.) %>% distinct(Category) %>% pull(Category),
			plots = (.) %>% 
			{
				
				category_name = (.) %>% pull(Category) %>% unique %>% gsub(" |/", "_", .)
				
				# Separate in time fractions
				(.) %>%
					group_by(symbol, symbol_ct, `Effect size`, b1, positively_associated, `Cell type`, Category, is_membrane) %>%
					do(
						(.) %>% 
							summarise(
								`CAPRA 0-2 / Low grade` = (.) %>% filter(inflection <= 2) %>% nrow(),
								`CAPRA 3-5 / Medium grade` = (.) %>% filter(inflection > 2 & inflection <= 5) %>% nrow(),
								`CAPRA 6-8 / High grade` = (.) %>% filter(inflection > 5) %>% nrow()
							)
					) %>%
					ungroup() %>% 
					{
						
						# Save temp variable
						tbl_de = (.)
						
						# radius of the circle/square
						# http://www.wolframalpha.com/input/?i=solve+pi+*+r%5E2+%2F+4+%3D+x+for+r
						{
							ceiling(
								(
									2 * 
										sqrt( 
											
											# Max number of genes in this plot for non membrane proteins
											(.) %>% 
												filter(!is_membrane) %>%
												group_by(`Cell type`) %>% 
												summarise(n=n()) %>% 
												pull(n) %>% 
												max 
										) 
								) / 
									sqrt(pi)
							)
						} %>%
							
							# Do combinations
						{
							foreach(
								el = list( 
									list(`Cell type` = "E", x = (-(.):0) - 1, y = (0:(.)) + 1, radius = (.)),
									list(`Cell type` = "T", x = (0:(.)) + 1,  y = (0:(.)) + 1, radius = (.)),
									list(`Cell type` = "M", x = (0:(.)) + 1,  y = (-(.):0) - 1, radius = (.)),
									list(`Cell type` = "F", x = (-(.):0) - 1, y = (-(.):0) - 1, radius = (.))
								),
								.combine = bind_rows
							) %do% {
								set.seed(654321)
								
								
								bind_rows(
									# Add coordinates for secreted
									expand.grid(`Cell type` = el$`Cell type`, x = el$x, y = el$y, r = el$radius) %>%
										as_tibble() %>%
										filter(sqrt(x^2 + y^2) < (r +1)) %>%
										sample_frac(1) %>%
										mutate(gene_idx = 1:n()) %>%
										inner_join( 
											tbl_de %>% 
												filter(!is_membrane) %>%
												filter(`Cell type` == el$`Cell type`) %>%
												{
													# To avoid error if n() == 0
													if((.) %>% nrow > 0) (.) %>% mutate(gene_idx = 1:n())
													else (.) %>% mutate(gene_idx = 0)
												}
										) %>%
										dplyr::select(-gene_idx),
									
									# Add coordinates for membrane
									seq(
										0, 
										2*pi, 
										2*pi/
											(
												tbl_de %>%
													filter(is_membrane) %>%
													group_by(`Cell type`) %>% 
													summarise(n=n()) %>% 
													pull(n) %>% 
													max * 4 + 1
											)
									) %>% 
										
										# find coord for radiants
									{ foreach(t = (.), .combine = bind_rows) %do% { c(x =(el$radius + 2)*cos(t), y =(el$radius + 2)*sin(t)) } } %>%
										filter(x != 0 & y != 0) %>%
										# Filter right quadrant
										filter(x * el$x[1] > 0 & y * el$y[1] > 0 ) %>%
										sample_frac(1) %>%
										mutate(gene_idx = 1:n()) %>%
										inner_join( 
											tbl_de %>% 
												filter(is_membrane) %>%
												filter(`Cell type` == el$`Cell type`) %>%
												{
													# To avoid error if n() == 0
													if((.) %>% nrow > 0) (.) %>% mutate(gene_idx = 1:n())
													else (.) %>% mutate(gene_idx = 0)
												}
										) %>%
										dplyr::select(-gene_idx)
									
								)
							}
						} 
					}	%>%
					# Separate through disease progression and filter if density less 5%
					gather(when, count, c("CAPRA 0-2 / Low grade", "CAPRA 3-5 / Medium grade",  "CAPRA 6-8 / High grade")) %>%
					mutate(when = factor(when, levels = c("CAPRA 0-2 / Low grade", "CAPRA 3-5 / Medium grade",  "CAPRA 6-8 / High grade"))) %>%
					
					# Set alpha if already present in the past
					group_by(symbol_ct) %>%
					do(
						rbind(
							(.) %>% filter(count <= 1000*0.3) %>% mutate(alpha = "0"),
							(.) %>% filter(count > 1000*0.3) %>% mutate(alpha = c("1", rep("0.5", (.) %>% nrow() -1))) 
						)
					) %>%
					ungroup() %>%
					
					{
						
						r = (.) %>% filter(!is.na(r)) %>% pull(r) %>% unique 
						im <- png::readPNG('figure-07.png')
						img <- grid::rasterGrob(im, width=unit(1,"npc"), height=unit(1,"npc"))
						my_max = max(abs((.) %>% pull(`Effect size`)))
						
						list(
							(.)  %>% 
								ggplot(
									aes(
										x = x, 
										y = y, 
										label = symbol,
										alpha = alpha,
										frame = when
									)
								) + 
								
								# Annotate the plot
								annotation_custom(img) +
								geom_vline(xintercept = 0, linetype  = "dashed", alpha = 0.2, size=0.2) +
								geom_hline(yintercept = 0, linetype  = "dashed", alpha = 0.2, size=0.2) +
								annotate("point", x = -r-2, y = r+2, shape = 21, colour = "#c2e2d8", fill = "white", size = 7) +
								annotate("text", x = -r-2, y = r+2, label = "E", size = 4) +
								annotate("point", x = r+2, y = r+2, shape = 21, colour = "#f6c4dd", fill = "white", size = 7) +
								annotate("text", x = r+2, y = r+2, label = "T", size = 4) +
								annotate("point", x = -r-2, y = -r-2, shape = 21, colour = "#f4d3bc", fill = "white", size = 7) +
								annotate("text", x = -r-2, y = -r-2, label = "F", size = 4) +
								annotate("point", x = r+2, y = -r-2, shape = 21, colour = "#d9d7ea", fill = "white", size = 7) +
								annotate("text", x = r+2, y = -r-2, label = "M", size = 4) +
								
								# Plot data
								geom_point(aes(
									fill=`Effect size`, 
									size = sqrt(abs(`Effect size`)), 
									shape=is_membrane
								)) +
								ggrepel::geom_label_repel(
									size = 1.6, 
									point.padding = 0.3, 
									fontface = 'bold', 
									label.padding = 0.1, 
									label.size = 0,
									segment.size = 0.2,
									seed = 654321
								) +
								
								# Facet
								facet_grid(positively_associated ~ when) +
								
								# Scales
								scale_shape_manual(values = c("TRUE" = 22, "FALSE" = 21) , guide = FALSE) +
								xlim( - r -2, r + 2) +
								ylim( - r -2, r + 2) +
								scale_size(trans = signed_log, range = c(1,6), guide = FALSE) +
								scale_fill_distiller(
									palette = "Spectral",
									na.value = 'white',
									direction = 1,
									trans = signed_log,
									breaks = signed_log_breaks,
									labels = scales::scientific,
									limits=c(	-my_max,my_max)
								) +
								scale_alpha_manual(values = c("0" = 0, "0.5" = 1, "1" = 1) , guide=FALSE) +
								
								# Theme
								theme_bw() +
								theme(
									panel.border = element_blank(), 
									axis.line = element_line(),
									panel.grid.major = element_line(size = 0.2),
									panel.grid.minor = element_line(size = 0.1),
									#text = element_text(size=12),
									legend.position="bottom",
									aspect.ratio=1,
									axis.text =element_blank(),
									axis.ticks=element_blank(),
									axis.title = element_blank(),
									legend.text = element_text(size = 4, angle = 60, vjust = 0.5),
									legend.title = element_text(size = 4),
									strip.background = element_blank(),
									strip.text = element_text(size = 10)
								) +
								guides(	fill = guide_colorbar(barwidth = 5, barheight = 0.3, title.vjust = 1)	) +
								ggtitle( (.) %>% pull(Category) %>% unique )
							
							
							
						)}
			}
		)
		
	}) %>%
	{
		
		num_rows = (.) %>% nrow()
		
		# Integrate the three plots
		(.) %>%
			pull(plots) %>%
			#gridExtra::arrangeGrob( grobs=., ncol = 1 ) %>%
			cowplot::plot_grid(plotlist = ., align = "v", ncol = 1, axis="b", rel_widths = 1 ) %>%
			
			# Save plots
			ggsave(
				"temporal.pdf",
				plot = .,
				useDingbats=FALSE,
				units = c("mm"),
				width = 183 ,
				height = (122 + 30) * num_rows,
				limitsize = FALSE
			)
	}



#############################################################
# Illustration plot CI ######################################

summary_df %>% 
	mutate(symbol_ct = as.character(symbol_ct)) %>%
	
	# Uify the two categories for plotting purposes
	mutate(Category = ifelse(Category %in% c("EMT", "Tissue remodelling/migration"), "Epithelial cell migration", Category)) %>%
	group_by(Category) %>%
	
	# Plot
	do({
		
		df1 = (.) %>% 	
			distinct(symbol, `Cell type`, Category, `Effect size`, `b1.2.5%`, `b1.97.5%`, symbol_ct, log_y_cross) %>%
			mutate(`Safe estimate` = ifelse(`b1.2.5%` > 0, `b1.2.5%`, `b1.97.5%`)) %>%
			arrange(desc(`Safe estimate`)) %>%
			mutate(
				symbol_ct = 
					factor(
						symbol_ct, 
						levels = (.) %>% pull(symbol_ct)
					)
			) 
		
		my_max = df1 %>% pull(`Effect size`) %>% max()
		Category  = df1 %>% pull(Category) %>% unique()
		
		df2 = (.) %>%
			mutate(
				symbol_ct = 
					factor(
						symbol_ct, 
						levels = df1 %>% pull(symbol_ct) %>% levels
					)
			) 
		
		df1 %>%
		{
			
			list(
				
				# Slope
				(.) %>%
					ggplot(aes(x=symbol_ct, color = `Cell type`)) +
					geom_hline(yintercept = 0, linetype  = "dashed", alpha = 0.2, size=0.2) +
					geom_errorbar(aes(ymin=`b1.2.5%`, ymax=`b1.97.5%`), width = 0) +
					facet_grid(~`Cell type`, scales = "free_x", space = "free_x") +
					theme_bw() +
					theme(
						panel.border = element_blank(), 
						axis.line = element_line(),
						panel.grid.major = element_line(size = 0.2),
						panel.grid.minor = element_line(size = 0.1),
						text = element_text(size=12),
						legend.position="none",
						strip.background = element_blank(),
						axis.ticks.x = element_blank(),
						axis.text.x = element_blank(),
						axis.title.x  = element_blank(),
						axis.text.y  = element_text(size=8),
						axis.title.y  = element_text(margin = margin(t = 10, r = 5, b = 10, l = 10), size=10, angle=60, vjust = 0), 
						rect = element_blank()
					) +
					scale_color_brewer(palette="Dark2") +
					ylab("Slope") +
					ggtitle(Category),
				
				# Effect size
				(.) %>%
					mutate(`Baseline gene count` = exp(log_y_cross)) %>%
					gather(Label, Estimate, c("Baseline gene count", "Effect size")) %>%
					ggplot(aes(x=symbol_ct, y = Estimate, color = `Cell type`)) +
					geom_hline(yintercept = 1, linetype  = "dashed", alpha = 0.2, size=0.2) +
					geom_point(aes(shape = Label)) +
					facet_grid(~`Cell type`, scales = "free_x", space = "free_x") +
					theme_bw() +
					theme(
						panel.border = element_blank(), 
						axis.line = element_line(),
						panel.grid.major = element_line(size = 0.2),
						panel.grid.minor = element_line(size = 0.1),
						text = element_text(size=12),
						legend.position="none",
						strip.background = element_blank(),
						strip.text.x = element_blank(),
						axis.text.x = element_blank(),
						axis.title.x  = element_blank(),
						axis.ticks.x = element_blank(),
						axis.text.y  = element_text(size=8),
						axis.title.y  = element_text(margin = margin(t = 10, r = 5, b = 10, l = 10), size=10, angle=60, vjust = 0), 
						rect = element_blank()
					) +
					scale_color_brewer(palette="Dark2") +
					scale_shape_manual(values = c("Effect size" = 16, "Baseline gene count" = 1)) +
					scale_y_continuous(	
						breaks = signed_log_breaks,
						labels = scales::scientific,
						limits=c(	-my_max,my_max),
						trans = signed_log
					) +
					ylab("Magnitude"),
				
				# Iflection
				df2 %>%
					ggplot(aes(x = symbol_ct, y = inflection, fill = `Cell type`)) +
					geom_hline(yintercept = 0, linetype  = "dashed", alpha = 0.2, size=0.2) +
					geom_hline(yintercept = 7, linetype  = "dashed", alpha = 0.2, size=0.2) +
					geom_boxplot(outlier.shape=NA, lwd=0.1) +
					scale_fill_brewer(palette="Dark2") +
					facet_grid(~`Cell type`, scales = "free_x", space = "free_x") +
					theme_bw() +
					theme(
						panel.border = element_blank(), 
						axis.line = element_line(),
						panel.grid.major = element_line(size = 0.2),
						panel.grid.minor = element_line(size = 0.1),
						text = element_text(size=12),
						legend.position="none",
						strip.background = element_blank(),
						strip.text.x = element_blank(),
						axis.text.y  = element_text(size=8),
						axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=6),
						axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10), size=10),
						axis.title.y  = element_text(margin = margin(t = 10, r = 5, b = 10, l = 10), size=10, angle=60, vjust = 0), 
						rect = element_blank()
					) +
					ylim(c(-15, 15)) +
					ylab("Inflection") +
					xlab("Genes") 
				
			) %>%
				
				# manage grid and legens
			{
				cowplot::plot_grid(
					plotlist = (.), 
					align = "v", 
					ncol = 1,
					axis="b", rel_widths = 1, rel_heights = c(0.35, 0.18, 0.47)
				) 
			} %>%
				
				# Save plots
				ggsave(
					sprintf("slope_effectSize_plot_%s.pdf", Category %>% gsub("[ \\/]", "_", .)),
					plot = .,
					useDingbats=FALSE,
					units = c("mm"),
					width = 183 ,
					height = 90,
					limitsize = FALSE
				)
			
			tibble(	Category = Category	)
			
		}
	})

#############################################################
# Animatiton plot ###########################################

summary_df %>% 
	
	# Add back membrane association
	left_join(
		tabi_res_E.annot %>% 
			bind_rows(tabi_res_F.annot) %>% 
			bind_rows(tabi_res_T.annot) %>% 
			bind_rows(tabi_res_M.annot) %>%
			filter(grepl("membrane", `Protein class`)) %>%
			distinct(symbol,  `Protein class`)
	) %>%
	mutate(is_membrane = !is.na(`Protein class`)) %>%
	
	# Correct automatic annotation - in the future use COMPARTMENTS
	mutate(is_membrane = ifelse(symbol %in% c("CSF1"), FALSE, is_membrane)) %>%
	
	# Analyse each category
	filter(Category == "Angiogenesis") %>%
	
	# Normalise inflection density
	filter(inflection >= 0 & inflection <= 7) %>%
	mutate(gr=cut(inflection, breaks= seq(0, 7, by = 0.5)) ) %>%
	mutate(gr = factor(gr)) %>%
	droplevels() %>%
	do({
		
		tt = 	# Create table of smoothed density
			table((.)$symbol_ct, (.)$gr) %>% 
			data.frame() %>%
			as_tibble() %>%
			setNames(c("symbol_ct", "gr", "count")) %>%
			group_by(symbol_ct) %>%
			mutate(count_norm = count/max(count)) %>%
			ungroup() %>%
			
			# Add numerical CAPRA
			left_join(
				(.) %>% distinct(gr) %>% mutate(gr_idx = approx(c(0,7), n = n())$y)
			) %>%
			
			# Smooth
			group_by(symbol_ct) %>%
			do({
				symbol_ct = (.) %>% pull(symbol_ct) %>% unique %>% as.character
				bezierCurve(
					(.) %>% pull(gr_idx), 
					(.) %>% pull(count_norm),
					100
				) %>%
					as_tibble() %>%
					setNames(c("gr_idx", "count_norm")) %>%
					mutate(symbol_ct = symbol_ct)
			})
		
		(.) %>% 
			distinct(symbol, symbol_ct, `Effect size`, b1, positively_associated, `Cell type`, Category, is_membrane) %>%
			right_join(tt) 
		
	}) %>%
	
	do({
		
		category_name = (.) %>% pull(Category) %>% unique %>% gsub(" |/", "_", .)
		
		# Separate in time fractions
		(.) %>%
		{
			# Save temp variable
			tbl_de = (.)
			
			# radius of the circle/square
			# http://www.wolframalpha.com/input/?i=solve+pi+*+r%5E2+%2F+4+%3D+x+for+r
			{
				ceiling(
					(
						2 * 
							sqrt( 
								
								# Max number of genes in this plot for non membrane proteins
								(.) %>% 
									filter(!is_membrane) %>%
									distinct(symbol, `Cell type`) %>%
									group_by(`Cell type`) %>% 
									summarise(n=n()) %>% 
									pull(n) %>% 
									max 
							) 
					) / 
						sqrt(pi)
				)
			} %>%
				
				# Do combinations
			{
				foreach(
					el = list( 
						list(`Cell type` = "E", x = (-(.):0) - 1, y = (0:(.)) + 1, radius = (.)),
						list(`Cell type` = "T", x = (0:(.)) + 1,  y = (0:(.)) + 1, radius = (.)),
						list(`Cell type` = "M", x = (0:(.)) + 1,  y = (-(.):0) - 1, radius = (.)),
						list(`Cell type` = "F", x = (-(.):0) - 1, y = (-(.):0) - 1, radius = (.))
					),
					.combine = bind_rows
				) %do% {
					set.seed(654321)
					
					
					bind_rows(
						# Add coordinates for secreted
						expand.grid(`Cell type` = el$`Cell type`, x = el$x, y = el$y, r = el$radius) %>%
							as_tibble() %>%
							filter(sqrt(x^2 + y^2) < (r +1)) %>%
							sample_frac(1) %>%
							mutate(gene_idx = 1:n()) %>%
							right_join( 
								tbl_de %>% 
									filter(!is_membrane) %>%
									filter(`Cell type` == el$`Cell type`) %>%
									left_join(
										(.) %>% distinct(symbol) %>% mutate(gene_idx = 1:n())
									)
							) %>%
							dplyr::select(-gene_idx),
						
						# Add coordinates for membrane
						seq(
							0, 
							2*pi, 
							2*pi/
								(
									tbl_de %>%
										filter(is_membrane) %>%
										group_by(`Cell type`) %>% 
										summarise(n=n()) %>% 
										pull(n) %>% 
										max * 4 + 1
								)
						) %>% 
							
							# find coord for radiants
						{ 
							
							foreach(t = (.), .combine = bind_rows) %do% { c(x =(el$radius + 2)*cos(t), y =(el$radius + 2)*sin(t)) } } %>%
							filter(x != 0 & y != 0) %>%
							# Filter right quadrant
							filter(x * el$x[1] > 0 & y * el$y[1] > 0 ) %>%
							sample_frac(1) %>%
							mutate(gene_idx = 1:n()) %>%
							right_join( 
								tbl_de %>% 
									filter(is_membrane) %>%
									filter(`Cell type` == el$`Cell type`) %>%
									left_join(
										(.) %>% distinct(symbol) %>% mutate(gene_idx = 1:n())
									)
							) %>%
							dplyr::select(-gene_idx)
						
					)
					
					
				}
			} 
		}	%>%
			
			# Plot
		{
			
			r = (.) %>% filter(!is.na(r)) %>% pull(r) %>% unique 
			im <- png::readPNG('figure-07.png')
			img <- grid::rasterGrob(im, width=unit(1,"npc"), height=unit(1,"npc"))
			my_max = max(abs((.) %>% pull(`Effect size`)))
			
			(.) %>%
				group_by(gr_idx) %>%
				do({
					gr_idx = (.) %>% pull(gr_idx) %>% unique
					
					{
						ggplot(
							(.),
							aes(
								x = x, 
								y = y, 
								label = symbol,
								alpha = count_norm
							)
						) + 
							
							# Annotate the plot
							annotation_custom(img) +
							geom_vline(xintercept = 0, linetype  = "dashed", alpha = 0.2, size=0.2) +
							geom_hline(yintercept = 0, linetype  = "dashed", alpha = 0.2, size=0.2) +
							annotate("point", x = -r-2, y = r+2, shape = 21, colour = "#c2e2d8", fill = "white", size = 7) +
							annotate("text", x = -r-2, y = r+2, label = "E", size = 4) +
							annotate("point", x = r+2, y = r+2, shape = 21, colour = "#f6c4dd", fill = "white", size = 7) +
							annotate("text", x = r+2, y = r+2, label = "T", size = 4) +
							annotate("point", x = -r-2, y = -r-2, shape = 21, colour = "#f4d3bc", fill = "white", size = 7) +
							annotate("text", x = -r-2, y = -r-2, label = "F", size = 4) +
							annotate("point", x = r+2, y = -r-2, shape = 21, colour = "#d9d7ea", fill = "white", size = 7) +
							annotate("text", x = r+2, y = -r-2, label = "M", size = 4) +
							
							# Plot data
							geom_point(aes(
								fill=`Effect size`, 
								size = count_norm, #  sqrt(abs(`Effect size`)), 
								shape=is_membrane
							)) +
							ggrepel::geom_label_repel(
								size = 1.6, 
								point.padding = 0.3, 
								fontface = 'bold', 
								label.padding = 0.1, 
								label.size = 0,
								segment.size = 0.2,
								seed = 654321
							) +
							
							# Scales
							scale_shape_manual(values = c("TRUE" = 22, "FALSE" = 21) , guide = FALSE) +
							xlim( - r -2, r + 2) +
							ylim( - r -2, r + 2) +
							scale_size(range = c(0,6), guide = FALSE) +
							scale_fill_distiller(
								palette = "Spectral",
								na.value = 'white',
								direction = 1,
								trans = signed_log,
								breaks = signed_log_breaks,
								labels = scales::scientific,
								limits=c(	-my_max,my_max)
							) +
							scale_alpha_continuous(range = c(0, 1), guide = FALSE, trans = cusotm_root_trans()) +
							#scale_alpha_manual(values = c("0" = 0, "1" = 1) , guide=FALSE ) +
							
							# Theme
							theme_bw() +
							theme(
								panel.border = element_blank(), 
								axis.line = element_line(),
								panel.grid.major = element_line(size = 0.2),
								panel.grid.minor = element_line(size = 0.1),
								text = element_text(size=12),
								legend.position="bottom",
								aspect.ratio=1,
								axis.text =element_blank(),
								axis.ticks=element_blank(),
								axis.title = element_blank(),
								legend.text = element_text(angle = 60, vjust = 0.5),
								axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
								axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
								
								
								#,
								#legend.position="right"
							) +
							#guides(	fill = guide_colorbar(barwidth = 5, barheight = 0.3, title.vjust = 1)	) +
							ggtitle(
								sprintf(
									"%s - CAPRA risk score %s",  
									(.) %>% pull(Category) %>% unique,
									round(gr_idx)
								)
							) 
					} %>% 
						# Save plots
						ggsave(
							sprintf("temporal_animation_%s%s.png", ifelse(gr_idx<10, "0", ""), format(gr_idx, nsmall = 20)),
							plot = .,
							width = 4 ,
							height = 4.5,
							limitsize = FALSE
						)
					
					tibble(gr_idx = gr_idx, done = "yes")
				}
				)}
	}) 





#############################################################
# Descriptive statistics ####################################

# DE
foreach(
	ct_list = 
		list( 
			E = list(tabi_res_E.annot, tabi_gen_E),
			F = list(tabi_res_F.annot, tabi_gen_F),
			T = list(tabi_res_T.annot, tabi_gen_T),
			M = list(tabi_res_M.annot, tabi_gen_M)
		),
	.combine = bind_rows,
	ct = c("E", "F", "T", "M")
) %do% {
	
	ct_list[[1]] %>%
		filter_too_sparse(ct_list[[2]]) %>%
		filter_out_CI(ct_list[[2]]) %>%
		mutate(is_interfaced = grepl("secreted|membrane", `Protein class`)) %>%
		distinct(symbol, is_interfaced) %>%
		group_by(is_interfaced) %>%
		summarise(n = n()) %>%
		mutate(	prop = 	n/(.) %>% pull(n) %>% sum()) %>%
		mutate(`Cell type` = ct) 
} %>%
{
	print((.))
	
	(.) %>% 
		filter(is_interfaced == T) %>% 
		summarise(m = mean(prop))
}

# Cancer genes
summary_df %>% 
	filter(Category %in% c("Tumor gene", "Prostate tumor gene")) %>% 
	distinct(symbol, `Cell type`, positively_associated) %>%
	count(`Cell type`, positively_associated) %>% 
	spread(positively_associated, n) %>%
	mutate(tot = `1`+`-1`) %>%
	mutate(consistent = `1`/(`1`+`-1`))

# PC genes
summary_df %>% 
	filter(Category %in% c("Prostate tumor gene")) %>%
	distinct(symbol, `Cell type`, positively_associated) %>%
	count(`Cell type`, positively_associated) %>% 
	spread(positively_associated, n) %>%
	mutate(tot = `1`+`-1`) %>%
	mutate(consistent = `1`/(`1`+`-1`))

# Enrichment between categories
summary_df %>%
	dplyr::distinct(symbol, `Cell type`, Category, `positively_associated`) %>%
	mutate_if(is.character, as.factor) %>%
	group_by(Category, positively_associated) %>%
	count() %>%
	spread(positively_associated ,n) %>%
	mutate(
		`P-value` = prop.test(
			table(
				c ( 
					rep("1", `1`),
					rep("-1", `-1`) 
				)
			), 
			alternative = "greater"
		)$p.value %>% 
			scientific() %>%
			as.numeric()
	) %>%
	ungroup() %>%
	mutate("Adj. p-value" = `P-value` * n()) %>%
	mutate(`Adj. p-value` = ifelse(`Adj. p-value`>1, 1, `Adj. p-value`)) %>%
	mutate(sig = ifelse(`Adj. p-value`<0.05, "*", "")) 

# Cancer genes for each cell type
summary_df %>%
	dplyr::distinct(symbol, `Cell type`, Category, `positively_associated`) %>%
	filter(Category == "Tumor gene") %>%
	mutate_if(is.character, as.factor) %>%
	group_by(`Cell type`, positively_associated) %>%
	count() %>%
	spread(positively_associated ,n) %>%
	mutate(
		`P-value` = prop.test(
			table(
				c ( 
					rep("1", `1`),
					rep("-1", `-1`) 
				)
			), 
			alternative = "greater"
		)$p.value %>% 
			scientific() %>%
			as.numeric()
	) %>%
	ungroup() %>%
	mutate("Adj. p-value" = `P-value` * n()) %>%
	mutate(`Adj. p-value` = ifelse(`Adj. p-value`>1, 1, `Adj. p-value`)) %>%
	mutate(sig = ifelse(`Adj. p-value`<0.05, "*", "")) 

# Antiinflammation in higher for edvanced stages of the disease
summary_df %>%
	filter(Category == "Immune modulator") %>%
	group_by(symbol, `Cell type`, positively_associated) %>%
	summarise(inflection_mean = median(inflection)) %>%
	filter(positively_associated == 1) %>%
	group_by(late = inflection_mean > 2) %>%
	count() %>%
	spread(late ,n) %>%
	mutate(
		`P-value` = prop.test(
			table(
				c ( 
					rep("TRUE", `TRUE`),
					rep("FALSE", `FALSE`) 
				)
			), 
			alternative = "less"
		)$p.value %>% 
			scientific() %>%
			as.numeric()
	)

# Angiogenic late enrichment
summary_df %>%
	filter(Category == "Angiogenesis") %>%
	group_by(symbol, `Cell type`, positively_associated) %>%
	summarise(inflection_mean = median(inflection)) %>%
	group_by(positively_associated, late = inflection_mean > 2) %>%
	count() %>%
	spread(positively_associated ,n) %>%
	mutate(
		`P-value` = prop.test(
			table(
				c ( 
					rep("1", `1`),
					rep("-1", `-1`) 
				)
			), 
			alternative = "greater"
		)$p.value %>% 
			scientific() %>%
			as.numeric()
	) %>%
	ungroup() %>%
	mutate("Adj. p-value" = `P-value` * n()) %>%
	mutate(`Adj. p-value` = ifelse(`Adj. p-value`>1, 1, `Adj. p-value`)) %>%
	mutate(sig = ifelse(`Adj. p-value`<0.05, "*", ""))

#############################################################
# Signatures ################################################

d_edgeR = d
d_edgeR$samples$norm.factors = d_edgeR$counts %>% 
	as_tibble(rownames="gene_idx") %>%
	gather(sample, value, -gene_idx) %>%
	norm_RNAseq( 
		sample_column = "sample", 
		gene_column = "gene_idx", 
		value_column = "value"
	) %>% 
	filter(!filt_for_calc) %>%
	distinct(sample, TMM) %>%
	arrange(match(colnames(d_edgeR), sample)) %>%
	pull(TMM)

foreach(ct = c("E", "F", "M", "T"), .combine = bind_rows) %do% {
	design <- 	model.matrix(
		~
			annot %>% filter(cell_type_formatted == ct) %>% pull(CAPRA_TOTAL) %>% `>`(2) + 
			annot %>% filter(cell_type_formatted == ct) %>% pull(batch) 
	)
	colnames(design) = c("(Intercept)", "high_grade", "batch")
	
	d_edgeR_E = d_edgeR[,which(annot$cell_type_formatted == ct)]
	# EGSEA
	# If file does NOT exist
	d_edgeR_E %>%
		voom(design, plot=FALSE) %>%
		
		# Run gene enrichment
		{
			v = (.)
			colnames(v$genes) = c("ENTREZID", "length", "SYMBOL")
			v$genes = v$genes[,c(1,3,2)]
			
			library(EGSEA)
			library(EGSEAdata)
			
			idx = buildIdx(entrezIDs=rownames(v), species="human")
			
			v %>%
				egsea(
					contrasts=2, 
					gs.annots=idx, 
					symbolsMap=
						v %$% 
						genes %>% 
						dplyr::select(1:2) %>%
						setNames(c("FeatureID", "Symbols")),
					baseGSEAs = egsea.base()[-c(6, 7, 8, 9)],
					sort.by="med.rank",
					num.threads = 10
				)
		} %>%
		{
			egsea.results = (.)
			save(egsea.results, file=sprintf("EGSEA_%s_results.RData", ct))
			(.)
		}
}


library(GSEABase)
library(GSVA)

# EMT ########################################################
# https://www.ncbi.nlm.nih.gov/pubmed/29921693


d_edgeR_E <- estimateDisp(d_edgeR_E,design)
fit <- glmQLFit(d_edgeR_E,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf, n = 999999) %>% 
	as.data.frame() %>%
	as_tibble() %>% 
	#filter(PValue < 0.05) %>%
	arrange(desc(abs(logFC))) %>%
	dplyr::select(symbol, logFC) %>%
	drop_na() %>%
	write_delim("DE_TME_GSEA.rnk", col_names = F, delim = "\t")

system(sprintf("java -cp gsea2-2.2.2.jar -Xmx8g xtools.gsea.GseaPreranked -gmx %s -rnk %s -collapse false -mode Max_probe -norm meandiv -nperm 5000 -scoring_scheme weighted -include_only_symbols true -make_sets true -plot_top_x 100 -rnd_seed timestamp -set_max 500 -set_min 1 -zip_report false -out gsea -gui false", "signature_EMT_3752.gmt", "DE_TME_GSEA.rnk"))

#############################################################
# Plot gene category ########################################

tabi_gen_E %>% 
	
	# Select category
	left_join(
		tabi_res_E.annot %>% 
			filter(namespace_1003 == "biological_process") %>% 
			filter(name_1006 == "cell adhesion") %>% 
			distinct(symbol, `Gene description`) %>%
			filter(grepl("Protocadherin", `Gene description` )) %>%
			mutate(plot=T) %>%
			dplyr::rename(gene=symbol)
	) %>%
	filter(plot) %>% 
	distinct() %>%
	
	# Plot
	ggplot(aes(x=CAPRA_TOTAL, y=`read count`)) +
	geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width = 0) +
	geom_jitter(color="red") +
	facet_wrap(~gene, scale="free") +
	scale_y_log10() +
	theme_bw() + theme(strip.background = element_blank()) + ggtitle("Protocadherin")

#############################################################
# Check bacteria viruses ####################################

stats_metgen = do.call("rbind", lapply(dir(path = "alignment_hg38", recursive=T, pattern="Unmapped.out.fixed.kraken.report", full.names = T), function(i){
	df <- suppressMessages(suppressWarnings(
		read_delim(i,  "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
	))
	colnames(df) = c("proportion_of_unmapped", "count_clade", "count_taxon", "rank_code", "taxa_id", "taxa_name")
	df$sample = strsplit(i, "/|_")[[1]][5]
	df$`N adapter` = strsplit(i, "_|-")[[1]][4]
	df$`S adapter` = strsplit(i, "_|-")[[1]][5]
	df$reference_genome = strsplit(i, "/")[[1]][1]
	df$belongs_to = df$belongs_to_count = NA
	
	# annotate with Phylum
	phila = (df %>% filter(rank_code=="D"))$taxa_name 
	positions = which(df$taxa_name%in%phila)
	for(i in 1:(length(phila)-1)){
		df$belongs_to[positions[i]:(positions[i+1]-1)] = phila[i]
		df$belongs_to_count[positions[i]:(positions[i+1]-1)] = (df %>% filter(taxa_name==phila[i]))$count_clade
	}
	
	df
})) %>% 
	filter(sample!="S1") %>% 
	left_join(
		annot %>% dplyr::select(`N adapter`, `S adapter`, CAPRA_TOTAL, cell_type_formatted, `Number of input reads`
		),
		by = c("N adapter", "S adapter")
	) %>%
	mutate(proportion_of_all = count_clade / `Number of input reads`)%>%
	mutate(proportion_of_domain = count_clade/belongs_to_count)


library(ggpmisc)

p1 = ggplot(
	stats_metgen %>% filter(rank_code=="D") %>% dplyr::select(cell_type_formatted, proportion_of_all, CAPRA_TOTAL, taxa_name), 
	aes(x = CAPRA_TOTAL, y = proportion_of_all, fill=cell_type_formatted)) +
	geom_smooth(colour="grey20", size=0.5, method="glm", method.args = list(family =quasibinomial()), se=F) +
	geom_jitter(size=2, shape=21, color="grey20") + 
	facet_grid(cell_type_formatted~taxa_name) +
	scale_fill_brewer(palette="Dark2") +
	theme_bw() +
	theme(
		panel.border = element_blank(), 
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size=12),
		legend.position="bottom",
		axis.text.x = element_text(angle = 60, hjust = 1), 
		strip.text.y = element_text(angle = 0), 
		strip.text.x = element_text(angle = 90, vjust = 0),
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
		
		
	) +
	xlab("CAPRA score") + ylab("Proportion of sequenced reads")

p2 = ggplot(
	stats_metgen %>% filter(belongs_to =="Bacteria" & rank_code=="P") %>% filter(proportion_of_domain>0.01), aes(x = CAPRA_TOTAL, y = proportion_of_domain, fill=cell_type_formatted)) +
	geom_smooth(colour="grey20", size=0.5, method="glm", method.args = list(family =quasibinomial()), se=F) +
	geom_jitter(size=2, shape=21, color="grey20") + 
	facet_grid(cell_type_formatted~taxa_name) +
	scale_fill_brewer(palette="Dark2") +
	theme_bw() +
	theme(
		panel.border = element_blank(), 
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size=12),
		legend.position="bottom",
		axis.text.x = element_text(angle = 60, hjust = 1), 
		strip.text.y = element_text(angle = 0), 
		strip.text.x = element_text(angle = 90, vjust = 0),
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
		
		
	) +
	xlab("CAPRA score") + ylab("Proportion of bacteria phila")

p3 = ggplot(
	stats_metgen %>% filter(belongs_to =="Bacteria" & rank_code=="C") %>% filter(proportion_of_domain>0.01)  %>% group_by(cell_type_formatted, taxa_name) %>% filter(n() > 5) %>% ungroup(), aes(x = CAPRA_TOTAL, y = proportion_of_domain, fill=cell_type_formatted)) +
	geom_smooth(colour="grey20", size=0.5, method="glm", method.args = list(family =quasibinomial()), se=F) +
	geom_jitter(shape=21, color="grey20") + 
	facet_grid(cell_type_formatted~taxa_name) +
	scale_fill_brewer(palette="Dark2") +
	theme_bw() +
	theme(
		panel.border = element_blank(), 
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size=12),
		legend.position="bottom",
		axis.text.x = element_text(angle = 60, hjust = 1), 
		strip.text.y = element_text(angle = 0), 
		strip.text.x = element_text(angle = 90, vjust = 0),
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
		
		
	) +
	xlab("CAPRA score") + ylab("Proportion of bacteria classes")

p4 = ggplot(
	stats_metgen %>% filter(belongs_to =="Bacteria" & rank_code=="O")%>% filter(proportion_of_domain>0.01) %>% group_by(cell_type_formatted, taxa_name) %>% filter(n() > 5) %>% ungroup()	, aes(x = CAPRA_TOTAL, y = proportion_of_domain, fill=cell_type_formatted)) +
	geom_smooth(colour="grey20", size=0.5, method="glm", method.args = list(family =quasibinomial()), se=F) +
	geom_jitter( shape=21, color="grey20") + 
	facet_grid(cell_type_formatted~taxa_name) +
	scale_fill_brewer(palette="Dark2") +
	theme_bw() +
	theme(
		panel.border = element_blank(), 
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size=12),
		legend.position="bottom",
		axis.text.x = element_text(angle = 60, hjust = 1), 
		strip.text.y = element_text(angle = 0), 
		strip.text.x = element_text(angle = 90, vjust = 0),
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
		
		
	) +
	xlab("CAPRA score") + ylab("Proportion of bacteria orders")

ggsave("microorganism_component.pdf", useDingbats=FALSE, p1 + theme(legend.position="bottom"), units = c("mm"), width = 183 , height = 140)


require(cowplot)
pcol = plot_grid(
	p2 + theme(legend.position="none"),
	p3 + theme(legend.position="none"),
	p4 + theme(legend.position="bottom"),
	labels = "AUTO",
	ncol = 1, align = 'v', axis = 'l'
)

ggsave("bacterial_component.pdf", useDingbats=FALSE, pcol, units = c("mm"), width = 183 , height = 183*2)

annot = annot %>% left_join(
	stats_metgen %>% 
		filter(rank_code=="D") %>% 
		dplyr::select(`N adapter`, `S adapter`, cell_type_formatted, proportion_of_all, taxa_name) %>%
		filter(taxa_name == "Bacteria") %>%
		rename(`Bacterial reads proportion` = proportion_of_all) %>%
		select(-taxa_name)
)

make_pair_plot_features()


#############################################################
# Validation ################################################

tibble(
	symbol = 
		c(
			# Mimicry
			"FCER1G",
			"FCGR1A",
			"FCER1G",
			"HLA-DRB5",
			"CSF1",
			"CSF1",
			
			# S1P
			"SPP2",
			"SPNS2",
			
			# Fibro-Mac
			"CXCL10",
			"CXCL14",
			"SLAMF1",
			"COL1A2",
			"CYR61",
			
			"SPINK5",
			
			# Hormone
			"CYP51A1",
			"STARD3NL",
			"SC5D",
			"CES3",
			"CH25H",
			"AIG1"
		),
	"Cell type" = c("E", "E", "E", "E", "M", "T", "M", "E", "F", "F", "F", "M", "M", "E", "M", "M", "T", "E", "T", "M")
)


tabi_res_E.annot %>%
	filter(
		symbol %in%
			c(
				# Mimicry
				"FCER1G",
				"FCGR1A",
				"FCER1G",
				"HLA-DRB5",
				"SPNS2",
				"SPINK5",
				
				# Hormone
				"CES3"
			)
	) %>%
	print_gene_regression(
		tabi_gen_E %>%
			ungroup() %>%
			filter(gene %in% (tabi_res_E.annot %>% pull(symbol) %>% unique()) ),
		tabi_inflection_E
	)


tabi_res_T.annot %>%
	filter(
		symbol %in%
			c(
				"CSF-1",
				"SC5D",
				"CH25H"	
			)
	) %>%
	print_gene_regression(
		tabi_gen_T %>%
			ungroup() %>%
			filter(gene %in% (tabi_res_T.annot %>% pull(symbol) %>% unique()) ),
		tabi_inflection_T
	)


tabi_res_M.annot %>%
	filter(
		symbol %in%
			c(
				"CSF-1", "SPP2",	"COL1A2",
				"CYR61",			"CYP51A1",
				"STARD3NL",			"AIG1"
			)
	) %>%
	print_gene_regression(
		tabi_gen_M %>%
			ungroup() %>%
			filter(gene %in% (tabi_res_M.annot %>% pull(symbol) %>% unique()) ),
		tabi_inflection_M
	)

tabi_res_F.annot %>%
	filter(
		symbol %in%
			c(
				"CXCL10",
				"CXCL14",		
				"SLAMF1"
			)
	) %>%
	print_gene_info(
		tabi_gen_F %>%
			ungroup() %>%
			filter(gene %in% (tabi_res_F.annot %>% pull(symbol) %>% unique()) ),
		tabi_inflection_F
	)


#############################################################
# Osteo signature ###########################################

foreach(
	ct_list = 
		list( 
			E = list(annot = tabi_res_E.annot, gen = tabi_gen_E),
			F = list(annot = tabi_res_F.annot, gen = tabi_gen_F),
			T = list(annot = tabi_res_T.annot, gen = tabi_gen_T),
			M = list(annot = tabi_res_M.annot, gen = tabi_gen_M)
		),
	.combine = bind_rows,
	ct = c("E", "F", "T", "M")
) %do% 
{
	ct_list %$% annot %>%
		left_join(
			read_sheet(ss = my_url, sheet = ct) %>%
				#filter(!is.na(`Osteogenesis`)) %>%
				rename(symbol = gene)
		) %>%
		filter_too_sparse(ct_list %$% gen, prop = ifelse(ct == "M", 13/4, 13/3)) %>%
		filter_out_CI(ct_list %$% gen, prop = ifelse(ct == "M", 13/4, 13/3)) 	%>%
		mutate(`Cell type` = ct)
}  %>%
	
	# merge curated and non curated
{
	
	bind_rows(
		
		# intracellular automatically annotated
		(.) %>%
			filter(grepl("bone|osteo|cartilage|chondro", name_1006)) %>%
			filter(!symbol %in% ( (.) %>% filter(grepl("secreted|membrane", `Protein class`)) %>% pull(symbol))) %>%
			mutate(Curated = "no"),
		
		# extracellular manually curated
		(.) %>%
			filter(!is.na(`Osteogenesis`)) %>%
			mutate(Curated = "yes")
		
	)
} %>%
	
	# Print violin
{
	(.) %>%
		distinct(`Cell type`, symbol, name_1006, description, Curated, `Cell type`) %>% 
		write_csv("tabi_osteogenesis.csv")
	
	my_tbl =		(.) %>%
		distinct(
			symbol, `Cell type`,
			#log_y_cross,
			`log_y_cross.2.5%`, `log_y_cross.97.5%`, 
			b0, `b0.2.5%`, `b0.97.5%`, 
			b1, `b1.2.5%`, `b1.97.5%`
		) %>%
		# Add inflection data sets
		inner_join(
			
			foreach(
				ct_list = 
					list( 
						E = tabi_inflection_E,
						F = tabi_inflection_F,
						T = tabi_inflection_T,
						M = tabi_inflection_M
					),
				.combine = bind_rows,
				ct = c("E", "F", "T", "M")
			) %do% {
				ct_list %>% 
					mutate(`Cell type` = ct) %>%
					rename(symbol = gene)
			}
		) %>%
		
		# Calculate Effect size
		rowwise() %>%
		mutate(`Effect size` = effect_size(b0, b1, log_y_cross)) %>%
		ungroup() %>%
		
		# Flip effect size if inhibitor
		#mutate(`Effect size` = `Effect size` * value)	%>%
		
		# Correct inflection
		rowwise() %>%
		mutate(inflection = get_log_inflection(log_y_cross, inflection, b1)) %>%	 
		mutate(inflection = rescale_inflection(inflection)) %>%	 
		ungroup() %>%
		filter(inflection != "NaN") %>%
		
		# Factorise
		mutate(`Cell type` = factor(`Cell type`)) %>%
		
		# Distinguish same symbol of different cell types
		mutate(symbol_ct = interaction(symbol, `Cell type`) ) %>%
		
		# Arrange on inflection
		do({
			my_order = (.) %>% 
				distinct(symbol_ct, inflection) %>%
				group_by(symbol_ct) %>%
				summarise(mean_inflection = median(inflection)) %>%
				arrange(mean_inflection) %>% 
				pull(symbol_ct)
			(.) %>% mutate(	symbol_ct = factor(symbol_ct, levels = my_order))
		}) %>%
		
		# See if pro or anti tumor
		mutate(positively_associated = `Effect size` / abs(`Effect size`) ) %>%
		mutate(positively_associated = ifelse(abs(positively_associated)>1, 2, positively_associated)) %>%
		mutate(positively_associated = 1) %>%
		
		# mutate effect size because buggy
		group_by(symbol_ct) %>%
		mutate(`Effect size` = mean(`Effect size`)) %>%
		ungroup()
	
	# Plot
	(
		ggplot(my_tbl, aes(x =symbol_ct , y=inflection, fill=`Effect size`)) +
			geom_violin(
				scale = "width", adjust = 1, width = 1
				
			) +
			scale_fill_distiller(
				palette = "Spectral",
				na.value = 'white',
				direction = 1,
				trans = signed_log,
				breaks = signed_log_breaks,
				labels = scales::scientific,
				limits=c(
					-max(abs(my_tbl %>% pull(`Effect size`))),
					max(abs(my_tbl %>% pull(`Effect size`)))
				)
			) +
			geom_point( 
				data = my_tbl %>% 
					group_by(symbol_ct, positively_associated) %>% 
					summarise(
						i = quantile(inflection, 0.5), 
						`Effect size` = mean(`Effect size`, na.rm = T)
					),
				aes(y = i, x = symbol_ct),
				size = 0.5
			) +
			geom_boxplot(
				data = my_tbl %>% 
					distinct(
						symbol_ct, 
						`Effect size`,
						`Cell type`, 
						inflection, 
						positively_associated
					) %>% 
					mutate(lower = -20, upper = -19),
				aes(x=symbol_ct, ymin = lower, lower = lower, middle = lower, upper = upper, ymax = upper, color=`Cell type`),
				stat = "identity"
			) +
			scale_color_brewer(palette="Dark2") +
			coord_flip(
				ylim=c(-25,17),
				xlim=c(1, my_tbl %>% distinct(symbol_ct, positively_associated) %>% group_by(positively_associated) %>% summarise(n = n()) %>% pull(n) %>% max )
			)+
			facet_wrap(
				~ positively_associated, 
				scales = "free", 
				#space = "free",
				drop = T
			) +
			scale_y_continuous(
				trans = magnify_trans(interval_low = 0, interval_high = 7,  reducer = 20),
				#trans = magnify_trans(intercept =  7,  reducer = 20),
				breaks = c(-10, 0, 1:7, 17)
			) +
			theme_bw() +
			theme(
				panel.border = element_blank(), 
				axis.line = element_line(),
				panel.grid.major = element_line(size = 0.2),
				panel.grid.minor = element_line(size = 0.1),
				text = element_text(size=12),
				legend.position="bottom",
				legend.text = element_text(angle = 60, vjust = 0.5),
				axis.text.x = element_text(angle = 60, vjust = 0.5),
				strip.background = element_blank(),
				axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
				axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
				
				
			) + 
			ylab("CAPRA score") + xlab("Genes") +
			guides(
				color = guide_legend(keywidth = 0.5, keyheight = 0.5, direction = "vertical"),
				fill = guide_colorbar(barwidth = 5, barheight = 0.3, title.vjust = 1)
			) + 
			ggtitle( "Osteogenesis") 
	) %>%
		ggsave(
			file = "violin_osteogenesis_for_Chris.pdf",
			device = "pdf",
			useDingbats=FALSE,
			units = c("mm"),
			width = 183 ,
			height = 50 + ( 4 * my_tbl %>% distinct(symbol_ct) %>% count() %>% pull(n) %>% max )
		)
}

#############################################################
# Survival analyses #########################################

my_seed = 1234567

# Gather TCGA information
mycgds = "http://www.cbioportal.org/" %>% cgdsr::CGDS() 


G_list <- 
	biomaRt::getBM(
		filters= "ensembl_gene_id", 
		attributes= c("ensembl_gene_id","hgnc_symbol"),
		values=genes,
		mart= biomaRt::useDataset(
			"hsapiens_gene_ensembl", 
			biomaRt::useMart("ensembl")
		)
	) %>%
	as_tibble()
#merge(df,G_list,by.x="gene",by.y="ensembl_peptide_id")


eee = read_csv("Prostate_Adenocarcinoma_Primary_Tumor.csv") %>%
	separate(sample, c("t1", "t2", "t3"), sep = "-") %>%
	unite(sample, c("t1", "t2", "t3"), sep = "-")  %>%
	
	# Add gene symbol
	left_join(
		(.) %>% 
			distinct(ens_iso) %>%
			separate(ens_iso, c("ens", "iso"), sep="\\.", remove = F) %>%
			mutate(
				symbol = 
					AnnotationDbi:::mapIds(
						org.Hs.eg.db::org.Hs.eg.db,
						keys=ens,
						column="SYMBOL",
						keytype="ENSEMBL",
						multiVals="first"
					)
			) %>%
			distinct(symbol, ens_iso)
	) %>%
	drop_na() %>%
	
	# Average duplicated genes
	left_join( 
		(.) %>% 
			distinct(symbol) %>% 
			mutate(
				part = 1:n() %>% 
					divide_by(length((.))) %>%
					multiply_by(44) %>% 
					ceiling 
			) 
	) %>%
	multidplyr::partition(part) %>%
	do({
		`%>%` = magrittr::`%>%`
		library(tidyverse)
		library(magrittr)
		
		(.) %>%
			group_by(symbol, sample) %>%
			summarise(
				`read count` =  
					`read count` %>% 
					`+` (1) %>% 
					log %>% 
					mean %>% 
					exp %>% 
					`+` (-1) %>% 
					as.integer 
			) %>%
			ungroup() 
		
	}) %>%
	collect() %>%
	ungroup() %>%
	select(-part) %>%
	
	# Normalise
	norm_RNAseq( 
		sample_column = "sample", 
		gene_column = "symbol", 
		value_column = "read count"
	) %>%
	mutate(`read count` = `read count normalised` %>% round) %>%
	select(symbol, sample, `read count`) %>%
	
	# Attach clinical annotation
	inner_join( get_TCGA_prostate_clinical_annotaton() )  %>%
	
	# Add Marker labels
	left_join(
		
		foreach(
			ct_list = 
				list( 
					E = list(tabi_res_E.annot, tabi_gen_E),
					F = list(tabi_res_F.annot, tabi_gen_F),
					T = list(tabi_res_T.annot, tabi_gen_T),
					M = list(tabi_res_M.annot, tabi_gen_M)
				),
			.combine = bind_rows,
			ct = c("E", "F", "T", "M")
		) %dopar% {
			ct_list[[1]] %>%
				filter_too_sparse(ct_list[[2]]) %>%
				filter_out_CI(ct_list[[2]]) %>%
				mutate(is_interfaced = grepl("secreted|membrane", `Protein class`)) %>%
				mutate(`Cell type` = ct) 
		} %>%
			
			# Add more steimates
			mutate(`Safe estimate` = ifelse(`b1.2.5%` > 0, `b1.2.5%`, `b1.97.5%`)) %>%
			rowwise() %>%
			mutate(`Effect size` = effect_size(b0, b1, log_y_cross, T)) %>%
			ungroup() %>%
			rowwise() %>%
			mutate(inflection = get_log_inflection(log_y_cross, inflection, b1)) %>%	 
			mutate(inflection = rescale_inflection(inflection)) %>%	 
			ungroup() %>%
			filter(inflection > 2) %>%
			
			# filter 'benign' genes that are not so expressed in epithelial
			left_join(
				tabi_res_E %>% 
					ungroup() %>% 
					distinct(gene, log_y_cross) %>% 
					rename(log_y_cross_E = log_y_cross, symbol = gene)
			) %>%
			mutate(`Good log ratio to epithelial` = ifelse(`Cell type` == "E", T, log_y_cross > log_y_cross_E)) %>%
			filter(`Good log ratio to epithelial`) %>%
			
			# Calculate the unique score
			mutate(
				`Safe estimate scaled` = abs(`Safe estimate`/max(`Safe estimate`)),
				`Effect size scaled` = sqrt( abs( `Effect size`/max(`Effect size`) ) ) 
			) %>%
			mutate(`Unique score` = sqrt( `Safe estimate scaled`^2 + `Effect size scaled`^2 )) %>%
			
			distinct(`Cell type`, symbol, `Unique score`, `Safe estimate`, `Effect size`, is_interfaced)
	) %>%
	
	# Run deconvolution estimation and merge
	left_join({
		
		my_df = (.)
		
		detach("package:signatureExtractionAlgorithm", unload=TRUE, force=T); library(signatureExtractionAlgorithm)

		deconvoluted_signatures = signatureExtractionAlgorithm_tc(
			mix =
				my_df %>%
				mutate(symbol = as.factor(symbol)) %>%
				rename(gene = symbol) %>%
				distinct(gene , sample,  `read count`) %>%
				spread(sample, `read count`),# %>% select(1:10),
			do_debug = T, save_fit = T,
			TABI_signatures =
				my_df %>%
				filter( `Cell type` %>% is.na %>% `!` ) %>%
				distinct(symbol, `Cell type`)# %>% arrange(symbol) %>% head(n=50)
		)
		
		# save(deconvoluted_signatures, file="fit.obj_ARMET_TABI_signature_just_critical_mean_14_11_2018.RData")
		# load("fit.obj_ARMET_TABI_signature_just_critical_mean_14_11_2018.RData")
		
		deconvoluted_signatures %$% 
			tree %$% 
			signatures_TABI %>%
			bind_rows(
				deconvoluted_signatures %$% 
					tree %$% 
					children %>%
					`[[`(4) %$%
					signatures_TABI
			) %>% 
			#spread( `.variable`, mean) %>%
			select(`Cell type`, symbol, sample, mean) %>%
			right_join(my_df) %>%
			rename(`Deconvoluted` = mean) %>%
			mutate(`Deconvoluted` = exp(`Deconvoluted`) ) 
		
	}) %>%
	
	# Add markers from naive DE analysis
	bind_rows(
		
		(.) %>% 
			inner_join(
				(.) %>%
					mutate(`Is CAPRA-S high` = `CAPRA-S` > 2) %>%
					as_edgeR(
						sample_column = "sample",
						gene_column = "symbol",
						value_column = "read count",
						design_column = "Is CAPRA-S high"
					) %>% 
					filter(FDR < 0.05) %>%
					arrange(FDR) %>%
					head(n=1000) %>%
					select(symbol, FDR)
			) %>%
			mutate(`Unique score` = FDR %>% log %>% multiply_by(-1)) %>%
			select(-FDR) %>%
			mutate(`Cell type` = "ZZ-edgeR")
		
	)  %>%
	
	# Annotate samples with tumor purity
	left_join({
		#load("fit.obj_ARMET_TABI_signature_just_critical_mean_14_11_2018.RData")
		
		deconvoluted_signatures %$%
			proportions %>% 
			gather(`Cell type`, proportion, -c(1:10)) %>% 
			drop_na %>%
			mutate(`Cell type` = ifelse(`Cell type` == "epithelial", "E", `Cell type`)) %>%
			mutate(`Cell type` = ifelse(`Cell type` == "fibroblast", "F", `Cell type`)) %>%
			mutate(`Cell type` = ifelse(`Cell type` == "t_cell", "T", `Cell type`)) %>%
			mutate(`Cell type` = ifelse(`Cell type` %in% c("mono_derived", "granulocyte"), "M", `Cell type`)) %>%
			filter(`Cell type` %in% c("E", "F", "M", "T")) %>%
			group_by(`Cell type`, sample) %>%
			summarise(proportion = proportion %>% sum) %>%
			distinct(sample, `Cell type`, proportion)
	}) %>%
	group_by(`Cell type`) %>%
	mutate(`median proportion` = median(proportion, na.rm=T)) %>%
	ungroup %>%
	mutate(proportion = ifelse(`Cell type` == "ZZ-edgeR", 1, proportion)) %>%
	mutate(`median proportion`  = ifelse(`Cell type` == "ZZ-edgeR", 0, `median proportion` )) %>%
	
	# Annotate with sign of estimate for CAPRA for TCGA
	left_join(
		(.) %>% 
			as_edgeR(
				sample_column = "sample",
				gene_column = "symbol",
				value_column = "read count",
				design_column = "CAPRA-S"
			) %>% 
			select(symbol, logFC, FDR) %>%
			setNames(c("symbol", "edgeR TCGA logFC", "edgeR TCGA FDR"))
	) %>%
	
	# filter and gather
	filter(`Cell type` %>% is.na %>% `!`) %>% 
	
	# Filter interfaced
	#filter(is_interfaced) %>%
	select(-is_interfaced) %>%
	distinct %>%
	
	# Filter cell too rare
	filter(proportion > `median proportion`) %>%
	
	# Filter genes that contradict TCGA trends
	filter(`Safe estimate` * `edgeR TCGA logFC` > 0) %>%
	
	# Correlation
	left_join({
		
		(.) %>%
			filter(`Cell type` %>% is.na %>% `!`) %>%
			group_by(`Cell type`) %>%
			#multidplyr::partition(`Cell type`) %>%
			do({
				
				source("N52_TABI.functions.R")
				source("https://gist.githubusercontent.com/stemangiola/dd3573be22492fc03856cd2c53a755a9/raw/bbcdbfd96bc3a33ea7e35d7e9a6b9ae2e0b6fdce/tidy_extensions.R")
				`%>%` = magrittr::`%>%`
				library(tidyverse)
				library(magrittr)
				library(foreach)
				
				my_seed = 1234567
				set.seed(my_seed)
				tcga_tbl = (.)
				
				tcga_tbl %>%
					
				#######################################################################
				###### CORRELATION ####################################################
				
				left_join(
					(.) %>%
						get_redundant_genes_cor(
							replicates_column = "sample",
							cluster_by_column = "symbol",
							value_column = "read count",
							log_transform = T,
							cor_threshold = 0.6
						) %>%
						rename(is_redundant_TCGA_cor = is_redundant) %>%
						mutate(`Cell type` = tcga_tbl %>% head(n=1) %>% pull(`Cell type`)) 
				) %>%
					
					#######################################################################
				#######################################################################
				
				select(symbol ,  `Cell type`, is_redundant_TCGA_cor) %>%
					distinct()
				
			}) %>%
			collect() %>%
			ungroup()
		
	}) %>%
	
	ungroup() %>%
	
	
	# Label genes based on rank
	left_join(
		(.) %>%
			distinct(
				`Cell type`, symbol, 
				is_redundant_TCGA_cor,
				`Unique score`
			) %>%
			group_by(`Cell type`) %>%
			
			# Correlation rank
			arrange( is_redundant_TCGA_cor,  `Unique score` %>% desc) %>%
			mutate(`Is marker TCGA cor` = row_number() < 43) %>%

			ungroup() %>%
			select(
				`Cell type`, 
				symbol, 
				`Is marker TCGA cor`
			)
	)  %>%
	
	rename(Observed = `read count`) %>%
	gather(Source, `Read count`, c("Observed")) %>%
	gather(
		`Marker selection method`,
		`Is marker`, 
		c( "Is marker TCGA cor"	)
	) %>%
	
	# Filter non sensical combinations 
	filter(`Is marker`) %>%
	
	# Perform tests and return probabilities
	group_by(`Cell type`, Source, `Marker selection method`) %>%
	do({
		
		# Print group
		print((.) %>% head(n=1) %>% select(`Cell type`, Source, `Marker selection method`) %>% as.character %>% paste(collapse=" "))
		
		# Is needed in nested blocks
		input = (.) %>% 
			mutate(`Read count` = ifelse(`Read count`==0, 1, `Read count`)) %>%
			mutate(`Read count` = sign(`Read count`)*log(abs(`Read count`))) %>%
			mutate(symbol       = gsub("-", "_", symbol)) %>%
			
			# Eliminate when too many 0s because script fails
			{
				if((.) %>% filter(`Read count` == 0) %>% nrow() > 0)
					(.) %>%
					inner_join(
						(.) %>% 
							mutate(is_zero = `Read count` == 0) %>% 
							group_by(symbol) %>% 
							count(is_zero) %>% 
							spread(is_zero, n) %>%
							replace_na(list(`FALSE` = 0, `TRUE` = 0)) %>%
							filter(`FALSE` / `TRUE` > 0.3) %>%
							select(symbol)
					)
				else (.)
			} %>%
			
			# NEW ERROR THAT SHOULD NOT BE HERE
			filter(`Read count` %>% is.na %>% `!`) 
		###################################
		

		proportion_covariate = 
			ifelse(
				(.)  %>% head(n=1) %>% pull(`Cell type`) %in% c("E", "F", "M", "T"), 
				"+ proportion", "")
		
		survival::coxph(
			formula = formula(
				sprintf(
					"survival::Surv(DFS_MONTHS,is_recurred)~ CAPRA_S %s + %s", 
					proportion_covariate,
					input %>%
						distinct(symbol) %>%
						pull(symbol) %>%
						paste(collapse=" + ")
					
				)
			),
			data = input %>%
				rename(CAPRA_S = `CAPRA-S`) %>%
				distinct(sample, symbol, `Read count`, DFS_MONTHS,is_recurred, CAPRA_S, proportion) %>%
				group_by(symbol) %>%
				mutate(`Read count` = `Read count` %>% scale) %>%
				ungroup() %>%
				mutate(proportion = proportion %>% scale, CAPRA_S = CAPRA_S %>% scale)  %>%
				spread(symbol, `Read count`)
		) %>%
			
			# Parse results
			print_summary_return() %>%
			
			# Stratify
			summary() %$% 
			coefficients %>%
			as_tibble(rownames = "symbol") %>% 
			
			# Filter genes non consistent between TCGA and N52 but include CAPRA
			{
				if(input %>% head(n=1) %>% pull(`Cell type`) == "ZZ-edgeR") (.)
				else
					(.) %>% 
					left_join(input %>% distinct(symbol, `Safe estimate`) %>% rename(`TABI estimate` = `Safe estimate`)) %>%
					# Filter not for ZZ-edgeR
					filter(coef * `TABI estimate` > 0 | symbol == "CAPRA_S")
				
			}  %>%
			
			mutate(`Is bad` = ifelse(coef>0, T, F)) %>%
			select(symbol, `Is bad`, `Pr(>|z|)`, coef, `se(coef)`) %>%
			rename(`Estimate cox` = coef) %>%
			left_join(input) %>%
			
			# Add annotation for CAPRA-s
			mutate(`Cell type` = ifelse(
				symbol == "CAPRA_S",
				input %>% head(n=1) %>% pull(`Cell type`), 
				`Cell type`
			)
			) %>%
			mutate(`Marker selection method` = ifelse(
				symbol == "CAPRA_S",
				input %>% head(n=1) %>% pull(`Marker selection method`), 
				`Marker selection method`
			)
			) %>%
			mutate(Source = ifelse(
				symbol == "CAPRA_S",
				input %>% head(n=1) %>% pull(Source), 
				Source
			)
			) %>%
			filter(symbol != "proportion")
		
	}) %>%
	
	# Filter significant genes
	filter(`Pr(>|z|)` < 0.1 | grepl("CAPRA", symbol)) %>%
	filter(`Cell type` %>% is.na %>% `!`) %>%
	filter(symbol != "CAPRA_S") %>%
	
	do({	
			
			# Print the three genes
			print_return %>%
			
			# Select thresholds
			{
				#browser()
				set.seed(my_seed)
				
				input_3 = (.) %>%
					
					# Eliminate when too many 0s because script fails
				{
					if((.) %>% filter(`Read count` == 0) %>% nrow() > 0)
						(.) %>%
						inner_join(
							(.) %>% 
								mutate(is_zero = `Read count` == 0) %>% 
								group_by(symbol) %>% 
								count(is_zero) %>% 
								spread(is_zero, n) %>%
								replace_na(list(`FALSE` = 0, `TRUE` = 0)) %>%
								filter(`FALSE` / `TRUE` > 0.3) %>%
								select(symbol)
						)
					else (.)
				}
				
				genes = input_3 %>% distinct(symbol) %>% pull(symbol)
				
				
				#samples_training = (.) %>% distinct(sample) %>% sample_frac(0.5) %>% pull(sample)
				
				# Ignore emptydata sets
				if((.) %>% nrow == 0) (.)
				else 
					
					input_3 %>% 
					#filter(sample %in% samples_training %>% `!`) %>%
					inner_join( (.) %>% distinct(sample) %>% sample_frac(0.5)) %>%
					# Add category for genes calculating cut points
					left_join(
						input_3 %>%
							#filter(sample %in% samples_training) %>%
							inner_join( (.) %>% distinct(sample) %>% sample_frac(0.5)) %>%
							distinct(symbol, `Read count`, sample, DFS_MONTHS, is_recurred) %>%
							spread(symbol, `Read count`) %>% 
							
							# Calculate cutpoints
							survminer::surv_cutpoint(
								time = "DFS_MONTHS",
								event = "is_recurred",
								variables = genes
							) %$% 
							cutpoint %>% 
							as_tibble(rownames="symbol") 
					) %>%
					
					# Label
					mutate(category = ifelse(`Read count` > cutpoint, "high", "low")) %>%
					
					# classify whether is bad and count
					mutate(`Bad event` = ifelse( 
						(`Is bad` & category == "high") | 
							(!`Is bad` & category == "low"), 
						1, 0
					)) %>%
					
					# Add category for CAPRA
					mutate(`CAPRA-S category` = 
								 	ifelse(
								 		`CAPRA-S` >
								 			input_3 %>% 
								 			inner_join( (.) %>% distinct(sample) %>% sample_frac(0.5)) %>%
								 			rename(CAPRA_S = `CAPRA-S`) %>%
								 			# Calculate cutpoints
								 			survminer::surv_cutpoint(
								 				time = "DFS_MONTHS",
								 				event = "is_recurred",
								 				variables = "CAPRA_S"
								 			) %$% 
								 			cutpoint %>% 
								 			as_tibble %>%
								 			pull(cutpoint),
								 		"high",
								 		"low"
								 	)
					) %>%
					
					# classify whether is bad and count
					mutate(`CAPRA-S Bad event` = ifelse(`CAPRA-S category` == "high" , 1, 0))
			}
		
	}) %>%
	
	# Add plot titles
	do(
		(.) %>%
			# Add title for later
			mutate(
				`Genes bad/good` = 
					(.) %>%
					ungroup() %>%
					distinct(symbol, `Is bad`) %>% 
					mutate(is = ifelse(`Is bad`, "bad", "good")) %>% 
					select(-`Is bad`) %>% 
					unite(name) %>% 
					pull(name) %>% 
					paste(collapse=" ")
			) 
	) %>%
	mutate(`\n` = "\n") %>%
	mutate(`\n 2` = "\n") %>%
	unite(title, c("Cell type", "Source", "\n", "Marker selection method", "\n 2", "Genes bad/good"), sep=" - ", remove = F) %>%
	ungroup()  %>%
	
	#Trend plots
  do_and_return(function(){

		(.) %>%
			ungroup() %>%
			filter(Source == "Observed" & `Marker selection method` == "Is marker TCGA cor") %>%
			distinct(`Cell type`, symbol, `Estimate cox` ) %>%
			left_join(
				tabi_res_E.annot %>% mutate(`Cell type` = "E") %>%
					bind_rows(tabi_res_F.annot %>% mutate(`Cell type` = "F")) %>%
					bind_rows(tabi_res_M.annot %>% mutate(`Cell type` = "M")) %>%
					bind_rows(tabi_res_T.annot %>% mutate(`Cell type` = "T"))
			) %>%
			unite(symbol, c("Cell type", "symbol"), remove = F) %>%

			print_gene_regression(

				(.) %>%
					distinct(`Cell type`, symbol) %>%
					left_join(
						tabi_gen_E %>% mutate(`Cell type` = "E") %>%
							bind_rows(tabi_gen_F %>% mutate(`Cell type` = "F")) %>%
							bind_rows(tabi_gen_M %>% mutate(`Cell type` = "M")) %>%
							bind_rows(tabi_gen_T %>% mutate(`Cell type` = "T")) %>%
							ungroup()
					) %>%
					unite(symbol, c("Cell type", "symbol"), remove = F),

				tabi_inflection_E %>% mutate(gene = sprintf("E_%s", gene)) %>%
					bind_rows(tabi_inflection_F %>% mutate(gene = sprintf("F_%s", gene)) ) %>%
					bind_rows(tabi_inflection_M %>% mutate(gene = sprintf("M_%s", gene)) ) %>%
					bind_rows(tabi_inflection_T %>% mutate(gene = sprintf("T_%s", gene)) ) %>%
					mutate(symbol = gene),

				multiple_ct = T,
				num_cols = 6,
				do_filter = F,
				gene_order = 	(.) %>%
					distinct(symbol, `Cell type`, `Estimate cox`) %>%
					arrange(`Cell type`, `Estimate cox` %>% abs %>% desc) %>%
					pull(symbol)
			) %>%
			ggsave(
				"trend_plots_for_classifier_survival.pdf",
				plot = .,
				device = "pdf",
				useDingbats=FALSE,
				units = c("mm"),
				width = 223 ,
				height = 223
			)

		(.)


	}) %>%
	
	# Count how many bad
	group_by(`Cell type`, Source,  `Marker selection method`, sample, `CAPRA-S`, `CAPRA-S category`, DFS_MONTHS, is_recurred, title) %>%
	summarise( `How_many_bad` = `Bad event` %>% sum	) %>%
	ungroup() %>%
	
	# Plot survival just with CAPRA-S
	do_and_return(function(){
		
		input_4 = (.) %>% 
			rename(CAPRA_S_cat = `CAPRA-S category`) %>%
			distinct(sample, CAPRA_S_cat, DFS_MONTHS,is_recurred)
		
		survival::survfit(
			survival::Surv(DFS_MONTHS,is_recurred) ~ CAPRA_S_cat, 
			data=input_4
		) %>%
			survminer::ggsurvplot(
				fit=., 
				input_4, 
				risk.table = FALSE, conf.int = T,
				palette = c("#ed6f68",  "#5366A0" ),
				legend = "none",
				pval = T,
				#title = "Only CAPRA-S", 
				ggtheme = my_theme + theme(text = element_text(size=8))
			) %>% 
			list(.) %>%
			survminer::arrange_ggsurvplots(print = FALSE, ncol = 1, nrow=1) %>%
			ggsave(
				"survival_plot_CAPRA_S_publication.pdf", 
				plot = .,
				device = "pdf",
				useDingbats=FALSE,
				units = c("mm"),
				width = 183 ,
				height = 183
			)
		
	}) %>%
	
	# Plot survival gradual
	do_and_return(function(){
		
		(.) %>%
			group_by(`Cell type`, Source, `Marker selection method`) %>% 
			do({
				
				input_2 = (.) %>%
					
					# Eliminate groups with n == 1 because causes error 
					group_by(How_many_bad) %>%
					mutate(n = n()) %>%
					ungroup() %>%
					filter(n > 1)
				
				tibble(
					title = input_2 %>% head(n=1) %>% pull(title),
					plots = list(
						
						survival::survfit(
							survival::Surv(DFS_MONTHS,is_recurred) ~ How_many_bad, 
							data=input_2
						) %>%
							survminer::ggsurvplot(
								fit=., 
								input_2, 
								risk.table = FALSE, conf.int = T,# pval = T,
								title = input_2 %>% head(n=1) %>% pull(title), 
								#	pval = T,
								test.for.trend = T,
								ggtheme = my_theme + theme(text = element_text(size=8))
								
							) 
					)
				)
				
			}) %>%
			ungroup() %>%
			pull(plots) %>%
			survminer::arrange_ggsurvplots(ncol = 5, nrow = 9) %>%
			ggsave(
				filename = sprintf("survival_plot_gradual.pdf"),
				device = "pdf",
				useDingbats=FALSE,
				units = c("mm"),
				width = 483 ,
				height = 583
			)
		
	}) %>%
	
	# Plot survival gradual publication
	do_and_return(function(){
		
		(.) %>%
			filter(Source == "Observed" & `Marker selection method` == "Is marker TCGA cor") %>%
			group_by(`Cell type`, Source, `Marker selection method`) %>% 
			do({
				
				input_2 = (.) %>%
					
					# Eliminate groups with n == 1 because causes error 
					group_by(How_many_bad) %>%
					mutate(n = n()) %>%
					ungroup() %>%
					filter(n > 1)
				
				palette = colorRampPalette(
					c("white",
						switch(
							input_2 %>% head(n=1) %>% pull(`Cell type`),
							"E" = "#199e77",
							"F" = "#d86012",
							"T" = "#e42e89",
							"M" = "#7570b2"
						), 
						"black")
				)(
					(.) %>% distinct(How_many_bad) %>% nrow %>% sum(4)
				) %>% 
					`[` (-c(1:2)) %>% rev %>% `[` (-c(1:2)) %>% rev
				
				tibble(
					title = input_2 %>% head(n=1) %>% pull(title),
					plots = list(
						survival::survfit(
							survival::Surv(DFS_MONTHS,is_recurred) ~ How_many_bad, 
							data=input_2
						) %>%
							survminer::ggsurvplot(
								fit=., 
								input_2, 
								risk.table = FALSE, 
								conf.int = F,
								legend = "none",
								#title = input_2 %>% head(n=1) %>% pull(title), 
								#pval = T,
								palette = palette,
								test.for.trend = T,
								ggtheme = my_theme + theme(text = element_text(size=8), legend.position='none')
								
							) 
					)
				)
				
			}) %>%
			ungroup() %>%
			pull(plots) %>%
			survminer::arrange_ggsurvplots(ncol = 1, nrow = 4) %>%
			ggsave(
				filename = sprintf("survival_plot_gradual_publication.pdf"),
				device = "pdf",
				useDingbats=FALSE,
				units = c("mm"),
				width = 400/4 ,
				height = 183
			)
		
	} ) %>%
	
	# Further stratify
	group_by(`Cell type`, Source,  `Marker selection method`) %>%
	do({
		
		(.) %>%
			
			# Add thresholds
			mutate(
				min = (.) %>% count(How_many_bad) %>% filter(n > 5) %>% arrange(How_many_bad) %>% head(n=(.)%>%nrow %>% divide_by(3) %>%ceiling) %>% tail(n=1) %>% pull(How_many_bad) ,
				max = (.) %>% count(How_many_bad) %>% filter(n > 5) %>% arrange(How_many_bad) %>%  tail(n=(.)%>%nrow %>% divide_by(3) %>% ceiling) %>% head(n=1) %>% pull(How_many_bad)
			) %>%
			
			mutate(
				`How_many_bad category` = 
					ifelse(
						How_many_bad <= min,
						"low",
						ifelse(
							How_many_bad >= max,
							"high",
							NA
						)
					)
			) %>%
			
			# If a category includes 1 sample drop the data frame
			{
				if( ((.) %>% count(`How_many_bad category`) %>% drop_na %>% pull(n) %>% min) < 2) (.) %>% head(n=0)
				else (.)
			}
		
		# mutate(
		# 	`How_many_bad category` =
		# 		# This is for avoiding strange error of the function
		# 		tryCatch({
		# 			(.) %>%
		# 				survminer::surv_cutpoint(
		# 					time = "DFS_MONTHS",
		# 					event = "is_recurred",
		# 					variables = "How_many_bad"
		# 				) %>%
		# 				survminer::surv_categorize() %>%
		# 				as_tibble() %>%
		# 				pull(How_many_bad)
		# 		},
		# 		error = function(e) {
		# 			sapply(
		# 				( (.) %>% pull(How_many_bad) > ( (.) %>% pull(How_many_bad) %>% median() %>% round) ) %>% as.numeric %>% `+` (1),
		# 				switch,
		# 				"low",
		# 				"high"
		# 			)
		# 		}
		# 		)
		# )
	}) %>%
	
	# Plot just genes categorised
	do_and_return(function(){
		
		(.) %>% do({
			
			input_3 = (.) %>% 
				rename(CAPRA_S_cat = `CAPRA-S category`) %>%
				rename(how_many_bad_genes = `How_many_bad category`) %>%
				mutate(
					title = 
						(.) %>% 
						head(n=1) %>% 
						select(`Cell type`, Source,  `Marker selection method`) %>% 
						as.character %>% 
						paste(collapse=" ")
				) %>% mutate(how_many_bad_genes = ifelse(how_many_bad_genes %>% is.na, "middle", how_many_bad_genes))
			
			tibble(
				title = input_3 %>% head(n=1) %>% pull(title),
				plots = list(
					survival::survfit(
						survival::Surv(DFS_MONTHS,is_recurred) ~ how_many_bad_genes, 
						data=input_3 
					) %>%
						survminer::ggsurvplot(
							fit=., 
							input_3, 
							risk.table = FALSE, conf.int = T, pval = T,
							title = input_3 %>% head(n=1) %>% pull(title), 
							ggtheme = my_theme + theme(text = element_text(size=8))
						) 
				)
			)
			
		}) %>%
			ungroup() %>%
			pull(plots) %>%
			survminer::arrange_ggsurvplots(ncol = 5, nrow = 9) %>%
			ggsave(
				filename = sprintf("survival_plot_stratified.pdf"),
				device = "pdf",
				useDingbats=FALSE,
				units = c("mm"),
				width = 483 ,
				height = 583
			)
		
	}) %>%
	
	ungroup() %>%
	do_and_return(function(){
		(.) %>% 
			filter(Source == "Observed" & `Marker selection method` == "Is marker TCGA cor") %>%
			#filter(`Cell type` != "ZZ-edgeR") %>%
			group_by(`Cell type`) %>%
			do({
				
				input_3 = (.) %>% 
					rename(CAPRA_S_cat = `CAPRA-S category`) %>%
					rename(how_many_bad_genes = `How_many_bad category`) %>%
					mutate(
						title = 
							(.) %>% 
							head(n=1) %>% 
							select(`Cell type`, Source,  `Marker selection method`) %>% 
							as.character %>% 
							paste(collapse=" ")
					) 
				#mutate(how_many_bad_genes = ifelse(how_many_bad_genes %>% is.na, "middle", how_many_bad_genes)) %>%
				
				palette = colorRampPalette(
					c("white",
						switch(
							input_3 %>% head(n=1) %>% pull(`Cell type`),
							"E" = "#199e77",
							"F" = "#d86012",
							"T" = "#e42e89",
							"M" = "#7570b2"
						), 
						"black")
				)(
					(.) %>% distinct(How_many_bad) %>% nrow %>% sum(4)
				) %>% 
					`[` (-1) %>% rev %>% `[` (-1) %>% rev %>%
					`[` (c(4,2))
				
				tibble(
					title = input_3 %>% head(n=1) %>% pull(title),
					plots = list(
						survival::survfit(
							survival::Surv(DFS_MONTHS,is_recurred) ~ how_many_bad_genes, 
							data=input_3 
						) %>%
							survminer::ggsurvplot(
								fit=., 
								input_3, 
								risk.table = FALSE, conf.int = T,
								palette = palette,
								legend = "none",
								pval = T,
								#title = input_3 %>% head(n=1) %>% pull(title), 
								ggtheme = my_theme + theme(text = element_text(size=8))
							) 
					)
				)
			}) %>%
			ungroup() %>%
			pull(plots) %>%
			survminer::arrange_ggsurvplots(ncol = 2, nrow = 4) %>%
			ggsave(
				filename = sprintf("survival_plot_stratified_publication.pdf"),
				device = "pdf",
				useDingbats=FALSE,
				units = c("mm"),
				width = 400/4 ,
				height = 183
			)
		
	}) %>%
	
	# Plot survival combination
	#do_and_return(function() 
	do({
		
		plot_list = list()
		
		input_5 = 
			
			(.) %>%
			filter(Source == "Observed" & `Marker selection method` == "Is marker TCGA cor") %>%
			ungroup() %>% 
			filter(`Cell type` != "ZZ-edgeR") %>%
			
			# Eliminate non significant Cell types
			######################################
		
		distinct(sample, `Cell type`, `How_many_bad category`, DFS_MONTHS, is_recurred, `CAPRA-S category`)  %>%
			mutate(`How_many_bad category` = ifelse(`How_many_bad category`=="low", -1, 1)) %>%
			mutate(`How_many_bad category` = as.numeric(`How_many_bad category`)) %>%
			spread(`Cell type`, `How_many_bad category`) %>%
			replace(is.na(.), 0) 
		
		
		# For epithelial
		input_7 = input_5 %>%
			mutate(`Category tot` = E) %>%
			mutate(`Category tot` = ifelse(`Category tot`>0, 1, `Category tot`)) %>%
			mutate(`Category tot` = ifelse(`Category tot`<0, -1, `Category tot`)) %>%
			unite(cat, c("CAPRA-S category", "Category tot")) %>% 
			mutate(cat = ifelse(grepl("0", cat) | cat =="low_1" | cat=="high_-1", "middle", cat)) 
		
		plot_list[[1]] =
			input_7 %>%
			survival::survfit(
				survival::Surv(DFS_MONTHS,is_recurred) ~ cat, 
				data=. 
			) %>%
			survminer::ggsurvplot(
				fit=., 
				input_7, 
				risk.table = FALSE, conf.int = T,
				palette = c("#ed6f68",  "#5366A0", "#fdcb68"  ),
				legend = "none",
				pval = T,
				#title = "CAPRA + Genes", 
				ggtheme = my_theme + theme(text = element_text(size=8), legend.position='none')
			) 
		
		# For epithelial
		input_8 = input_5 %>%
			mutate(`Category tot` = M + F + T) %>%
			mutate(`Category tot` = ifelse(`Category tot`>0, 1, `Category tot`)) %>%
			mutate(`Category tot` = ifelse(`Category tot`<0, -1, `Category tot`)) %>%
			unite(cat, c("CAPRA-S category", "Category tot")) %>% 
			mutate(cat = ifelse(grepl("0", cat) | cat =="low_1" | cat=="high_-1", "middle", cat)) 
		
		plot_list[[2]] =
			input_8 %>%
			survival::survfit(
				survival::Surv(DFS_MONTHS,is_recurred) ~ cat, 
				data=. 
			) %>%
			survminer::ggsurvplot(
				fit=., 
				input_8, 
				risk.table = FALSE, conf.int = T,
				palette = c("#ed6f68",  "#5366A0", "#fdcb68"  ),
				legend = "none",
				pval = T,
				#title = "CAPRA + Genes", 
				ggtheme = my_theme + theme(text = element_text(size=8), legend.position='none')
			) 

		# For edgeR
		input_9 = 
			(.) %>%
			filter(Source == "Observed" & `Marker selection method` == "Is marker TCGA cor") %>%
			ungroup() %>% 
			filter(`Cell type` == "ZZ-edgeR") %>%

			distinct(sample, `Cell type`, `How_many_bad category`, DFS_MONTHS, is_recurred, `CAPRA-S category`)  %>%
			mutate(`How_many_bad category` = ifelse(`How_many_bad category`=="low", -1, 1)) %>%
			mutate(`How_many_bad category` = as.numeric(`How_many_bad category`)) %>%
			spread(`Cell type`, `How_many_bad category`) %>%
			replace(is.na(.), 0) %>%
			mutate(`Category tot` = `ZZ-edgeR`) %>%
			mutate(`Category tot` = ifelse(`Category tot`>0, 1, `Category tot`)) %>%
			mutate(`Category tot` = ifelse(`Category tot`<0, -1, `Category tot`)) %>%
			unite(cat, c("CAPRA-S category", "Category tot")) %>% 
			mutate(cat = ifelse(grepl("0", cat) | cat =="low_1" | cat=="high_-1", "middle", cat)) 
		
		plot_list[[3]] =
			input_9 %>%
			survival::survfit(
				survival::Surv(DFS_MONTHS,is_recurred) ~ cat, 
				data=. 
			) %>%
			survminer::ggsurvplot(
				fit=., 
				input_9, 
				risk.table = FALSE, conf.int = T,
				palette = c("#ed6f68",  "#5366A0", "#fdcb68"  ),
				legend = "none",
				pval = T,
				#title = "CAPRA + Genes", 
				ggtheme = my_theme + theme(text = element_text(size=8), legend.position='none')
			) 
		
		plot_list %>%
			survminer::arrange_ggsurvplots(print = FALSE, ncol = 3, nrow=1) %>%
			ggsave(
				"survival_plot_publication.pdf", 
				plot = .,
				device = "pdf",
				useDingbats=FALSE,
				units = c("mm"),
				width = 183 ,
				height = 183
			)
	})

#############################################################
# ARMET in TCGA #############################################

mycgds = "http://www.cbioportal.org/" %>% cgdsr::CGDS() 

#library('biomaRt')
mart <- 
	genes <- c("ENSG00000000457", "ENSG00000000003")

G_list <- 
	biomaRt::getBM(
		filters= "ensembl_gene_id", 
		attributes= c("ensembl_gene_id","hgnc_symbol"),
		values=genes,
		mart= biomaRt::useDataset(
			"hsapiens_gene_ensembl", 
			biomaRt::useMart("ensembl")
		)
	) %>%
	as_tibble()
#merge(df,G_list,by.x="gene",by.y="ensembl_peptide_id")


TCGA_tbl = read_csv("Prostate_Adenocarcinoma_Primary_Tumor.csv") %>%
	separate(sample, c("t1", "t2", "t3"), sep = "-") %>%
	unite(sample, c("t1", "t2", "t3"), sep = "-")  %>%
	
	# Add gene symbol
	left_join(
		(.) %>% 
			distinct(ens_iso) %>%
			separate(ens_iso, c("ens", "iso"), sep="\\.", remove = F) %>%
			mutate(
				symbol = 
					AnnotationDbi:::mapIds(
						org.Hs.eg.db::org.Hs.eg.db,
						keys=ens,
						column="SYMBOL",
						keytype="ENSEMBL",
						multiVals="first"
					)
			) %>%
			distinct(symbol, ens_iso)
	) %>%
	drop_na() %>%
	
	# Average duplicated genes
	do_parallel_start(40, "symbol") %>%
	do({
		`%>%` = magrittr::`%>%`
		library(tidyverse)
		library(magrittr)
		
		(.) %>%
			group_by(symbol, sample) %>%
			summarise(
				`read count` =  
					`read count` %>% 
					`+` (1) %>% 
					log %>% 
					mean %>% 
					exp %>% 
					`+` (-1) %>% 
					as.integer 
			) %>%
			ungroup() 
	}) %>%
	do_parallel_end() %>%
	
	# Normalise
	norm_RNAseq( 
		sample_column = "sample", 
		gene_column = "symbol", 
		value_column = "read count"
	) %>%
	mutate(`read count` = `read count normalised` %>% round) %>%
	select(symbol, sample, `read count`) %>%
	
	# Attach clinical annotation
	inner_join( get_TCGA_prostate_clinical_annotaton()	) 



library(ARMET)

# Model on CAPRA-S
fit_TCGA_CAPRA = ARMET_tc(
	mix =
		TCGA_tbl %>% 
		distinct(sample, symbol, `read count`) %>%
		spread(sample, `read count`) %>%
		filter(symbol %>% is.na %>% `!`) %>%
		rename(gene = symbol), 
	my_design = 
		TCGA_tbl %>% 
		distinct(sample, `CAPRA-S`) %>% 
		mutate(`(Intercept)` = 1) %>%
		select(1,3,2), 
	cov_to_test = "CAPRA-S",
	save_fit = F
)
save(fit_TCGA_CAPRA, file="ARMET_fit_TCGA_CAPRA.RData")

# Survival model
# fit_TCGA = ARMET_tc(
# 	mix =
# 		TCGA_tbl %>% distinct(sample, symbol, `read count`) %>%
# 		spread(sample, `read count`) %>%
# 		filter(symbol %>% is.na %>% `!`) %>%
# 		rename(gene = symbol),
# 	save_fit = F
# )
# save(fit_TCGA, file="ARMET_fit_TCGA.RData")
load("ARMET_fit_TCGA.RData")

prop = 
	fit_TCGA %$%
	proportions %>%
	filter(ct %in% c("epithelial", "endothelial", "fibroblast")) %>%
	bind_rows(
		fit_TCGA %$%
			proportions %>%
			filter(ct %in% c("macrophage_M0", "macrophage_M1", "macrophage_M2")) %>%
			group_by(sample) %>% 
			summarise(absolute_proportion = absolute_proportion %>% sum) %>%
			mutate(ct = "mono_derived")
	) %>%
	bind_rows(
		fit_TCGA %$%
			proportions %>%
			filter(grepl("t_", ct)) %>%
			group_by(sample) %>% 
			summarise(absolute_proportion = absolute_proportion %>% sum) %>%
			mutate(ct = "t")
	) %>%
	bind_rows(
		fit_TCGA %$%
			proportions %>%
			filter(ct %in% c("eosinophil", "neutrophil")) %>%
			group_by(sample) %>% 
			summarise(absolute_proportion = absolute_proportion %>% sum) %>%
			mutate(ct = "granulocyte")
	) %>%
	bind_rows(
		fit_TCGA %$%
			proportions %>%
			filter(grepl("b_", ct)) %>%
			group_by(sample) %>% 
			summarise(absolute_proportion = absolute_proportion %>% sum) %>%
			mutate(ct = "b")
	) %>%
	bind_rows(
		fit_TCGA %$%
			proportions %>%
			filter(grepl("nk_", ct)) %>%
			group_by(sample) %>% 
			summarise(absolute_proportion = absolute_proportion %>% sum) %>%
			mutate(ct = "nk")
	) 

# Survival analyses
survival::coxph(
	formula = formula(
		sprintf(
			"survival::Surv(DFS_MONTHS,is_recurred)~ %s", 
			prop %>%
				distinct(ct) %>%
				pull(ct) %>%
				paste(collapse=" + ")
			
		)
	),
	data =prop %>% 
		select(sample, ct, absolute_proportion) %>% 
		spread(ct, absolute_proportion) %>% 
		mutate_if(is.numeric, scale) %>%
		left_join(
			TCGA_tbl %>% distinct(`CAPRA-S`, sample,DFS_MONTHS,is_recurred) %>% rename(CAPRA_S = `CAPRA-S`)
		)
) %>%
	
	# Parse results
	#print_summary_return() %>%
	
	# Stratify
	summary() %$% 
	coefficients %>%
	as_tibble(rownames = "symbol") %>% 
	mutate(`Is bad` = ifelse(coef>0, T, F)) %>%
	select(symbol, `Is bad`, `Pr(>|z|)`, coef, `se(coef)`) %>%
	rename(`Estimate cox` = coef)



prop %>% 
	left_join(
		TCGA_tbl %>% distinct(`CAPRA-S`, sample,DFS_MONTHS,is_recurred) %>% rename(CAPRA_S = `CAPRA-S`)
	) %>%
	group_by(ct) %>%
	do({
		set.seed(78537)
		my_prop = (.)
		
		samples_training = my_prop %>% distinct(sample) %>% sample_frac(0.5) %>% pull(sample)
		
		input = my_prop %>%
			filter(sample %in% samples_training %>% `!`) %>%
			#inner_join( (.) %>% distinct(sample) %>% sample_frac(0.5)) %>%
			# Add category for genes calculating cut points
			left_join(
				my_prop %>%
					filter(sample %in% samples_training) %>%
					#inner_join( (.) %>% distinct(sample) %>% sample_frac(0.5)) %>%
					distinct(sample,absolute_proportion, DFS_MONTHS, is_recurred, ct) %>%
					
					# Calculate cutpoints
					survminer::surv_cutpoint(
						time = "DFS_MONTHS",
						event = "is_recurred",
						variables = "absolute_proportion"
					) %$% 
					cutpoint %>% 
					as_tibble(rownames="ct") %>% mutate(ct = my_prop %>% head(n=1) %>% pull(ct))
			) %>%
			
			# Label
			mutate(category = ifelse(absolute_proportion > cutpoint, "high", "low")) 
		
		
		tibble(
			ct = my_prop %>% head(n=1) %>% pull(ct),
			plot = list(
				survival::survfit(
					survival::Surv(DFS_MONTHS,is_recurred) ~ category, 
					data=input 
				) %>%
					survminer::ggsurvplot(
						fit=., 
						input, 
						risk.table = FALSE, conf.int = T,
						palette = c("#CE2D29",  "#5265A0" ),
						legend = "none",
						title = my_prop %>% head(n=1) %>% pull(ct), 
						ggtheme = my_theme + theme(text = element_text(size=8)), 
						pval = T
						#, 
						#ggtheme = theme(axis.title.x=element_blank(), axis.title.y=element_blank())
					)
			)
		)
	}) %>% 
	pull(plot) %>%
	survminer::arrange_ggsurvplots(print = F, ncol = 3, nrow=3) %>%
	ggsave(
		"survival_plot_all_cells_TCGA_publication.pdf", 
		plot = .,
		device = "pdf",
		useDingbats=FALSE,
		units = c("mm"),
		width = 183 ,
		height = 183
	)

#############################################################
# Print examples of trends ##################################

tabi_res_E.annot %>%
	filter_too_sparse(tabi_gen_E) %>%
	filter_out_CI(tabi_gen_E, 13/3) %>%
	mutate(`Safe estimate` = ifelse(`b1.2.5%` > 0, `b1.2.5%`, `b1.97.5%`)) %>%
	mutate(inflection = get_log_inflection(log_y_cross, inflection, b1)) %>%
	filter(`Safe estimate` * estimate > 0) %>%
	arrange(`Safe estimate` %>% abs %>% desc) %>%
	distinct(symbol, estimate, inflection) %>%
	inner_join(summary_df %>% filter(`Cell type` == "E") %>% distinct(symbol)) %>%
	View

tabi_res_E.annot %>%
	
	filter(symbol %in% c("AVPR1A", "LUM", "SPNS2", "SPINK5")) %>%
	
	print_gene_regression(
		tabi_gen_E %>%
			ungroup() %>%
			filter(gene %in% (tabi_res_E.annot %>% pull(symbol) %>% unique()) ) %>% 
			mutate(symbol = gene),
		tabi_inflection_E %>% mutate(symbol = gene), 
		do_filter = F, num_cols = 5
	)


summary_df %>%
	filter(symbol %in% c("AVPR1A", "LUM", "SPNS2", "SPINK5") & `Cell type` == "E") %>%
	
	{
		
		# Arrange symbol_ct
		(.) %>% 
			
			# Plot
			ggplot( aes(x =symbol_ct , y=inflection, fill=`Effect size`)) +
			geom_violin(
				scale = "width", adjust = 1, width = 1
				
			) +
			scale_fill_distiller(
				palette = "Spectral",
				na.value = 'white',
				direction = 1,
				trans = signed_log,
				breaks = signed_log_breaks,
				labels = scales::scientific,
				limits=c(
					-max(abs((.) %>% pull(`Effect size`))),
					max(abs((.) %>% pull(`Effect size`)))
				)
			) +
			geom_point( 
				data = (.) %>% 
					group_by(symbol_ct) %>% 
					summarise(
						i = quantile(inflection, 0.5), 
						`Effect size` = unique(`Effect size`)
					),
				aes(y = i),
				size = 0.5
			) +
			geom_boxplot(
				data = (.) %>% 
					distinct(
						symbol_ct, 
						`Effect size`,
						`Cell type`, 
						inflection, 
						Category,
						positively_associated
					) %>% 
					mutate(lower = -20, upper = -19),
				aes(x=symbol_ct, ymin = lower, lower = lower, middle = lower, upper = upper, ymax = upper, color=`Cell type`),
				stat = "identity"
			) +
			scale_color_brewer(palette="Dark2") +
			coord_flip(
				ylim=c(-25,17),
				xlim=c(1, (.) %>% distinct(symbol_ct, positively_associated) %>% group_by(positively_associated) %>% summarise(n = n()) %>% pull(n) %>% max )
			)+
			facet_wrap(
				~ symbol, 
				scales = "free", 
				#space = "free",
				drop = T, ncol = 4
			) +
			scale_y_continuous(
				trans = magnify_trans(interval_low = 0, interval_high = 7,  reducer = 20),
				#trans = magnify_trans(intercept =  7,  reducer = 20),
				breaks = c(-10, 0, 1:7, 17)
			) +
			theme_bw() +
			theme(
				panel.border = element_blank(), 
				axis.line = element_line(),
				panel.grid.major = element_line(size = 0.2),
				panel.grid.minor = element_line(size = 0.1),
				text = element_text(size=12),
				legend.position="bottom",
				legend.text = element_text(angle = 60, vjust = 0.5),
				axis.text.x = element_text(angle = 60, vjust = 0.5),
				strip.background = element_blank(),
				axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
				axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
				
				
			) + 
			ylab("CAPRA score") + xlab("Genes") +
			guides(
				color = guide_legend(keywidth = 0.5, keyheight = 0.5, direction = "vertical"),
				fill = guide_colorbar(barwidth = 5, barheight = 0.3, title.vjust = 1)
			) 
	} 

#############################################################
# Creation of patient clinical table ########################



#############################################################
# Creation of EGA files for publication #####################

#java -jar ~/third_party_sofware/EGA_submission/EgaCryptor/EgaCryptor/EgaCryptor.jar -file *fastq.gz
	

foreach(file = dir("EGA_encripted_files", pattern="gpg.md5$", full.names = T), .combine = bind_rows) %do% {
	tibble(`Fastq.File` = !!file, Checksum = read_file(file))
} %>%
	separate(`Fastq.File`, c(sprintf("dummy%s", 1:7), "Sample", "Direction"), remove = F) %>%
	select(-contains("dummy")) %>%
	left_join(
		foreach(file = dir("EGA_encripted_files", pattern="gz.md5", full.names = T), .combine = bind_rows) %do% {
			tibble(`Fastq.File` = !!file, Unencrypted.checksum = read_file(file))
		} %>%
			separate(`Fastq.File`, c(sprintf("dummy%s", 1:7), "Sample", "Direction")) %>%
			select(-contains("dummy"))
	) %>%
	mutate(`Fastq.File` = gsub(".gpg.md5", "", `Fastq.File`, fixed = T)) %>%
	mutate(Direction = ifelse(Direction=="R1", "First", "Second")) %>%
	gather(variable, value, -(Sample:Direction)) %>%
	unite(temp, Direction, variable, sep=".") %>%
	spread(temp, value) %>%
	select(1,3, 2, 4, 6, 5, 7) %>%
	write_csv("EGA_samples_MD5.csv")


annot %>% 
	select(			
		UBR, 
		file,
		PSAResult.x, PathGG1, PathGG2, ClinStageT, PathStageT, gleason_total, `N adapter`, `S adapter`, 
		cohort , sample , `position in plate` , `cell type` , batch , cell_no, `truseq adapter`, age
	) %>%
	rename(subjectId = UBR) %>%
	mutate(
		title = "Dissection of prostate tumour, stroma and immune transcriptional components reveals a key contribution of the microenvironment for disease progression", 
		description = "Analysis of four key cell types (epithelial, fbroblast, myeloid and T cells)",
		caseOrControl = NA,
		gender = "male",
		organismPart = "prostate",
		cellLine = NA,
		region = "prostate",
		phenotype="prostate cancer"
	) %>% 
	separate(`file`, c("dummy1", "dummy2", "dummy3", "dummy4", "dummy5", "dummy6" , "alias"), remove = F)  %>% 
	select(-contains("dummy")) %>%
	select(title, alias, description, subjectId, caseOrControl, gender, organismPart, cellLine, region, phenotype, everything()) %>%
	write_csv("EGA_samples_clinic.csv")

