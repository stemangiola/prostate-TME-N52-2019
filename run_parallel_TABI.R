library(tidyverse)
library(foreach)
if("package:TABI" %in% search()) detach("package:TABI", unload=TRUE, force=TRUE)
library(TABI)

args = commandArgs(trailingOnly=TRUE)

load("input_parallel_TABI.RData")

# Filter data sets
ct_label = "E"
ct_label = args[1]

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

ct.annot = 
	annot %>% 
	filter(cell_type_formatted==ct_label)%>% arrange(UBR)

ct.ex = ex  %>%
	select(symbol, as.character(ct.annot$file)) %>%
	filter((.) %>% select(-symbol) %>% rowSums() > 0) %>%
	average_duplicated_genes_tibble_spreaded() %>%
	gather(sample, count, -symbol) %>%
	spread(symbol, count) %>%
	select(sample, everything())

# Run TABI_glm
which_chunk = as.numeric(args[2])
genes_to_analyse = 
	(
		2:ncol(ct.ex) %>%
		split(., ceiling(seq_along(.)/( length(.)%/%5 ))) %>%
		{
			if(length((.))>5) {
				temp = (.)
				temp[[5]] = c(temp[[5]], temp[[6]])  
				temp[1:5]
			} else { (.) }
		} 
	) [[which_chunk]]


tabi_res = TABI_glm(
	
	formula = ~ CAPRA_TOTAL + batch,
	
	data = ct.ex[,c(1,genes_to_analyse)] %>%
		#prostate_df = ct.ex %>% 
		left_join(
			ct.annot %>% 
				select(file, CAPRA_TOTAL, batch), by = c("sample" = "file" )
		) %>% 
		select(-sample), 
	
	prior = list( 
		prop_DE =0.05,
		scale_DE = 5, 	
		nu_global = 40,
		slab_df = 40
	),
	
	model=rstan::stan_model("~/PhD/TABI/src/stan_files/DE_sigmoid.stan"),
	control=list(
		adapt_delta=0.9,
		stepsize = 0.1
		#,
		# max_treedepth =13
	)
)

save(tabi_res, file=sprintf("tabi_res_%s_chunk_%s.RData", ct_label, which_chunk))