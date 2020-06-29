source("setup.R")
args = commandArgs(trailingOnly=TRUE)
experiment_name = args[1]
expr_index = as.numeric(args[2])
rep_index = as.numeric(args[3])

config = strsplit(experiment_name, "_")[[1]]
#sample pairing define experiment configuration
if (config[1] == "pair") {
  expr_func = experiment_pair_side
} else if (config[1] == "unpair") {
  expr_func = experiment_unpair
} else if (config[1] == "mixed") {
  expr_func = experiment_mixed
}
#set up parameter values
if (config[3] == "shift") {
  para_vary = list(list(name = "dat_type", value = config[2]),
                   list(name = "t", value = expr_index/5),
                   list(name = "C_delta", value = 2))
} else if (config[3] == "incSignal") {
  para_vary = list(list(name = "dat_type", value = config[2]),
                   list(name = "t", value = 2),
                   list(name = "C_delta", value = expr_index/5*2))
}

result = expr_func(para_vary)
save(result, file=paste("results/", experiment_name,"/",
                        expr_index, rep_index,".Rdata", sep = ""))
