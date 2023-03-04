library(dplyr)
library(tidyr)
library(readr)
library(arules)
library(stringr)

message("Starting script...")
gene_mappings <- read_csv("/specific/elhanan/PROJECTS/ARULES_ALPHA_EM/DATA/gene_mappings.csv")
message("Finished reading gene_mappings.csv")
gene_mappings <- gene_mappings %>%
  filter(!is.na(ko)) %>%
  distinct(bac_name, ko, .keep_all = TRUE) %>%
  select(genome, ko)
message("Filtered out rows with no KO")
genome_per_ko <- gene_mappings %>%
  group_by(ko) %>%
  summarize(n_genome = n())

genome_per_ko_updated <- genome_per_ko %>%
  filter(n_genome < n_distinct(gene_mappings$genome) * 0.7)

gene_mappings <- gene_mappings %>%
  left_join(genome_per_ko_updated, by = c("ko" = "ko")) 

gene_mappings <- gene_mappings %>%
  filter(!is.na(n_genome)) %>%
  select(genome, ko)
message("Removed KO's that are too common (above 70%)")
gene_mappings <- gene_mappings %>%
  mutate(is_true = TRUE) 

gene_mappings <- gene_mappings %>%
  pivot_wider(names_from = "ko", values_from = "is_true", values_fill = FALSE)
gene_mappings <- gene_mappings %>% tibble::remove_rownames() %>%
  tibble::column_to_rownames(var = "genome")
message("Created two dimensional array with genoms as row names and KO's as column names, values T/F")
gene_mappings_transactions <- as(gene_mappings, "transactions")
message("Finished creating transactions")

mine_hyperedge_sets <- function(sup, confd) { 
  #print out support and confidence every time i run to keep track of current iteration
  message("Support: ", sup)
  message("Confidence: ", confd)
  #begin counting how long the function is running
  start_time <- Sys.time()
  #run apriori algorithm with the correct parameters, looking for hyperedge sets
  rules <- apriori(gene_mappings_transactions, parameter = list(support = sup, confidence = confd, "target" = "hyperedgesets", maxtime = 0, maxlen = 4000, ext = F))
  #save original amount of rules before filtering out irrelevant rules
  amount_rules <- length(rules)
  #print it out 
  message("Amount of rules mined: ", amount_rules)
  #if there are no rules, returns basic values as there's no need to waste runtime
  if (amount_rules == 0) { 
    return(list(run_time = -1,
                amount_rules = 0,
                maximal_rules = 0,
                max_rule = 0,
                rules = NULL))
  }
  
  # Convert rules into a list of vectors that i can use to filter out irrelevant rules
  #first turn rules into the dataframe
  rules_df <- DATAFRAME(rules)
  #then pull out only the relevant items, which are the rules themselves
  rules_list <- as.character(pull(rules_df, items))
  #reverse the order, as automatically the rules are sorted in ascending size and I want them to be in descending size 
  rules_list <- rev(rules_list)
  #remove all curly braces and commas
  rules_list <- sapply(rules_list, 
                       function(x) {
                         gsub('\\{|\\}', '', x)
                       }, USE.NAMES = F)
  rules_list <- sapply(rules_list, 
                       function(x){
                         strsplit(x, split = ",")[[1]]
                       }, USE.NAMES = F)
  #save size of the largest rule, which is the length of the first vector in the reverse array
  max_rule <- length(rules_list[[1]])
  #create super_set, which will save all our final rules in the end
  super_set <- list()
  #while loop that runs until there are no more rules in rules_list (have been removed becuase they're either subsets, or supersets and saved somewhere else)
  while (length(rules_list) > 0) {
    #save the first rule as a string, add it to super_set (it must be a superset as it's the largest rule in the list)
    string_i <- rules_list[[1]]
    super_set <- append(super_set, list(string_i))
    #remove the rule from the list of rules
    rules_list[[1]] <- NULL
    #if the last rule has just been removed, break the loop and skip to the end
    if (length(rules_list) == 0) { break }
    #create a mask which is a logical with true where I want to keep the values, and false where they're subsets and need to be removed
    mask <- sapply(rules_list, 
                   function(x){
                     for (i in x) {
                       if (! i %in% string_i) {return(TRUE)}
                     }
                     return(FALSE)
                   })
    rules_list <- rules_list[mask]
  }
  #save all relevant values (runtime, largest rule, final rules etc)
  maximal_rules <- length(super_set)
  end_time <- Sys.time()
  run_time <- end_time - start_time
  message("Overall runtime: ", run_time)
  
  message("Amount of maximal rules: ", maximal_rules)
  message("Largest rule: ", max_rule)
  message("Size of rules_length: ", length(rules_list))
  #change the rules back into the original string format
  rules_final <- lapply(super_set, 
                        function(x){
                          x <- paste0(x, collapse = ",")
                          x <- paste0("{", x, "}")
                          return(x)
                        })
  #filter out all rules that aren't the ones we want out of rule_df
  rules_df <- rules_df %>%
    filter(items %in% rules_final)
  #return all relevant values
  return(list(run_time = run_time,
              amount_rules = amount_rules,
              maximal_rules = maximal_rules,
              max_rule = max_rule,
              rules = rules_df))
}
message("Start running mine_hyperedge_sets() 27 times!")
support <- c(0.6, 0.6, 0.6, 0.55, 0.55, 0.55, 0.5, 0.5, 0.5, 0.45, 0.45, 0.45, 0.4, 0.4, 0.4, 0.35, 0.35, 0.35, 0.3, 0.3, 0.3, 0.25, 0.25, 0.25, 0.2, 0.2, 0.2)
confidence <- c(0.8, 0.85, 0.9, 0.8, 0.85, 0.9, 0.8, 0.85, 0.9, 0.8, 0.85, 0.9, 0.8, 0.85, 0.9, 0.8, 0.85, 0.9, 0.8, 0.85, 0.9, 0.8, 0.85, 0.9, 0.8, 0.85, 0.9)
run_time <- rep(0, length(support))
amount_rules <- rep(0, length(support))
maximal_rules <- rep(0, length(support))
largest_rule <- rep(0, length(support))
#create variable that holds the beginning of the path to which I'll save the file with the rules of each call, and the final results of all calls
output_folder <- "/specific/elhanan/PROJECTS/ARULES_ALPHA_EM/ANALYSIS"
#call mine_hyperedge_sets 27 times in accordance with support and confidence 
for (i in 1:27) {
  message("Run number: ", i, "...")
  #set current values for support and confidence
  sup <- support[i]
  confd <- confidence[i]
  #create the output file name
  output_file <- file.path(output_folder, paste0("rules_s", sup, "_c", confd, ".txt"))
  #run algorithm, save to res
  res <- mine_hyperedge_sets(sup, confd)
  #save the results in the seperate vectors for each result
  run_time[i] <- res$run_time
  amount_rules[i] <- res$amount_rules
  maximal_rules[i] <- res$maximal_rules
  largest_rule[i] <- res$max_rule
  #create the txt file that holds the rules mined
  write_delim(x = res$rules, file = output_file, delim = "\t")
}
#save the final results of all values into one dataframe, and create a txt file holding them
results <- data.frame(support, confidence, run_time, amount_rules, maximal_rules, largest_rule)
write_delim(x = results, file = file.path(output_folder, paste0("final_results.txt")), delim = "\t")
message("Execution finished!")