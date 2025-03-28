## ---- random subsets ----

seeds <- list()

seeds$cl_subsets <- rbindlist(list(
  data.table(name = str_c("r", 10, ".", 1:10), 
             seed = c(76051, 27686, 60881, 96010, 21636,  
                      95722, 52190, 32418, 78809, 16460)),
  data.table(name = str_c("r", 20, ".", 1:10), 
             seed = c(55832, 33141, 29411, 47977, 47556, 
                      42260, 49859, 56298, 68241, 26665)),
  data.table(name = str_c("r", 30, ".", 1:10), 
             seed = c(39017, 16539, 23580, 78021, 98185, 
                      40101, 16154, 43546, 35636, 67284)),
  data.table(name = str_c("r", 50, ".", 1:10), 
             seed = c(74730, 37992, 99880, 60579, 55543, 
                      71326, 73079, 58715, 28121, 25770)),
  data.table(name = str_c("r", 100, ".", 1:10), 
             seed = c(26327, 42079, 87450, 69377, 15247, 
                      27114, 76848, 85076, 91155, 61657)),
  data.table(name = str_c("r", 200, ".", 1:10), 
             seed = c(42790, 53560, 76948, 99320, 95753, 
                      74973, 55784, 51426, 43000, 34779)),
  data.table(name = str_c("r", 500, ".", 1:10), 
             seed = c(24022, 78870, 31236, 67505, 39942, 
                      28989, 54639, 52702, 33779, 83893))))

# The random cell line subsampling is coupled to individual seeds

for (.ssq in c(100, 200, 500)) {
  for (.i in 1:10) {
    .n <- str_c("r", .ssq, ".", .i)
    .s <- seeds$cl_subsets[name == .n, seed]
    set.seed(.s)
    cl_subsets$random[[str_c(.n, "_", .s)]] <- sample(cl_subset_default, .ssq)}}
