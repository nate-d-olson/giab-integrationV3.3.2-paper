# Generating summary analysis and figures for giab integration pipeline V3.3.2 manuscript 
library(drake)
source("R/packages.R")  # Load all the packages you need.
source("R/functions.R") # Load all the functions into your environment.
source("R/plan.R")      # Build your workflow plan data frame.

# Optionally plot the graph of your workflow.
#config <- drake_config(plan) # nolint
#vis_drake_graph(config)         # nolint

# Now it is time to actually run your project.
make(plan,jobs = 6) # Or make(my_plan, jobs = 2), etc.
