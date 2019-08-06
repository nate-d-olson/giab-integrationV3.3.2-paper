source("R/packages.R")  # Load all the packages you need.
source("R/functions.R") # Load all the functions into your environment.
source("R/plan.R")      # Build your workflow plan data frame.

# Optionally plot the graph of your workflow.
# config <- drake_config(plan) # nolint
# vis_drake_graph(config)         # nolint

# Now it is time to actually run your project.
make(plan, jobs = 12) # Or make(my_plan, jobs = 2), etc.
