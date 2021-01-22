#code edited from https://github.com/noahabe/pulp/blob/master/examples/test1.py

from pulp import *

# A new LP problem
prob = LpProblem("ex1", LpMaximize)
#prob = LpProblem("ex1", LpMinimize)

# Variables
x1 = LpVariable("x1", 0, None)
x2 = LpVariable("x2", 0, None)

# Objective
prob += 2 * x1 + 5 * x2 

# Constraints
prob += 3 * x1 + 2 * x2 <= 6  
prob += -2 * x1 + 4 * x2 <= 8  
prob += x1 + x2 >= 1  

# Write the problem as an LP file
prob.writeLP("test1.lp")

# Solve the problem using the default solver
prob.solve()
# Use prob.solve(GLPK()) instead to choose GLPK as the solver
# Use GLPK(msg = 0) to suppress GLPK messages
# If GLPK is not in your path and you lack the pulpGLPK module,
# replace GLPK() with GLPK("/path/")
# Where /path/ is the path to glpsol (excluding glpsol itself).
# If you want to use CPLEX, use CPLEX() instead of GLPK().
# If you want to use XPRESS, use XPRESS() instead of GLPK().
# If you want to use COIN, use COIN() instead of GLPK(). In this last case,
# two paths may be provided (one to clp, one to cbc).

# Print the status of the solved LP
print("Status:", LpStatus[prob.status])

# Print the value of the variables at the optimum
for v in prob.variables():
	print(v.name, "=", v.varValue)

# Print the value of the objective
print("objective=", value(prob.objective))
