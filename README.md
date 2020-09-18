# ses-opt-instrument
 SESYNC-RFF Two Instrument Model

### objective.m
  - Objective function as described in manuscript_2020_09_15.tex
  - No discounting - a "two step" model
  - Damages calculated at end of step 2 only, not considering damages incurred in step 1
  
### objective_homogControls.m
  - Same as the above, but no patch-specific controls.
  - X1 = X2, I1 = I2.

### constraints.m
  - Budget constraint
  - Works for either objective function

### main.m
  - Parameters set at the top carry throughout the script
  - Figure 1 - Contour graph showing budget percent to I as a function of p and q
  - Figure 2 - Share of budget to I as a function of budget
  - Figure 3 - Expected damages as a function of budget
  - Figure 4/5 - Don't exist
  - Figure 6 - Share of budget to X and I as a function of budget for many efficacies
  - Figure 7 - Expected damages as a function of budget for many efficacies
  - Figure 8 - Four-panel graph, optimal strategy as a function of budget for 4 parameter spaces of efficacy
