Installing BSTZINB
---------------------------------------------------------------------

# Step 1: Devtools

You would require the R package <tt>devtools</tt> to proceed with installation. If you do not have <tt>devtools</tt>, please install <tt>devtools</tt> first as below. If you have <tt>devtools</tt> installed in your machine, proceed to the next step.

<li><tt>install.packages("devtools")</tt></li>

# Step 2: Installing BSTZINB 

With <tt>devtools</tt> on your machine, you can load it and install the <tt>BSTZINB</tt> in the following way:

<li><tt>library(devtools)</tt></li>
<li><tt>devtools::install_github("SumanM47/BSTZINB",dependencies=TRUE, build_vignettes=TRUE)</tt></li>


Using BSTZINB
------------------------------------------------------------------------

Please see the vignette fordifferent functions available in the package and their usage. You can access the vignette by using the <tt>vignette<\tt> function: <tt>vignette("BSTZINB_vignette",package="BSTZINB")</tt>
