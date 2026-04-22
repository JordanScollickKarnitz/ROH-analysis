🦌 Welcome to the Elk ROH Analysis Toolkit
This repository contains everything you need to complete the Runs of Homozygosity (ROH) lab module for our course. No prior experience with genomics code is required — just follow the steps below in order.

📖 Start Here — Read This First
Before opening any code, download and read the instructor guide:

📄 README_teacher.docx
This Word document explains everything — what ROH are, what elk populations we're studying, what every file does, and how to interpret your results. Read this before you touch any code.


🗂️ File Guide — What to Open and When
StepFile to OpenWhat it Does1README_teacher.docxRead the background and full instructions2scripts/simulate_elk_data.RGenerates the elk SNP dataset you'll analyze3scripts/run_elk_roh_analysis.RRuns the complete analysis — ROH detection, statistics, and all figures4results/tables/Your output data tables (created after Step 3)5results/figures/Your output figures (created after Step 3)

💡 You do not need to open anything in the R/ folder. Those are the behind-the-scenes functions that run automatically when you source the scripts above. You are welcome to read them — they are heavily commented — but it is not required.


⚡ Quick Start (3 commands in RStudio)
Open RStudio, set your working directory to this folder, then run these three lines:
rsource("scripts/simulate_elk_data.R")
rsource("scripts/run_elk_roh_analysis.R")
That's it. Results will appear in results/tables/ and results/figures/.

📦 Install packages first (one time only)
If this is your first time, run this before anything else:
rinstall.packages(c("ggplot2", "ggridges", "scales", "vcfR"))
