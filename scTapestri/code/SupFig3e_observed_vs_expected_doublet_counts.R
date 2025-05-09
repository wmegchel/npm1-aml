###################################################################################
#
# Wout Megchelenbrink
# May 05, 2025
#
# Calculate the number of EXPECTED "self-doublets" (AA or BB) in a 2-sample pooled library
# and compare with the number of OBSERVED "self-doublets"; which is the number of cells with
# a suspiciously "blurred" protein profile.
###################################################################################


calc.expected.doublets <- function(lib, n.cells=17625, frac.A.obs=.413, frac.B.obs=.464, frac.mixed.obs=.1231)
{
  # let's assume that A and B contain perfect singlets
  # and the mixed only doublets, to estimate the fractions per donor
  # In reality, this is more difficult, because A and B contain different 
  # fractions of self doublets
  frac.A <- frac.A.obs + .5*frac.mixed.obs 
  frac.B <- frac.B.obs + .5*frac.mixed.obs 
  
  # AA mixture 
  AA <- frac.A^2 * n.cells
  
  # BB mixture
  BB <-frac.B^2 * n.cells
  
  # AB + BA mixture
  AB <- 2*frac.A * frac.B * n.cells
  
  real.dbl.frac <- (frac.mixed.obs * n.cells) / AB
  
  dbl.estimates <- round(c(AA, BB, AB) * real.dbl.frac)
  names(dbl.estimates) <- c("AA", "BB", "AB")
  
  DT <- data.table(t(dbl.estimates))
  DT[, library:=lib]
  
  
  DT <- DT[, .(library, p1=AA, p2=BB, mixed=AB)]
  
  return (DT)
  
}

# Get the annotated doublets per cell barcode
meta <- fread("processed_data/scTapestri_90K_cell_barcodes_annotated.tsv")

# Re-annotate the (small nr) of alloSCT cells to the correct "sample"
meta[sample.genotype == "P10R alloSCT donor", sample.genotype:="P10R"]
meta[sample.genotype == "P08R alloSCT donor", sample.genotype:="P08R"]

# Compute observed fractions
DT <- merge(meta[, .N, by=.(library, sample.1, sample.2, sample.genotype)], 
            meta[, .N, by=.(library, sample.1, sample.2)],
            by = c("library", "sample.1", "sample.2"),
            suffixes = c("", ".tot"))
DT[, frac:=N/N.tot]

# Annotate whether cell counts belongs to patient1, patient2 or had a mixed genotype
setorder(DT, library, sample.genotype)
DT[, id:=.I]
DT[id %% 3 == 1, sample:="p1" ]
DT[id %% 3 == 2, sample:="p2" ]
DT[id %% 3 == 0, sample:="mixed"]

libs <- dcast.data.table(DT, formula = library + N.tot ~ sample, value.var= "frac")

exp.doublets <- list()
for(i in 1:nrow(libs))
{
  exp.doublets[[i]] <- calc.expected.doublets(lib=libs[i, library], n.cells=libs[i, N.tot], frac.A.obs=libs[i, p1], frac.B.obs=libs[i, p2], frac.mixed.obs=libs[i, mixed])
}

exp.doublets <- rbindlist(exp.doublets)

DT.expected <- rbind(exp.doublets[, .(sample=str_sub(library, end=4), N = p1)],
                      exp.doublets[, .(sample=str_sub(library, start=6), N = p2)])
             
# For some libraries with low expected counts, we actually did not observe doublets
DT.observed <- meta[cell.adt.profile == "doublet", .(.N), by=.(sample=sample.genotype)]

DT <- merge(DT.expected, DT.observed, by="sample", suffixes=c(".exp", ".obs"))


# Correlation = 0.85, p=0.009
DT[, cor.test(N.exp, N.obs)]

# Scatterplot expected vs observed
ggplot(DT, aes(x=N.exp, y=N.obs, color=sample)) +
geom_point() + 
geom_abline(slope = 1, linetype = "dashed") +
coord_equal() + 
scale_x_continuous(limits = c(0, 2200)) +
scale_y_continuous(limits = c(0, 2200)) +
theme_minimal()

# Save as PDF
ggsave("scTapestri/img/SupFig3e_observed_vs_expected_doublet_counts.pdf", width = 4, height = 3)
