### Produce shadow plot ###
## Load in required packages - install if not
if (!require("circlize")) install.packages("circlize")
library(circlize)

## Set working directory
setwd(<>)

## Read in data
connect_genome <- read.delim('relative_genome_connections.tsv')

## Plot
# Initialise sag file
svg('quasisymmetric_shadow_plot.svg', 
    width = 15, height = 15, bg = 'transparent')

# Set circlize parameters
circos.par(start.degree = 90, gap.degree = 0, track.height = 0.055, track.margin = c(0,0))

# Set track parameters
circos.initialize(sectors = c('1','2'), sector.width = 1, xlim = c(0, 1))

# Draw track
circos.track(ylim = c(0, 0.001), bg.col = 'black', bg.border = 'white')

# Draw all shadow lines with 10% opacity
for (line_i in 1:nrow(connect_genome)){
  line_oi <- connect_genome[line_i,]
  circos.link('1', c(line_oi$relative_replichore_1_start, line_oi$relative_replichore_1_end), 
              '2', c(1-line_oi$relative_replichore_2_start, 1-line_oi$relative_replichore_2_end),
              col = '#0000000D')
}

# Close plotting 
dev.off()
# Clear circlize plot
circos.clear()

