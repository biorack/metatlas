# max_abundance_isotopes.R

# An R script to isolate the isotopes of each element 
#    that have maximum terrestrial abundance

require(ecipex)

nistiso$max_abundance = ave(nistiso$abundance, 
						    nistiso$element, 
						    FUN=max)

keeps = nistiso[nistiso$abundance == nistiso$max_abundance,]

write.table(keeps[, c('element', 'nucleons')], 
		    'max_abundance_isotopes.csv', 
		    row.names=F,
		    col.names=F,
		    sep = ','
		    )