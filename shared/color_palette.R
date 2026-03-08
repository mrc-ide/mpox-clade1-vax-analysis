
library(scales)

# Color palettes for the age groups and key populations 

age_groups_eq <- c('0-4' = "#4e6d58",'5-14' = "#749e89", '15+'= "#abccbe",'CSW'="#edd7d9",'ASW'= "#e3aba7",'HCW'="#d88782", 'PBS'="#b1615c")
age_groups_sk <- c("0-4"="#41507b", "5-14"= "#7d87b2", "15+" = "#c2cae3","CSW"= "#e9dcd8","ASW"= "#d7b9a7", "HCW"="#c79b82",'PBS'= "#a1755c")

# Color palette for the scenario modelling 

scenarios_eq <- c("#053c29", "#2b614e", "#418979","#2A675A","#A9845B") 
scenarios_sk <- c("#0F292B","#344f50", "#4E717B","#35525B", "#A9845B")
scenarios_expanded <- c("#0F292B", "#006B6B", "#1ABC9C", "#A7D8D8", "#87CEEB", "#B0C4DE", "#CCCCFF") 

cols_region <- MetBrewer::met.brewer("Demuth", 3)
names(cols_region) <- c("equateur", "sudkivu", "both")

# Create a color palette function
generate_palette_eq <- colorRampPalette(scenarios_eq)
generate_palette_sk <- colorRampPalette(scenarios_sk)
generate_palette_exp <- colorRampPalette(scenarios_expanded)
# Generate a palette with 10 colors
scenarios_palette_eq <- generate_palette_eq(10)  
scenarios_palette_sk <- generate_palette_sk(10)   
scenarios_palette_expanded <- generate_palette_exp(10)