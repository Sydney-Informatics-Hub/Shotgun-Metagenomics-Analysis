
########################################################################################
############ PLOT for Taxonomic association of ARGs id in DAA
########################################################################################
ARG_sp <- read.csv("DAA_ARG_Sp.csv")

ggplot(ARG_sp, aes(fill=Genus, y=TPM, x=Gene)) + 
  geom_bar(position = "fill", stat = "identity")+
  labs(x = "ARG",
       y = "Relative Abundance")+
  scale_fill_igv() +
  facet_wrap( ~ Severity_caries)+ 
  scale_y_continuous(labels = scales::percent)+theme_bw()


