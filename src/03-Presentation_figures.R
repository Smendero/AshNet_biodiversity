
ggplot(unique(Sample_metadata[!Sample_metadata$Treatment == "C",c("Site", "ash_amt")]), aes(Site, ash_amt)) +
        geom_point(size = 5) + 
        theme_minimal()+
        ylab(expression(kg~ha^-1~ash))


ggplot(unique(Sample_metadata[!Sample_metadata$Treatment == "C",c("Site", "ash_amt", "Ca_application_amt", "Ash.Na.kg.applied")]), aes(Ca_application_amt, Ash.Na.kg.applied, color = Site, shape = Site)) +
        geom_point(size = 10) + 
        theme_minimal()+
        ylab(expression(kg~ha^-1~ash~Na)) +
        xlab(expression(kg~ha^-1~ash~Ca)) +
        scale_shape_manual(values = c(1,2,3,4,5,6,7,8)) +
        geom_text(aes(label = ifelse(Ca_application_amt>500|Ash.Na.kg.applied>50, ash_amt, "")), hjust = -1, size = 10) +
        coord_cartesian(xlim = c(0, 3750)) + 
        theme(axis.title = element_text(size = 20))

metadata_char_cols <- colnames(metadata)[sapply(metadata, class) == "character" & colnames(metadata) != "Site_name"]

metadata[1:8, ] %>%
        dplyr::select(c("Site_name", metadata_char_cols)) %>%
        pivot_longer(cols = !c("Site_name"), names_to = "Parameter", values_to = "value") %>%
        group_by(Parameter) %>%
        summarize(groups = length(unique(value))) %>%
        filter(Parameter %in% c("Ash.feedstock", "Dominant.tree.species", "Plot.size", "Soil.parent.material")) %>%
        ggplot(aes(Parameter, groups)) + 
        geom_bar(stat = "identity", fill = cust.cols[2]) + 
        theme_minimal() +
        theme(axis.text = element_text(size = 15))

range(2017 - as.numeric(metadata$Year.ash.applied))

metadata$Site_name[metadata$Year.ash.applied %in% c("2011", "2012", "2013")]

