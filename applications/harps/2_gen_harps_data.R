library(dplyr)
library(readr)

data1 <- read_csv("data/harps/harps_gto_processed_with_cluster_label.csv")


selected_columns1 <- colnames(data1)[grepl("Mg|Al|Si|Ca|TiI|fe|Cu|Zn|Sr|Y|ZrII|Ce|Ba",colnames(data1))]

selected_columns_new1 <- selected_columns1[! selected_columns1 %in% c('Y_tsne_teffcut40', 'Y_tsne_teffcut40_old', 'Y_tsne_teffcut40_errlim_mc', "Y_tsne_teffcut40_errlim_mc_old", "Y_tsne_teffcut40_nofeh_mc", "Y_tsne_teffcut40_mc", "Yg_sig", "vYg", "vYg_sig", "[Fe/H]_C11", "Yg", "nCu", "nZn" , "nSr", "nY", "nZrII", "nBa", "nCe", "nAl", "nMg", "nSi", "nCa", "[Fe/H]_T13")]

selected_columns_new1 <- selected_columns_new1 %>% append("cluster") #"tsne_class_teffcut40",



#
selected_columns_new1 <- selected_columns_new1[selected_columns_new1 %in% c("feh", "AlMg", "TiIFe", "SiFe", "MgFe", "ZnFe", "CuFe", "YMg", "BaFe", "CeFe", "ZrBa", "SrBa", "ZrIIFe", "CaFe", "AlFe", "YFe", "ZrFe", "SrFe", "cluster")] ###Some what correct
#




df <- data1 %>% select(all_of(selected_columns_new1))
write_rds(df, "data/harps/harps_data.rds")
