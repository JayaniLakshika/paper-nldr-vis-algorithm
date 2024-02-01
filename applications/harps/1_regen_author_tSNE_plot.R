library(dplyr)
library(readr)
library(ggplot2)

data <- read_csv("data/harps/harps_data_by_author.csv")

data <- data %>% mutate(vXg_minus = vXg * (-1),
                        YMg = YFe - MgFe,
                        AlMg = AlFe - MgFe,
                        YAl = YFe - AlFe,
                        SrBa = SrFe - BaFe,
                        ZrBa = ZrIIFe - BaFe,
                        CFe = `[C/H]_SA17` - feh,
                        CO = `[C/H]_SA17` - `[O/H]_6158_BdL15`,
                        sqrt_Vxg_Vzg = sqrt(vXg**2. + vZg**2. ),
                        MgY = MgFe - YFe,
                        BaY = BaFe - YFe,
                        Rg_Rc = Rg - Rc,
                        Rc_rb = Rc - rb,
                        sqrt_errY_errMg = sqrt(errY**2. + errMg**2.),
                        sqrt_errAl_errMg = sqrt(errAl**2. + errMg**2.),
                        sqrt_errY_errAl = sqrt(errY**2. + errAl**2.),
                        sqrt_errSr_errBa = sqrt(errSr**2. + errBa**2.),
                        sqrt_errZrII_errBa = sqrt(errZrII**2. + errBa**2.),
                        sqrt_errCH_errfeh = sqrt(`e_[C/H]_SA17`**2. + erfeh**2.),
                        sqrt_errCH_errOH = sqrt(`e_[C/H]_SA17`**2. + `e_[O/H]_6158_BdL15`**2.),
                        sqrt_vXg_sig_vZg_sig = sqrt(vXg_sig**2. + vZg_sig**2.),
                        sqrt_errY_errBa = sqrt(errY**2. + errBa**2.),
                        sqrt_Rg_sig_Rc_sig = sqrt(Rg_sig**2. + Rc_sig**2.),
                        sqrt_rb_err_Rc_sig = sqrt(rb_err**2. + Rc_sig**2.),
                        abs_JR_st_U_JR_st_L = 0.5*abs(JR_st_U-JR_st_L),
                        abs_JZ_st_U_JZ_st_L = 0.5*abs(JZ_st_U-JZ_st_L))

## Assign outer cluster names
data <- data %>%
  mutate(outer_cluster_teffcut40 = ifelse(tsne_class_teffcut40 == "thin", "",
                                          ifelse(tsne_class_teffcut40 == "thick1", "",
                                                 ifelse(tsne_class_teffcut40 == "thick2", "",
                                                        ifelse(tsne_class_teffcut40 == "thick3", "",
                                                               ifelse(tsne_class_teffcut40 == "thick4", "",
                                                                      ifelse(tsne_class_teffcut40 == "mpthin", "",
                                                                             ifelse(tsne_class_teffcut40 == "mpthintrans", "Transition group",
                                                                                    ifelse(tsne_class_teffcut40 == "smr", "",
                                                                                           ifelse(tsne_class_teffcut40 == "t4trans", "",
                                                                                                  ifelse(tsne_class_teffcut40 == "youngthin", "Young local disc",
                                                                                                         ifelse(tsne_class_teffcut40 == "debris1", "",
                                                                                                                ifelse(tsne_class_teffcut40 == "debris2", "",
                                                                                                                       ifelse(tsne_class_teffcut40 == "debris3", "[s/Fe]-enhanced",
                                                                                                                              ifelse(tsne_class_teffcut40 == "debris4", "",
                                                                                                                                     ifelse(tsne_class_teffcut40 == "debris5", "",
                                                                                                                                            ifelse(tsne_class_teffcut40 == "smr2", "",
                                                                                                                                                   ifelse(tsne_class_teffcut40 == "t2trans1", "Debris candidate",
                                                                                                                                                          ifelse(tsne_class_teffcut40 == "highTi", "Extreme-Ti star",
                                                                                                                                                                 ifelse(tsne_class_teffcut40 == "lowMg", "Low-[Mg/Fe] star",
                                                                                                                                                                        ifelse(tsne_class_teffcut40 == "highAlMg?", "High-[Al/Mg] star",

                                                                                                                                                                               "")))))))))))))))))))))

data <- data %>%
  mutate(outer_cluster_teffcut40_mc = ifelse(tsne_class_mc == "thin", "Thin Disc",
                                             ifelse(tsne_class_mc == "thick1", "Thick Disc I",
                                                    ifelse(tsne_class_mc == "thick2", "Thick Disc II",
                                                           ifelse(tsne_class_mc == "thick3", "Inner Disc I",
                                                                  ifelse(tsne_class_mc == "thick4", "Inner Disc II",
                                                                         ifelse(tsne_class_mc == "mpthin", "Outer Thin Disc",
                                                                                ifelse(tsne_class_mc == "smr", "SMR",
                                                                                       ifelse(tsne_class_mc == "t4trans", "Inner Disc III",
                                                                                              ifelse(tsne_class_mc == "debris1", "",
                                                                                                     ifelse(tsne_class_mc == "debris2", "",
                                                                                                            ifelse(tsne_class_mc == "debris3", "Satellite debris",
                                                                                                                   ifelse(tsne_class_mc == "debris4", "",
                                                                                                                          ifelse(tsne_class_mc == "debris5?", "",
                                                                                                                                 ifelse(tsne_class_mc == "t2trans1", "TII/III",
                                                                                                                                        ifelse(tsne_class_mc == "t2trans2", "",
                                                                                                                                               ifelse(tsne_class_mc == "highTi", "Extreme-Ti star",
                                                                                                                                                      ifelse(tsne_class_mc == "thicklow", "Lowest-[Fe/H] star",


                                                                                                                                                             ""))))))))))))))))))


data <- data %>%
  mutate(outer_cluster_teffcut40_errlim_mc = ifelse(tsne_class_errlim_mc == "thin", "Thin Disc",
                                                    ifelse(tsne_class_errlim_mc == "thick1", "Thick Disc",
                                                           ifelse(tsne_class_teffcut40 == "thick2", "",
                                                                  ifelse(tsne_class_teffcut40 == "thick3", "",
                                                                         ifelse(tsne_class_teffcut40 == "thick4", "Metal-rich Thick Disc",
                                                                                ifelse(tsne_class_teffcut40 == "mpthin", "Metal-poor Thin Disc",
                                                                                       ifelse(tsne_class_teffcut40 == "youngthin", "Young Local Disc",
                                                                                              ifelse(tsne_class_teffcut40 == "smr", "SMR",
                                                                                                     ifelse(tsne_class_teffcut40 == "t4trans", "Transition",

                                                                                                            ifelse(tsne_class_teffcut40 == "debris1", "Debris",
                                                                                                                   ifelse(tsne_class_teffcut40 == "debris2", "",
                                                                                                                          ifelse(tsne_class_teffcut40 == "debris3", "",
                                                                                                                                 ifelse(tsne_class_teffcut40 == "debris4", "",
                                                                                                                                        ifelse(tsne_class_teffcut40 == "debris5?", "",
                                                                                                                                               ifelse(tsne_class_teffcut40 == "t2trans1", "",
                                                                                                                                                      ifelse(tsne_class_teffcut40 == "t2trans2", "",
                                                                                                                                                             ifelse(tsne_class_teffcut40 == "highTi", "",
                                                                                                                                                                    ifelse(tsne_class_teffcut40 == "t2trans3", "",
                                                                                                                                                                           ifelse(tsne_class_teffcut40 == "smr2", "",
                                                                                                                                                                                  ifelse(tsne_class_teffcut40 == "lowMg", "",

                                                                                                                                                                                         "")))))))))))))))))))))





data <- data %>%
  mutate(outer_cluster_teffcut40_simple = ifelse(tsne_class_teffcut40_simple == "thin", "Thin Disc",
                                                 ifelse(tsne_class_teffcut40_simple == "thick1", "Thick Disc I",
                                                        ifelse(tsne_class_teffcut40_simple == "thick2", "Thick Disc II",
                                                               ifelse(tsne_class_teffcut40_simple == "thick3", "Thick Disc III",

                                                                      ifelse(tsne_class_teffcut40_simple == "mpthin", "Metal-poor thin disc",

                                                                             ifelse(tsne_class_teffcut40_simple == "smr", "SMR",
                                                                                    ifelse(tsne_class_teffcut40_simple == "t1trans", "Transition I",

                                                                                           ifelse(tsne_class_teffcut40_simple == "debris", "Satellite debris",
                                                                                                  ifelse(tsne_class_teffcut40_simple == "highAlMg", "High-[Al/Mg]",
                                                                                                         ifelse(tsne_class_teffcut40_simple == "t3trans", "Transition III",
                                                                                                                ifelse(tsne_class_teffcut40 == "highTioutlier", "Extreme-Ti star",
                                                                                                                       ifelse(tsne_class_teffcut40_simple == "lowalphaoutlier", "Low-[Mg/Fe] star",


                                                                                                                              "")))))))))))))







### Actual result generation

actual_tsne <- data %>% select("Star","X_tsne_teffcut40", "Y_tsne_teffcut40", "tsne_class_teffcut40", "outer_cluster_teffcut40",
                               "tsne_class_teffcut40_simple", "outer_cluster_teffcut40_simple","X_tsne_teffcut40_errlim_mc", "Y_tsne_teffcut40_errlim_mc",
                               "tsne_class_errlim_mc","outer_cluster_teffcut40_errlim_mc", "X_tsne_teffcut40_nofeh_mc", "Y_tsne_teffcut40_nofeh_mc",
                               "X_tsne_teffcut40_mc", "Y_tsne_teffcut40_mc", "tsne_class_mc", "outer_cluster_teffcut40_mc")

actual_tsne1 <- actual_tsne %>% filter(outer_cluster_teffcut40 != "")#

actual_tsne1$cluster <- actual_tsne1$outer_cluster_teffcut40

actual_tsne2 <- actual_tsne %>% filter(outer_cluster_teffcut40 == "")

actual_tsne3 <- actual_tsne2 %>%  filter(outer_cluster_teffcut40_mc != "")

actual_tsne3 <- actual_tsne3 %>%
  mutate(outer_cluster_teffcut40_mc = ifelse(outer_cluster_teffcut40_mc == "SMR", "Thin Disc",
                                             ifelse(outer_cluster_teffcut40_mc == "Lowest-[Fe/H] star", "Thick Disc I", outer_cluster_teffcut40_mc))) #

actual_tsne3$cluster <- actual_tsne3$outer_cluster_teffcut40_mc


stars <- c(actual_tsne1$Star, actual_tsne3$Star)

actual_tsne4 <- actual_tsne[-which(actual_tsne$Star %in% stars),]

actual_tsne5 <- actual_tsne4 %>% filter(outer_cluster_teffcut40_simple != "") #

actual_tsne5$cluster <- actual_tsne5$outer_cluster_teffcut40_simple

actual_tsne6 <- actual_tsne4 %>% filter(outer_cluster_teffcut40_simple == "")

actual_tsne6$cluster <- rep("Debris candidate", nrow(actual_tsne6)) #

data_with_cluster_label <- bind_rows(actual_tsne1, actual_tsne3, actual_tsne5, actual_tsne6)

data_with_cluster_label <- data_with_cluster_label |>
  rename(c("tSNE1" = "X_tsne_teffcut40", "tSNE2" = "Y_tsne_teffcut40"))

write_rds(data_with_cluster_label, "data/harps/harps_tsne_40.rds")


# tSNE_plot_actual <- data_with_cluster_label %>%
#   ggplot(aes(x = X_tsne_teffcut40, y = Y_tsne_teffcut40, color = cluster)) +
#   geom_point(aes(shape=cluster, color=cluster)) +
#   scale_shape_manual(values=c(16, 5, 5, 13, 8, 15, 16, 15, 10, 17, 24, 25, 17, 17)) + ggtitle("Visualization from author's \n results by tSNE") + xlab("t-SNE X dimension") + ylab("t-SNE Y dimension") + coord_fixed()
