# Load required packages ----
library(tidyverse)
library(XLConnect)
library(ggpubr)
library(broom)


# Load Excel file with combined data sets ----

#Set working directory to the folder containing the excel file "Combined data 400 mutants"
setwd("Enter path here") 

#Load data
combined_data <- loadWorkbook("Combined data 400 mutants.xlsx")


# Figure 1 - Expression levels ----
expression_data <- readWorksheet(combined_data, sheet = "Expression levels", check.names = FALSE)

#>>> Plot expression level as a heat map ----
expression_data %>%
  filter(!str_detect(Variant, "control")) %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  group_by(aa112, aa121) %>%
  summarise(norm_mean = mean(Norm), .groups = "drop") %>%
ggplot(aes(x = aa121, y = aa112, fill = norm_mean)) +
  geom_tile() +
  colorspace::scale_fill_continuous_sequential(palette = "YlGn", limits = c(0, NA), 
                       guide = guide_colorbar(title.position = "right",
                                              label.position = "left",
                                              barwidth = 0.5,
                                              barheight = 6,
                                              label.hjust = 1,
                                              draw.ulim = FALSE,
                                              label.theme = element_text(margin = margin(0,-2,0,0, unit = "pt"), size = 7),
                                              title.theme = element_text(margin = margin(0,0,0,3, unit = "pt"), angle = -90, size = 7))) +
  labs(fill = expression("Expression level\nin 然 biotin biding sites/OD"[600]), 
       y = "Amino acid at position 112", 
       x = "Amino acid at position 121") +
  theme_classic() +
  coord_fixed() +
  theme(text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"), 
        axis.title.y = element_text(margin = margin(t = 0, r = 3, b = 0, l = 0, unit = "pt")),
        axis.line = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
        legend.margin = margin(0,0,0,-3, "pt"),
        legend.title.align = 1)

#>>> Standard deviation of expression levels ----
expression_data %>%
  filter(!str_detect(Variant, "control")) %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  group_by(aa112, aa121) %>%
  summarise(norm_sd = sd(Norm), .groups = "drop") %>%
ggplot(aes(x = aa121, y = aa112, fill = norm_sd)) +
  geom_tile() +
  colorspace::scale_fill_continuous_sequential(palette = "YlGn", limits = c(0, 3.17), 
                                               guide = guide_colorbar(title.position = "right",
                                                                      label.position = "left",
                                                                      barwidth = 0.5,
                                                                      barheight = 6,
                                                                      label.hjust = 1,
                                                                      draw.ulim = FALSE,
                                                                      label.theme = element_text(margin = margin(0,-2,0,0, unit = "pt"), size = 7),
                                                                      title.theme = element_text(margin = margin(0,0,0,3, unit = "pt"), angle = -90, size = 7))) +
  labs(fill = expression("Expression level\nin 然 biotin biding sites/OD"[600]), 
       y = "Amino acid at position 112", 
       x = "Amino acid at position 121") +
  theme_classic() +
  coord_fixed() +
  theme(text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"), 
        axis.title.y = element_text(margin = margin(t = 0, r = 3, b = 0, l = 0, unit = "pt")),
        axis.line = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
        legend.margin = margin(0,0,0,-3, "pt"),
        legend.title.align = 1)


#>>> Percentage of mutants with a certain expression level ----
binding_site_threshold <- function(x){
#Calculates the percentage of Sav variants with more than x 然 biotin-binding sites at the cofactor incubation step

  expression_data %>%
  filter(!str_detect(Variant, "control")) %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  mutate(BindingSites = OD*Norm*2.44) %>%
  group_by(aa112, aa121) %>%
  summarise(Sites_mean = mean(BindingSites), .groups = "drop") %>%
  filter(Sites_mean > x) %>%
  summarise(percent = n()/400*100)
}

#Calculate the percentage of Sav variants with an excess of binding sites for cofactor concentrations from 0 to 50 然, assuming 50 % uptake into the periplasm
threshold <- data.frame(concentration = seq(0, 50, 1)) %>%
  mutate(percentage = unlist(sapply(0.5*concentration, binding_site_threshold, USE.NAMES = FALSE)))

#Plot percentage as a function of the periplasmic cofactor concentration
threshold %>%
  ggplot(aes(concentration, percentage)) +
  geom_line(col = "darkgreen", size = 1) +
  geom_area(fill = "darkgreen", alpha = 0.4) +
  geom_segment(aes(x = 10, y = 0, xend = 10, yend = 91.75), col = "black", size = 0.8, lineend = "square", linetype = "dashed") +
  geom_segment(aes(x = 0, y = 91.75, xend = 10, yend = 91.75), col = "black", size = 0.8, linetype = "dashed") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 50), breaks = seq(0, 50, 1), 
                     labels = c("0", rep("", 4), "5", rep("", 4), "10", rep("", 4), "15", rep("", 4), "20", rep("", 4), "25", rep("", 4), "30", rep("", 4), "35", rep("", 4), "40", rep("", 4), "45", rep("", 4), "50")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = seq(0, 100, 10),
                     labels = c("0", "", "20", "", "40", "", "60", "", "80", "", "100")) +
  labs(x = "Cofactor concentration in 然", y = "Percentage of Sav variants\nwith free biotin binding sites") +
  theme_classic() +
  theme(text = element_text(size = 7),
        axis.text = element_text(size = 7, colour = "black"),
        axis.title.y = element_text(margin = margin(0,3,0,0,"pt"), size = 7),
        axis.title.x = element_text(margin = margin(3,0,1,0,"pt")),
        axis.ticks = element_line(colour = "black"))


# Process screening results for plotting ----

#>>> Metathesis ----
#Relative cell-specific activity of all Sav 112X K121X mutants
meta_A <- readWorksheet(combined_data, sheet = "Screening_metathesis") %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  group_by(aa112, aa121) %>%
  summarise(norm_wt_mean = mean(norm_wt),
            norm_wt_sd = sd(norm_wt), .groups = "drop")

#Validation of hits (periplasmic screening workflow)
meta_B <- readWorksheet(combined_data, sheet = "Validation_metathesis") %>%
  filter(!Variant %in% c("control_SL")) %>%
  group_by(Variant) %>%
  mutate(norm_wt_mean = mean(norm_wt)) %>%
  ungroup() %>%
  mutate(Variant = reorder(Variant, norm_wt_mean)) %>%
  mutate(Variant = recode(Variant, "control_empty" = "Empty\nvector", "control_wt" = "wt"))

#In vitro results
meta_C <- readWorksheet(combined_data, sheet = "invitro") %>%
  gather(Replicate, TON, -c("Reaction", "Variant")) %>%
  select(-c("Replicate")) %>%
  filter(Reaction == "Meta_sulfo") %>%
  mutate(Variant = factor(Variant, levels = c("cat", "wt", "IV", "WV", "II"), labels = c("Cofactor", "wt", "IV", "WV", "II")))


#>>> Deallylation (coumarin) ----
#Relative cell-specific activity of all Sav 112X K121X mutants
coumarin_A <- readWorksheet(combined_data, sheet = "Screening_deallylation_coumarin") %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  group_by(aa112, aa121) %>%
  summarise(norm_wt_mean = mean(norm_wt), 
            norm_wt_sd = sd(norm_wt), .groups = "drop")

#Validation of hits (periplasmic screening workflow)
coumarin_B <- readWorksheet(combined_data, sheet = "Validation_deallyl_coumarin") %>%
  group_by(Variant) %>%
  mutate(norm_wt_mean = mean(norm_wt),
         norm_wt = norm_wt) %>%
  ungroup() %>%
  filter(Variant != "control_MA", Variant != "YS") %>%
  mutate(Variant = reorder(Variant, norm_wt_mean)) %>%
  mutate(Variant = recode(Variant, "control_empty" = "Empty\nvector", "control_wt" = "wt"))

#In vitro results
coumarin_C <- readWorksheet(combined_data, sheet = "invitro") %>%
  gather(Replicate, TON, -c("Reaction", "Variant")) %>%
  select(-c("Replicate")) %>%
  filter(Reaction == "Ru_coumarin") %>%
  mutate(Variant = factor(Variant, levels = c("cat", "wt", "MR", "FR", "MW"), labels = c("Cofactor", "wt", "MR", "FR", "MW")))

#>>> Deallylation (indole) ----
#Relative cell-specific activity of all Sav 112X K121X mutants
RuIndole_A <- readWorksheet(combined_data, sheet = "Screening_deallylation_indole") %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  group_by(aa112, aa121) %>%
  summarise(norm_wt_mean = mean(norm_wt), 
            norm_wt_sd = sd(norm_wt), .groups = "drop")

#Validation of hits (periplasmic screening workflow)
RuIndole_B <- readWorksheet(combined_data, sheet = "Validation_deallyl_indole") %>%
  group_by(Variant) %>%
  mutate(norm_wt_mean = mean(norm_wt)) %>%
  ungroup() %>%
  filter(Variant != "control_MA", Variant != "YR") %>%
  mutate(Variant = reorder(Variant, norm_wt_mean)) %>%
  mutate(Variant = recode(Variant, "control_empty" = "Empty\nvector", "control_wt" = "wt"))

#In vitro results
RuIndole_C <- readWorksheet(combined_data, sheet = "invitro") %>%
  gather(Replicate, TON, -c("Reaction", "Variant")) %>%
  select(-c("Replicate")) %>%
  filter(Reaction == "Ru_indole") %>%
  mutate(Variant = factor(Variant, levels = c("cat", "wt", "MA", "LR", "MR"), labels = c("Cofactor", "wt", "MA", "LR", "MR")))

#>>> Hydroamination ----
#Relative cell-specific activity of all Sav 112X K121X mutants
AuIndole_A <- readWorksheet(combined_data, sheet = "Screening_Hydroamination") %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  group_by(aa112, aa121) %>%
  summarise(norm_wt_mean = mean(norm_wt), 
            norm_wt_sd = sd(norm_wt), .groups = "drop")

#Validation of hits (periplasmic screening workflow)
AuIndole_B <- readWorksheet(combined_data, sheet = "Validation_hydroamination") %>%
  filter(!Variant %in% c("control_SL")) %>%
  group_by(Variant) %>%
  mutate(norm_wt_mean = mean(norm_wt)) %>%
  ungroup() %>%
  mutate(Variant = reorder(Variant, norm_wt_mean)) %>%
  mutate(Variant = recode(Variant, "control_empty" = "Empty\nvector", "control_wt" = "wt"))

#In vitro results
AuIndole_C <- readWorksheet(combined_data, sheet = "invitro") %>%
  gather(Replicate, TON, -c("Reaction", "Variant")) %>%
  select(-c("Replicate")) %>%
  filter(Reaction == "Au_indole", Variant != "SL") %>%
  mutate(Variant = factor(Variant, levels = c("cat", "wt", "LL", "DL", "FQ"), labels = c("Cofactor", "wt", "LL", "DL", "FQ")))

#>>> Hydroarylation ----
#Relative cell-specific activity of all Sav 112X K121X mutants
AuFluo_A <- readWorksheet(combined_data, sheet = "Screening_hydroarylation") %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  group_by(aa112, aa121) %>%
  summarise(norm_wt_mean = mean(norm_wt), 
            norm_wt_sd = sd(norm_wt), .groups = "drop")

#Validation of hits (periplasmic screening workflow)
AuFluo_B <- readWorksheet(combined_data, sheet = "Validation_hydroarylation") %>%
  filter(!Variant %in% c("control_SL")) %>%
  group_by(Variant) %>%
  mutate(norm_wt_mean = mean(norm_wt)) %>%
  ungroup() %>%
  mutate(Variant = reorder(Variant, norm_wt_mean)) %>%
  mutate(Variant = recode(Variant, "control_empty" = "Empty\nvector", "control_wt" = "wt"))

#In vitro results
AuFluo_C <- readWorksheet(combined_data, sheet = "invitro") %>%
  gather(Replicate, TON, -c("Reaction", "Variant")) %>%
  select(-c("Replicate")) %>%
  filter(Reaction == "Au_fluo") %>%
  mutate(Variant = factor(Variant, levels = c("cat", "wt", "LG", "FG", "GA"), labels = c("Cofactor", "wt", "LG", "FG", "GA")))

# Figure 2 - Screening results ----

#>>> Heat maps ----

#Function for plotting screening results as a heat map
screening_figure_v5 <- function(dataA, breaks_p1 = waiver()){
  #Helper function for adding white space to 0s on the y axis (to align axes)
  add2x0 <- function(breaks){
    limits <- as.character(breaks)
    limits[1] <- "  0"
    return(limits)
  }
  
  p1 <- dataA %>%
    ggplot(aes(x = aa121, y = aa112, fill = norm_wt_mean)) +
    geom_tile() +
    scale_fill_viridis_c(option = "viridis", limits = c(0, NA), 
                         guide = guide_colorbar(direction = "horizontal",
                                                title.position = "top",
                                                barwidth = 4,
                                                barheight = 0.25,
                                                label.hjust = 0.5,
                                                draw.ulim = FALSE,
                                                label.theme = element_text(margin = margin(-1,0,0,0, unit = "pt"), size = 6),
                                                title.theme = element_text(angle = 0, size = 6))) +
    labs(fill = "Relative cell-specific activity", 
         y = NULL, 
         x = NULL) +
    coord_fixed() +
    theme_classic() +
    theme(text = element_text(size = 7, family = "sans", face = "bold"), 
          axis.text.x = element_text(margin = margin(-2,0,0,0,"pt"), size = 6, family = "sans", face = "plain", colour = "black"),
          axis.text.y = element_text(margin = margin(0,-2,0,0,"pt"), size = 6, family = "sans", face = "plain", colour = "black"),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(0,0,0,0,"pt"),
          legend.margin = margin(0,0,-10,0, "pt"),
          legend.position = "top",
          legend.title.align = 1)
}

#Plot heat maps for all five reactions
ggarrange(screening_figure_v5(meta_A),
          screening_figure_v5(coumarin_A), 
          screening_figure_v5(RuIndole_A),
          screening_figure_v5(AuIndole_A),
          screening_figure_v5(AuFluo_A),
          ncol = 5, nrow = 1)


#>>> Standard deviation of screening results ----
#Function for plotting standard deviations as a heat map
screening_sd <- function(dataA){
  max <- max(dataA$norm_wt_mean)
  #Helper function for adding white space to 0s on the y axis (to align axes)
  add2x0 <- function(breaks){
    limits <- as.character(breaks)
    limits[1] <- "  0"
    return(limits)
  }
  
  p1 <- dataA %>%
    ggplot(aes(x = aa121, y = aa112, fill = norm_wt_sd)) +
    geom_tile() +
    scale_fill_viridis_c(option = "viridis", limits = c(0, max), 
                         labels = add2x0,
                         guide = guide_colorbar(title.position = "right",
                                                label.position = "left",
                                                barwidth = 0.5,
                                                barheight = 6,
                                                label.hjust = 1,
                                                draw.ulim = FALSE,
                                                label.theme = element_text(margin = margin(0,-2,0,0, unit = "pt"), size = 7),
                                                title.theme = element_text(margin = margin(0,0,0,3, unit = "pt"), angle = -90, size = 7))) +
    labs(fill = "Standard deviation of\nrelative cell-specific activity", 
         y = "Amino acid at position 112", 
         x = "Amino acid at position 121") +
    coord_fixed() +
    theme_classic() +
    theme(text = element_text(size = 7, family = "sans", colour = "black"),
          axis.text = element_text(size = 7, family = "sans", colour = "black"), 
          axis.title.y = element_text(margin = margin(t = 0, r = 3, b = 0, l = 0, unit = "pt")),
          axis.line = element_blank(),
          axis.ticks = element_line(colour = "black"),
          plot.margin = margin(t = 0, r = 0, b = 6, l = 0, unit = "pt"), 
          legend.margin = margin(0,0,0,-8, "pt"),
          legend.title.align = 1)
}

#Plot standard deviations for all five reactions
ggarrange(screening_sd(meta_A),
          screening_sd(coumarin_A),
          screening_sd(RuIndole_A),
          screening_sd(AuIndole_A),
          screening_sd(AuFluo_A),
          labels = c("a", "b", "c", "d", "e"),
          ncol = 3, nrow = 2)


#>>> Screening results as a violin plot ----

data400_orig <- bind_rows("AuIndole" = AuIndole_A, "AuFluo" = AuFluo_A, "RuIndole" = RuIndole_A, "coumarin" = coumarin_A, "metathesis" = meta_A, .id = "reaction") %>%
  select(reaction, aa112,aa121, norm_wt_mean)

data400_orig %>%
  mutate(reaction = factor(reaction, levels = c("metathesis", "coumarin", "RuIndole", "AuIndole", "AuFluo"),
                           labels = c("Metathesis", "Deallylation\n(coumarin)", "Deallylation\n(indole)", "Hydroamination", "Hydroarylation"))) %>%
  ggplot(aes(reaction, norm_wt_mean)) +
  geom_violin(aes(fill = reaction, col = reaction), alpha = 0.8) +
  geom_jitter(data = . %>% group_by(reaction) %>% top_n(10, norm_wt_mean), aes(reaction, norm_wt_mean), 
              width = 0.05, shape = 16, size = 1, alpha = 0.8) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_fill_manual(values = colorspace::qualitative_hcl(palette = "Harmonic", n = 5)) +
  scale_colour_manual(values = colorspace::qualitative_hcl(palette = "Harmonic", n = 5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 15), breaks = seq(0,15,1), labels = c("0", rep("", 4), "5", rep("", 4), "10", rep("", 4), "15")) +
  labs(x = "Reaction", y = "Relative cell-specific activity") +
  theme_classic() +
  theme(text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"), 
        axis.title.y = element_text(margin = margin(t = 0, r = 3, b = 0, l = 0, unit = "pt")),
        axis.ticks = element_line(colour = "black"),
        plot.margin = margin(t = 5, r = 0, b = 0, l = 0, unit = "pt"),
        legend.position = "none")


#>>> Screening verification combined plot ----
tophits <- bind_rows(list("Metathesis" = select(meta_B, Variant, norm_wt),
                     "De-allylation (coumarin)" = select(coumarin_B, Variant, norm_wt),
                     "De-allylation (indole)" = select(RuIndole_B, Variant, norm_wt),
                     "Hydroamination" = select(AuIndole_B, Variant, norm_wt),
                     "Hydroarylation" = select(AuFluo_B, Variant, norm_wt)),
                     .id = "reaction")

plot_tophits <- function(reaction_id, text_or_blank = element_blank(), line_or_blank = element_blank(), n = 1){
#Plots the top three variants along with wild type and empty vector as a bar chart.
#Provide element_text() to text_or_blank and element_line() to line_or_blank to plot the y axis.
#(Meant to be used for various reactions in ggarrange, of which only the first should have an axis)
#Always uses top_variants as the data set.
  
top_variants <- tophits %>%
  filter(reaction == reaction_id) %>%
  group_by(Variant) %>%
  summarise(mean = mean(norm_wt), .groups = "drop") %>%
  top_n(3, mean) %>%
  pull(Variant) %>%
  as.character()

  tophits %>%
    filter(reaction == reaction_id,
           Variant %in% c(top_variants)) %>%
    mutate(Variant = fct_reorder(Variant, norm_wt, .fun = mean)) %>%
  ggplot(aes(x = Variant, y = norm_wt)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge2(preserve = "single"), 
               fill = colorspace::qualitative_hcl(palette = "Harmonic", n = 5)[n], alpha = 0.8) +
  geom_jitter(shape = 16, size = 1, alpha = 0.4, width = 0.1) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.35, size = 0.7) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_fill_viridis_c() +
  scale_y_continuous(limits = c(0, 17), expand = c(0,0), breaks = seq(0, 17, 1),
                     labels = c("0", rep("", 4), "5", rep("", 4), "10", rep("", 4), "15", "", "")) +
  labs(y = "Relative cell-specific activity", x = NULL, colour = "", legend = "", fill = "") +
  theme_classic() + 
  guides(fill = "none") +
  theme(text = element_text(size = 7, family = "sans", colour = "black"),
        axis.text = element_text(size = 7, family = "sans"), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_text(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
        axis.title.y = text_or_blank,
        axis.line.y = line_or_blank,
        axis.text.y = text_or_blank,
        axis.ticks.y = line_or_blank)
}

ggarrange(plot_tophits("Metathesis", text_or_blank = element_text(colour = "black"), line_or_blank = element_line(colour = "black")),
          plot_tophits("De-allylation (coumarin)", n = 2),
          plot_tophits("De-allylation (indole)", n = 3),
          plot_tophits("Hydroamination", n = 4),
          plot_tophits("Hydroarylation", n = 5),
          nrow = 1, widths = c(7, 5, 5, 5, 5))


# Figure 3 - In vitro data ----

#Combine in vitro results
invitro <- bind_rows(meta_C, coumarin_C, RuIndole_C, AuIndole_C, AuFluo_C)

#Function for plotting in vitro results
invitro_plot <- function(reaction_ID, n = 1){  
  add4x0 <- function(breaks){
    #Adds white space to the labels of the left y axis to align the axes
    limits <- as.character(breaks)
    limits[1] <- "    0"
    return(limits)
  }
  add2x0_end <- function(breaks){
    #Adds white space to the labels of the right y axis to align the axes
    limits <- as.character(breaks)
    limits[1] <- "0  "
    return(limits)
  }
  
  #Mean TON of free cofactor (for scaling the secondary axis)
  free_cat <-   invitro %>%
    filter(Reaction == reaction_ID, Variant == "Cofactor") %>%
    summarise(mean = mean(TON)) %>%
    pull()
    
  invitro %>%
    filter(Reaction == reaction_ID) %>%
    mutate(Variant = fct_reorder(Variant, TON, .fun = mean),
           Variant = fct_relevel(Variant, "Cofactor")) %>%
  ggplot(aes(Variant, TON)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge2(preserve = "single"),
               fill = colorspace::qualitative_hcl(palette = "Harmonic", n = 5)[n], alpha = 0.8) +
  geom_jitter(shape = 16, size = 1, alpha = 0.4, width = 0.2) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.25, size = 0.7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, NA), labels = add4x0,
                     sec.axis = sec_axis(~ . / free_cat, 
                                         name = "Fold change\nover free cofactor",
                                         labels = add2x0_end)) +
  labs(y = "TON", x = NULL, colour = "", legend = "", fill = "") +
  theme_classic() + 
  theme(text = element_text(size = 7, family = "sans", colour = "black"),
        axis.text = element_text(size = 7, family = "sans", colour = "black"),
        axis.ticks.x = element_blank(), 
        axis.text.x = element_text(), 
        axis.ticks.y = element_line(colour = "black"),
        plot.margin = margin(t = 0, r = 0, b = 3, l = 0, unit = "pt"), 
        axis.title.y = element_text(margin = margin(0,3,0,3,"pt"), size = 7))
}

#Plot in vitro results for all five reactions
ggarrange(invitro_plot("Meta_sulfo"),
          invitro_plot("Ru_coumarin", n = 2),
          invitro_plot("Ru_indole", n = 3),
          invitro_plot("Au_indole", n = 4),
          invitro_plot("Au_fluo", n = 5),
          nrow = 5, ncol = 1)


# Figure 4 - Data analysis ----

aminoacids <- c("F", "Y", "C", "W", "L", "P", "H", "Q", "R", "I", "M", "T", "N", "K", "S", "V", "A", "D", "E", "G") %>%
  sort() 

reaction_names <- c("AuIndole", "RuIndole", "AuFluo", "coumarin", "metathesis")
reactions <- c("AuIndole", "RuIndole", "AuFluo", "coumarin", "metathesis")


#>>> Clustering ----

#Scale values to have mean zero and standard deviation 1
data400_clust <- bind_cols("aa112" = AuIndole_A$aa112, 
                           "aa121" = AuIndole_A$aa121, 
                           "AuIndole" = AuIndole_A$norm_wt_mean, 
                           "AuFluo" = AuFluo_A$norm_wt_mean, 
                           "RuIndole" = RuIndole_A$norm_wt_mean, 
                           "coumarin" = coumarin_A$norm_wt_mean, 
                           "metathesis" = meta_A$norm_wt_mean) %>%
  mutate(AuIndole = scale(AuIndole), RuIndole = scale(RuIndole), AuFluo = scale(AuFluo), coumarin = scale(coumarin), metathesis = scale(metathesis)) 

#Create a data frame where each column contains the 20 scaled activity values associated with a specific mutation at position 112 or 121 for a specific reaction
data400_clust_112 <- bind_cols(select(data400_clust, aa112, aa121, AuIndole) %>% spread(aa121, AuIndole, sep = "AuIndole"),
                               select(data400_clust, aa112, aa121, AuFluo) %>% spread(aa121, AuFluo, sep = "AuFluo"),
                               select(data400_clust, aa112, aa121, RuIndole) %>% spread(aa121, RuIndole, sep = "RuIndole"),
                               select(data400_clust, aa112, aa121, coumarin) %>% spread(aa121, coumarin, sep = "coumarin"),
                               select(data400_clust, aa112, aa121, metathesis) %>% spread(aa121, metathesis, sep = "meta")) %>%
  select(-c("aa112...85", "aa112...22", "aa112...43", "aa112...64")) %>%
  rename("aa112" = "aa112...1")

data400_clust_121 <- bind_cols(select(data400_clust, aa112, aa121, AuIndole) %>% spread(aa112, AuIndole, sep = "AuIndole"),
                               select(data400_clust, aa112, aa121, AuFluo) %>% spread(aa112, AuFluo, sep = "AuFluo"),
                               select(data400_clust, aa112, aa121, RuIndole) %>% spread(aa112, RuIndole, sep = "RuIndole"),
                               select(data400_clust, aa112, aa121, coumarin) %>% spread(aa112, coumarin, sep = "coumarin"),
                               select(data400_clust, aa112, aa121, metathesis) %>% spread(aa112, metathesis, sep = "meta")) %>%
  select(-c("aa121...85", "aa121...22", "aa121...43", "aa121...64")) %>%
  rename("aa121" = "aa121...1")

data400_clust_all <- bind_cols(data400_clust_112, data400_clust_121) %>%
  select(-aa112) %>%
  rename(aa = aa121)

#Clustering

#Calculate distance matrix
clust_all <- data400_clust_all %>%
  column_to_rownames("aa") %>%
  dist()

#Hierarchical clustering with complete linkage
tree_all <- hclust(clust_all, method = "complete")

#Plot dendrogram
plot(tree_all, ylim = c(0, 26))


#>>> Amino acid main effects ----

#Calculate the average effect of an amino acid on activity

#Empty data frame for results
main_effects <- data.frame(aa = aminoacids, 
                           AuIndole_112 = 0, 
                           AuFluo_112 = 0, 
                           RuIndole_112 = 0,
                           coumarin_112 = 0, 
                           metathesis_112 = 0,
                           AuIndole_121 = 0, 
                           AuFluo_121 = 0, 
                           RuIndole_121 = 0, 
                           coumarin_121 = 0, 
                           metathesis_121 = 0) %>%
  gather(reaction_pos, effect, -1) %>%
  separate(reaction_pos, into = c("reaction", "position"), sep = "_")

for (reaction_ID in 1:5) { #Loop through reactions
  for (aa_ID in aminoacids) { #Loop through amino acids
    
    #Average scaled activity of all variants harbouring the respective amino acid at position 112
    mean_112 <- data400_clust %>%
      pivot_longer(names_to = "reaction", values_to = "value", cols = -c(1:2)) %>%
      filter(aa112 == aa_ID, reaction == reaction_names[reaction_ID]) %>%
      summarise(mean = mean(value)) %>%
      as.numeric()
    
    #Average scaled activity of all variants harbouring the respective amino acid at position 121
    mean_121 <- data400_clust %>%
      pivot_longer(names_to = "reaction", values_to = "value", cols = -c(1:2)) %>%
      filter(aa121 == aa_ID, reaction == reaction_names[reaction_ID]) %>%
      summarise(mean = mean(value)) %>%
      as.numeric()
    
    #Enter results into data frame
    main_effects <- main_effects %>%
      mutate(effect = ifelse((aa == aa_ID & reaction == reaction_names[reaction_ID] & position == "112"), mean_112, effect)) 
    
    main_effects <- main_effects %>%
      mutate(effect = ifelse((aa == aa_ID & reaction == reaction_names[reaction_ID] & position == "121"), mean_121, effect)) 
  }}


#Amino acid effects as spider plots

#Colours for amino acid labels
aa_colour <- colorspace::qualitative_hcl(palette = "Dark 3", 6)

aa_sorted <- c("A", "V", "I", "L", "M", "F", "Y", "W", "S", "T", "N", "Q", "C", "R", "H", "K", "D", "E", "P", "G")

main_effects %>%
  mutate(reaction = factor(reaction, levels = c("metathesis", "coumarin", "RuIndole", "AuIndole", "AuFluo"),
                           labels = c("Metathesis", "Deallylation (coumarin)", "Deallylation (indole)", "Hydroamination", "Hydroarylation")),
         position = factor(position, labels = c("Position 112", "Position 121"))) %>%
  mutate(aa = fct_relevel(aa, aa_sorted),
         aa_class = fct_collapse(aa, "Hydrophobic aliphatic" = c("A", "V", "I", "L", "M"), 
                                 "Hydrophobic aromatic" = c("Y", "F", "W"), 
                                 "Polar neutral" = c("C", "N", "Q", "T", "S"), 
                                 "Charged basic" = c("K", "R", "H"), 
                                 "Charged acidic" = c("D", "E"),
                                 "Unique" = c("G", "P"))) %>%
  arrange(aa) %>%
  ggplot(aes(aa, effect, col = position, group = position)) +
  geom_hline(yintercept = 0, col = "grey50") +
  geom_polygon(fill = NA, size = 1.0) +
  scale_y_continuous(limits = c(-1.5, 2), expand = c(0, 0), minor_breaks = NULL) +
  scale_color_manual(values = c("#4167B0", "#663695")) +
  labs(x = NULL, y = NULL, col = "Position") +
  coord_polar() +
  facet_wrap(~ reaction, ncol = 2) +
  theme_bw() +
  theme(text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = c(rep(aa_colour[1], 5), rep(aa_colour[2], 3), rep(aa_colour[3], 5), rep(aa_colour[4], 3), rep(aa_colour[5], 2), rep(aa_colour[6], 2))),
        strip.text = element_text(size = 7, colour = "black"),
        legend.position = "none",
        panel.grid = element_line())


#>>> Comparison of screening strategies ----

#Simulate which variants would have been found when testing single mutants, combining single mutations and using ISM

top_ISM <- data.frame(reaction = reaction_names, 
                      top112_121 = 0, #ISM - randomize position 112 first, then position 121
                      top121_112 = 0, #ISM - randomize position 121 first, then position 112
                      top400 = 0, #Best mutant of all 400
                      top112 = 0, #Single mutants (position 112)
                      top121 = 0, #Single mutants (position 121)
                      top112plus112 = 0) #Combine the amino acid exchanges found in the best single mutants

for (i in 1:5) {
  #ISM - 121 first
  #Keep position 112 constant at first and screen position 121
  top121 <- data400_orig %>%
    group_by(reaction, aa112) %>%
    filter(reaction == reactions[i]) %>%
    filter(aa112 == "S") %>% #Wild-type amino acid at position 112
    top_n(1, norm_wt_mean)
  
  #Keep best mutation at position 121 and screen position 112
  top_ISM[i,3] <- data400_orig %>%
    group_by(reaction, aa121) %>%
    filter(reaction == reactions[i]) %>%
    filter(aa121 == top121$aa121) %>%
    top_n(1, norm_wt_mean) %>%
    ungroup() %>%
    select(norm_wt_mean)
  
  
  #ISM - 112 first
  #Keep position 121 constant at first and screen position 112
  top112 <- data400_orig %>%
    group_by(reaction, aa121) %>%
    filter(reaction == reactions[i]) %>%
    filter(aa121 == "K") %>%
    top_n(1, norm_wt_mean)
  
  #Keep best mutation at position 112 and screen position 121
  top_ISM[i,2] <- data400_orig %>%
    group_by(reaction, aa112) %>%
    filter(reaction == reactions[i]) %>%
    filter(aa112 == top112$aa112) %>%
    top_n(1, norm_wt_mean) %>%
    ungroup() %>%
    select(norm_wt_mean)
  
  
  #Best variant of all 400
  top_ISM[i,4] <- data400_orig %>%
    group_by(reaction) %>%
    filter(reaction == reactions[i]) %>%
    top_n(1, norm_wt_mean) %>%
    ungroup() %>%
    select(norm_wt_mean)
  
  
  #Screening only single mutants - position 112
  top_ISM[i,5] <- top112 %>%
    ungroup() %>%
    select(norm_wt_mean)
  
  #Screening only single mutants - position 121
  top_ISM[i,6] <- top121 %>%
    ungroup() %>%
    select(norm_wt_mean)
  
  #Combine the best single mutations
  top_ISM[i,7] <- data400_orig %>%
    group_by(reaction, aa112, aa121) %>%
    filter(reaction == reactions[i]) %>%
    filter(aa112 == top112$aa112, aa121 == top121$aa121) %>%
    ungroup() %>%
    select(norm_wt_mean)
}

#For ISM, two directions are possible (112 -> 121 or 121 -> 112)
#Typically, one tests all positions of interest individually first and then proceeds from the position 
#that led to the highest activity to positions resulting in lower activity
top_ISM <- top_ISM %>%
  mutate(top_ISM = ifelse(top112 > top121, top112_121, top121_112))


#Simulate which mutants would have been found using reduced amino acid sets

NDT <- str_split(c("RNDCGHILFSYV"), pattern = "") %>%
  unlist()

NRT <- str_split(c("RNDCGHSY"), pattern = "") %>%
  unlist()

top_reduced_AA <- data.frame(reaction = reactions,
                             top_NDT = 0,
                             top_NRT = 0)

for (i in 1:5) {
  #NDT
  top_reduced_AA[i, "top_NDT"] <- data400_orig %>%
    group_by(reaction, aa112, aa121) %>%
    filter(reaction == reactions[i]) %>%
    filter(aa112 %in% NDT, aa121 %in% NDT) %>%
    ungroup() %>%
    top_n(1, norm_wt_mean) %>%
    select(norm_wt_mean)

  
  #NRT
  top_reduced_AA[i, "top_NRT"] <- data400_orig %>%
    group_by(reaction, aa112, aa121) %>%
    filter(reaction == reactions[i]) %>%
    filter(aa112 %in% NDT, aa121 %in% NRT) %>%
    ungroup() %>%
    top_n(1, norm_wt_mean) %>%
    select(norm_wt_mean)
}


#Combine all results and divide the activity by the activity of the best mutant among all 400 variants
top_DE <- left_join(top_ISM, top_reduced_AA, by = c("reaction")) %>%
  mutate(top112_121 = top112_121/top400, 
         top121_112 = top121_112/top400, 
         top112 = top112/top400, 
         top121 = top121/top400, 
         top112plus112 = top112plus112/top400, 
         top_ISM = top_ISM/top400,
         top_NDT = top_NDT/top400, 
         top_NRT = top_NRT/top400,
         top400 = 1)

#Plot results
top_DE %>%
  gather(method, relTON, -reaction) %>%
  mutate(strategy = method) %>%
  mutate(strategy = recode(strategy, top112_121 = "ISM both ways", top121_112 = "ISM both ways", top121 = "One position", top112 = "One position", top112plus112 = "Combining mutations", top_ISM = "ISM", top400 = "All 400", top_NDT = "Reduced AA set (NDT)", top_DBK = "Reduced AA set (DBK)", top_NRT = "Reduced AA set (NRT)")) %>%
  filter(strategy != "ISM both ways", strategy != "Reduced AA set (DBK)") %>%
  mutate(strategy = factor(strategy, levels = c("All 400", "One position", "Combining mutations", "ISM", "Reduced AA set (NDT)", "Reduced AA set (DBK)", "Reduced AA set (NRT)"),
                           labels = c("Full-factorial\nlibrary", "Single\nmutants", "Combining\nmutations", "ISM", "Reduced AA\nset (NDT)", "Reduced AA\nset (DBK)", "Reduced AA\nset (NRT)")),
         reaction = factor(reaction, levels = c("metathesis", "coumarin", "RuIndole", "AuIndole", "AuFluo"),
                           labels = c("Metathesis", "Deallylation (coumarin)", "Deallylation (indole)", "Hydroamination", "Hydroarylation"))) %>%
  ggplot(aes(strategy, relTON)) +
  stat_summary(fun = "mean", geom = "bar", alpha = 0.3) +
  geom_jitter(aes(strategy, relTON, col = reaction), size = 2, shape = 16, width = 0.1, alpha = 0.8) +
  scale_colour_manual(values = colorspace::qualitative_hcl(palette = "Harmonic", n = 5)) +
  theme_classic() +
  labs(y = "Relative activity compared to best 112X 121X mutant", x = "Search strategy", col = "Reaction:") +
  theme(text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"),
        axis.title.y = element_text(margin = margin(t = 0, r = 4, b = 0, l = 0, unit = "pt")),
        axis.title.x = element_text(margin = margin(t = 6, r = 0, b = 0, l = 0, unit = "pt")),
        legend.text = element_text(size = 7, face = "plain"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour = "black"), 
        plot.margin = margin(0,0,0,0,"pt"),
        legend.position = "bottom",
        legend.margin = margin(-5,0,0,0,"pt")) +
  scale_y_continuous(limits = c(0,1.1), expand = c(0,0)) 



# Supplementary figures and further analyses ----
#>>> Activity vs. expression level ----

#Summarise expression data
expression_data2 <- expression_data %>%
  filter(!str_detect(Variant, "control")) %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  group_by(aa112, aa121) %>%
  summarise(norm_mean = mean(Norm), .groups = "drop")

#Combine expression and activity data
expression_activity_list <- list(expr = expression_data2 %>% dplyr::rename("Expression" = norm_mean),
                                 meta = select(meta_A, aa112, aa121, norm_wt_mean) %>% dplyr::rename("Metathesis" = norm_wt_mean),
                                 coumarin = select(coumarin_A, aa112, aa121, norm_wt_mean) %>% dplyr::rename("Deallylation (coumarin)" = norm_wt_mean),
                                 RuIndole = select(RuIndole_A, aa112, aa121, norm_wt_mean) %>% dplyr::rename("Deallylation (indole)" = norm_wt_mean),
                                 AuIndole = select(AuIndole_A, aa112, aa121, norm_wt_mean) %>% dplyr::rename("Hydroamination" = norm_wt_mean),
                                 AuFluo = select(AuFluo_A, aa112, aa121, norm_wt_mean) %>% dplyr::rename("Hydroarylation" = norm_wt_mean))

expression_activity <- purrr::reduce(expression_activity_list, left_join, by = c("aa112", "aa121")) %>%
  pivot_longer(-c(1:3), names_to = "Reaction", values_to = "Activity") %>%
  mutate(Reaction = fct_relevel(Reaction, "Metathesis"),
         Reaction = fct_relevel(Reaction, "Hydroarylation", after = Inf))


#Plot activity vs. expression level (Supplementary Fig. 7)
ggplot(expression_activity, aes(Expression, Activity)) +
  geom_point(alpha = 0.6, shape = 16, size = 1.5) +
  facet_wrap(~ Reaction, scales = "free") +
  labs(x = expression("Expression level in 然 biotin binding sites/OD"[600]), y = "Relative cell-specific activity") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  theme_bw() +
  theme(text = element_text(size = 7, colour = "black", family = "sans"),
        axis.text = element_text(size = 7, colour = "black", family = "sans"),
        strip.text = element_text(size = 7, colour = "black"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(margin = margin(0,5,0,0, unit = "pt")),
        axis.title.x = element_text(margin = margin(3,0,0,0, unit = "pt")),
        plot.margin = margin(0,0,0,0,"pt"))


#Linear regression to check correlation between activity and expression level
expression_activity_list_lm <- list(expr = expression_data2 %>% rename("Expression" = norm_mean),
                                 meta = select(meta_A, aa112, aa121, norm_wt_mean) %>% rename("Metathesis" = norm_wt_mean),
                                 coumarin = select(coumarin_A, aa112, aa121, norm_wt_mean) %>% rename("DeallylationC" = norm_wt_mean),
                                 RuIndole = select(RuIndole_A, aa112, aa121, norm_wt_mean) %>% rename("DeallylationI" = norm_wt_mean),
                                 AuIndole = select(AuIndole_A, aa112, aa121, norm_wt_mean) %>% rename("Hydroamination" = norm_wt_mean),
                                 AuFluo = select(AuFluo_A, aa112, aa121, norm_wt_mean) %>% rename("Hydroarylation" = norm_wt_mean))

reduce(expression_activity_list_lm, left_join, by = c("aa112", "aa121")) %>%
  lm(data = ., Metathesis ~ Expression) %>%
  summary()

reduce(expression_activity_list_lm, left_join, by = c("aa112", "aa121")) %>%
  lm(data = ., DeallylationC ~ Expression) %>%
  summary()

reduce(expression_activity_list_lm, left_join, by = c("aa112", "aa121")) %>%
  lm(data = ., DeallylationI ~ Expression) %>%
  summary()

reduce(expression_activity_list_lm, left_join, by = c("aa112", "aa121")) %>%
  lm(data = ., Hydroamination ~ Expression) %>%
  summary()

reduce(expression_activity_list_lm, left_join, by = c("aa112", "aa121")) %>%
  lm(data = ., Hydroarylation ~ Expression) %>%
  summary()

#>>> Analysis of variance ----

#Metathesis
#Data set including replicates
lm_meta_res <- readWorksheet(combined_data, sheet = "Screening_metathesis") %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  select(aa112, aa121, norm_wt)

anova_meta <- aov(norm_wt ~ aa112 * aa121, data = lm_meta_res)
summary.lm(anova_meta)
summary(anova_meta)
round(summary(anova_meta)[[1]][2]/sum(summary(anova_meta)[[1]][2]), 2) #Variance explained

#Deallylation (coumarin)
#Data set including replicates
lm_coumarin_res <- readWorksheet(combined_data, sheet = "Screening_deallylation_coumarin") %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  select(aa112, aa121, norm_wt)

anova_coumarin <- aov(norm_wt ~ aa112 * aa121, data = lm_coumarin_res)
summary.lm(anova_coumarin)
summary(anova_coumarin)
round(summary(anova_coumarin)[[1]][2]/sum(summary(anova_coumarin)[[1]][2]), 2) #Variance explained

#Deallylation (indole)
#Data set including replicates
lm_RuIndole_res <- readWorksheet(combined_data, sheet = "Screening_deallylation_indole") %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  select(aa112, aa121, norm_wt)

anova_RuIndole <- aov(norm_wt ~ aa112 * aa121, data = lm_RuIndole_res)
summary.lm(anova_RuIndole)
summary(anova_RuIndole)
round(summary(anova_RuIndole)[[1]][2]/sum(summary(anova_RuIndole)[[1]][2]), 2) #Variance explained

#Hydroamination
#Data set including replicates
lm_AuIndole_res <- readWorksheet(combined_data, sheet = "Screening_Hydroamination") %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  select(aa112, aa121, norm_wt)

anova_AuIndole <- aov(norm_wt ~ aa112 * aa121, data = lm_AuIndole_res)
summary.lm(anova_AuIndole)
summary(anova_AuIndole)
round(summary(anova_AuIndole)[[1]][2]/sum(summary(anova_AuIndole)[[1]][2]), 2) #Variance explained

#Hydroarylation
#Data set including replicates
lm_AuFluo_res <- readWorksheet(combined_data, sheet = "Screening_hydroarylation")  %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  select(aa112, aa121, norm_wt)

anova_AuFluo <- aov(norm_wt ~ aa112 * aa121, data = lm_AuFluo_res)
summary.lm(anova_AuFluo)
summary(anova_AuFluo)
round(summary(anova_AuFluo)[[1]][2]/sum(summary(anova_AuFluo)[[1]][2]), 2) #Variance explained


#>>> Predicted activity based on additive model vs. measured activity ----

#Predict activity based on a purely additive model (no interaction between aa112 and aa121)
ANOVA_pred <- augment(lm(norm_wt ~ aa112 + aa121, data = lm_RuIndole_res), newdata = RuIndole_A) %>%
  dplyr::rename("prediction" = ".fitted") %>%
  mutate(reaction = "Deallylation (indole)") %>%
  bind_rows(augment(lm(norm_wt ~ aa112 + aa121, data = lm_AuIndole_res), newdata = AuIndole_A) %>%
              dplyr::rename("prediction" = ".fitted") %>%
              mutate(reaction = "Hydroamination")) %>%
  bind_rows(augment(lm(norm_wt ~ aa112 + aa121, data = lm_coumarin_res), newdata = coumarin_A) %>%
              dplyr::rename("prediction" = ".fitted") %>%
              mutate(reaction = "Deallylation (coumarin)")) %>%
  bind_rows(augment(lm(norm_wt ~ aa112 + aa121, data = lm_AuFluo_res), newdata = AuFluo_A) %>%
              dplyr::rename("prediction" = ".fitted") %>%
              mutate(reaction = "Hydroarylation")) %>%
  bind_rows(augment(lm(norm_wt ~ aa112 + aa121, data = lm_meta_res), newdata = meta_A) %>%
              dplyr::rename("prediction" = ".fitted") %>%
              mutate(reaction = "Metathesis")) %>%
  group_by(reaction) %>%
  mutate(limit = max(norm_wt_mean, prediction)) %>%
  ungroup()

#Limits for making the plots have identical axes
dummy_limits <- ANOVA_pred %>%
  select(reaction, limit, norm_wt_mean = limit, prediction = limit) %>%
  mutate(reaction = fct_relevel(reaction, "Metathesis"),
         reaction = fct_relevel(reaction, "Hydroarylation", after = Inf))

#Plot predictions vs. measured values (Supplementary Fig. 9)
ANOVA_pred %>%
  mutate(reaction = fct_relevel(reaction, "Metathesis"),
         reaction = fct_relevel(reaction, "Hydroarylation", after = Inf)) %>%
ggplot(aes(norm_wt_mean, prediction)) +
  geom_abline(intercept = 0, slope = 1, size = 1.25, col = "grey70") +
  geom_point(alpha = 0.6, shape = 16, size = 1.5) +
  facet_wrap(~ reaction, scales = "free") +
  geom_blank(data = dummy_limits) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0,0.1))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0,0.1))) + 
  labs(x = "Measured relative cell-specific activity", y = "Predicted activity based on additive model") +
  theme_bw() +
  theme(text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"),
        strip.text = element_text(size = 7, colour = "black"),
        strip.background = element_rect(fill = "white"),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
        axis.title.y = element_text(margin = margin(0,5,0,0,"pt")),
        axis.title.x = element_text(margin = margin(3,0,0,0,"pt")),
        panel.grid = element_blank())
