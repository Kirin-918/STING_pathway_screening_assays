IRF3_plate_results <- read_csv("~/Documents/masters/research_project/drug_screening_assay/IRF3_assay_1.csv")

IRF3_expression <- IRF3_plate_results |> 
  summarise(
    nuc_ratio = mean(`Nuclear translocation ratio`),
    mean_cell = mean(`Total cell count`),
    SEM = sd(`Nuclear translocation ratio`, na.rm = TRUE) / sqrt(n()),
    .by = c(Compound, Plate) #averaging the nuclear trans and cell count and organising by combinations of compound and plate
  ) |> 
  mutate(
    nuc_bycell = nuc_ratio / mean_cell,
    nuc_diff = case_when(
      Plate == 3 ~ (nuc_ratio - 56.33492),
      Plate == 4 ~ (nuc_ratio - 57.46977)
    ),
    percent_inhib = case_when(
      Plate == 3 ~ 100 * (56.33492	- nuc_ratio) / (56.33492 - 54.86843),
      Plate == 4 ~100 * (57.46977	- nuc_ratio) / (57.46977 - 54.92656)
    ), # finding % inhibition for each plate and difference in nuclear translocation
    SEM_percent_inhib = case_when(
      Plate == 3 ~ SEM * 100 / (56.33492 - 54.86843),
      Plate == 4 ~ SEM * 100 / (57.46977 - 54.92656)
    )) |> #applying the same gradient i.e -ve - +ve control, to the SEM. 
  summarise(
    nuc_ratio = mean(nuc_ratio),
    mean_cell = mean(mean_cell),
    mean_nuc_by_cell = mean(nuc_bycell),
    mean_percent_inhib = mean(percent_inhib),
    mean_nuc_diff = mean(nuc_diff),
    SEM_final = mean(SEM),
    SEM_inhib_final = mean(SEM_percent_inhib),
    .by = c(Compound)
  ) #averaging once comparisons to controls have already been done

final_compounds <- IRF3_expression |> #applying filters and removing odd compounds
  mutate(
    compound_of_int =  mean_percent_inhib > 50) |> 
  filter(Compound != c("H151"),
         Compound != c("H151 no diABZI"),
         Compound != c("8"),
         Compound != c("1"),
         Compound != c("DMSO"),
         mean_cell > 200,
         mean_nuc_diff < 3000)


###############################################

############################################### total translocation ratio
ggplot(IRF3_expression) +
  aes(
    x = factor(Compound, levels = unique(Compound)),
    y = nuc_ratio,
    fill = mean_cell
  ) +
  geom_col(width = 0.75) +
  scale_fill_gradient(name = "Mean cell count") +
  labs(
    x = "Compound",
    y = "IRF3 translocation ratio"
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 0)),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(
      size = 16,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 15),
    plot.title = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  coord_cartesian(ylim = c(45, 60))
###############################################

############################################### final compounds 
percentage_IRF3 <- ggplot(final_compounds) +
  aes(x = factor(Compound, levels = unique(Compound)), 
      y = mean_percent_inhib,
      fill = compound_of_int) +
  geom_col(width = 0.75) +
  scale_fill_gradient() +
  geom_hline(yintercept = 100) +
  geom_hline(yintercept = 0) +
  annotate(geom = "text", x = 1.5, y = 120, label = "H151") +
  annotate(geom = "text", x = 1.5, y = 20, label = "DMSO")+
  labs(
    x = "Compound",
    y = "Percentage inhibition of IRF3",
    title = "Percentage inhibition of IRF3",
  ) +
  theme_minimal(base_size = 18) + 
  theme(
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 0)),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(
      size = 16,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 15),
    plot.title = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  ) + 
  guides(fill="none") +
  geom_linerange(
    aes(
      x = Compound,
      ymin = mean_percent_inhib - SEM_inhib_final,
      ymax = mean_percent_inhib + SEM_inhib_final
    ),
    colour = "orange",
    linewidth = 0.7,
    alpha = 0.9
  )
###############################################


library(patchwork)
perccentage_pSTING + percentage_IRF3
