# data wrangling

psting_plate_result <- read_csv("~/Documents/masters/research_project/pSTING_assay_data_final.csv")

psting_expression <- psting_plate_result |> 
  summarise(
    mean_punctate = mean(`Puncta in entire cell - Number of Objects`),
    mean_cell = mean(`Total cell count`),
    SEM = sd(`Puncta in entire cell - Number of Objects`, na.rm = TRUE) / sqrt(n()),
    .by = c(Compound, Plate) #averaging the puncta and cell count and organising by combinations of compound and plate
  ) |> 
  mutate(
    punctate_by_cell = mean_punctate / mean_cell,
    punctate_diff = case_when(
      Plate == 1 ~ (mean_punctate - 2728.5),
      Plate == 2 ~ (mean_punctate - 446.333333)
    ),
    percent_inhib = case_when(
      Plate == 1 ~ 100 * (2728.5	- mean_punctate) / (2728.5 - 843),
      Plate == 2 ~100 * (446.333333	- mean_punctate) / (446.333333 - 183)
    ), # finding % inhibition for each plate and difference in punctate
    SEM_percent_inhib = case_when(
      Plate == 1 ~ SEM * 100 / (2728.5 - 843),
      Plate == 2 ~ SEM * 100 / (446.333333 - 183)
  )) |> 
  summarise(
    mean_punctate = mean(mean_punctate),
    mean_cell = mean(mean_cell),
    mean_punctate_by_cell = mean(punctate_by_cell),
    mean_percent_inhib = mean(percent_inhib),
    mean_punctate_diff = mean(punctate_diff),
    SEM_final = mean(SEM),
    SEM_inhib_final = mean(SEM_percent_inhib),
    .by = c(Compound)
  ) #averaging once comparisons to controls have already been done 

write_csv(psting_expression, "wrangled_data_pSTING.csv")

  final_compounds_pSTING <- psting_expression |> 
  mutate(
    compound_of_int =  mean_percent_inhib > 50) |> 
    filter(Compound != c("H151"),
           Compound != c("H151 no diABZI"),
           Compound != c("8"),
           Compound != c("DMSO"),
           mean_cell > 200,
           mean_punctate_diff < 3000)
     # removing outliers, 
  
  ############################ Plot for total pSTING puncta count 
  ggplot(psting_expression) +
    aes(
      x = factor(Compound, levels = unique(Compound)),
      y = mean_punctate,
      fill = mean_cell
    ) +
    geom_col(width = 0.75) +
    scale_fill_gradient(name = "Mean cell count") +
    labs(
      x = "Compound",
      y = "pSTING puncta count"
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
    guides(fill="none") +
    geom_linerange(
      aes(
        x = Compound,
        ymin = mean_punctate - SEM_final,
        ymax = mean_punctate + SEM_final
      ),
      colour = "orange",
      linewidth = 0.7,
      alpha = 0.9
    )
  ############################
  
  ############################ logged raw data
  
  ggplot(psting_expression) +
    aes(x = factor(Compound, levels = unique(Compound)), 
        y = mean_punctate, fill = mean_cell) +
    geom_col(width = 0.75) +
    scale_fill_gradient() +
    scale_y_continuous(trans = "log10") +
    labs(
      x = "Compound",
      y = "Log10 pSTING puncta count",
      title = "Total pSTING puncta count",
      fill = "Mean cell count"
    ) +
    theme_classic() +
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
        ymin = mean_punctate - SEM_final,
        ymax = mean_punctate + SEM_final
      ),
      colour = "orange",
      linewidth = 0.7,
      alpha = 0.9
    )
  
  ############################

perccentage_pSTING <- ggplot(final_compounds_pSTING) +
  aes(x = factor(Compound, levels = unique(Compound)), 
      y = mean_percent_inhib,
      fill = compound_of_int) +
  geom_col(width = 0.75) +
  scale_fill_gradient() +
  geom_hline(yintercept = 100) +
  geom_hline(yintercept = 0) +
  annotate(geom = "text", x = 1.5, y = 110, label = "H151") +
  annotate(geom = "text", x = 1.5, y = -10, label = "DMSO")+
  labs(
    x = "Compound",
    y = "Percentage inhibition of pSTING",
    title = "Percentage inhibition of pSTING",
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
