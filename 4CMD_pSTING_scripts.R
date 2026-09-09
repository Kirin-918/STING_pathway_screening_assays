psting_plate_result <- read_csv("~/Documents/masters/research_project/drug_screening_assay/pSTING_CMD_9_12_14_23_sum_data.csv")

psting_expression <- psting_plate_result |> 
  summarise(
    mean_punctate = mean(`Puncta in entire cell - Number of Objects`),
    mean_cell = mean(`Total cell count`),
    SEM = sd(`Puncta in entire cell - Number of Objects`, na.rm = TRUE) / sqrt(n()),
    .by = c(Compound, Concentration) #averaging the puncta and cell count and organising by combinations of compound and plate
  ) |> 
esquisser(psting_expression)

########################################
psting_expression %>%
  filter(!(Compound %in% c("DMSO/diABZI", "DMSO/DMSO"))) %>%
  ggplot() +
  aes(x = Concentration, y = mean_punctate, colour = Compound) +
  geom_point() +
  scale_color_hue(direction = 1) +
  scale_x_continuous(trans = "log10") +
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
  labs(
    x = "Concentration",
    y = "pSTING puncta count"
  ) +
  guides(fill="none") +
  geom_linerange(
    aes(
      x = Concentration,
      ymin = mean_punctate - SEM,
      ymax = mean_punctate + SEM
    ),
    colour = "black",
    linewidth = 0.7,
    alpha = 0.9
  ) +
  geom_smooth(
    se = FALSE
  )

psting_expression_ky <- psting_plate_result |> 
  summarise(
    mean_punctate = mean(`Puncta in entire cell - Number of Objects`),
    mean_cell = mean(`Total cell count`),
    SEM = sd(`Puncta in entire cell - Number of Objects`, na.rm = TRUE) / sqrt(n()),
    .by = c(Compound, Concentration) #averaging the puncta and cell count and organising by combinations of compound and plate
  ) |> 
  filter(
    Compound == 12
  )


psting_expression_ky %>%
  filter(!(Compound %in% c("DMSO/diABZI", "DMSO/DMSO"))) %>%
  ggplot() +
  aes(x = Concentration, y = mean_punctate, colour = Compound) +
  geom_point() +
  scale_color_hue(direction = 1) +
  scale_x_continuous(trans = "log10") +
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
  labs(
    x = "Concentration",
    y = "pSTING puncta count"
  ) +
  guides(fill="none") +
  geom_linerange(
    aes(
      x = Concentration,
      ymin = mean_punctate - SEM,
      ymax = mean_punctate + SEM
    ),
    colour = "black",
    linewidth = 0.7,
    alpha = 0.9
  ) +
  geom_smooth(
    se = FALSE
  )