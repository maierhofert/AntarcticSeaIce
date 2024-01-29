# analyze the permutation test results
res_data_red = readRDS("fogt_rebuild_results_lindetrend_permute1.RDS")
res_data_red2 = readRDS("fogt_rebuild_results_lindetrend_permute2.RDS")
res_data_red = rbind(res_data_red, res_data_red2)
res_data_red3 = readRDS("fogt_rebuild_results_lindetrend_permute3.RDS")
res_data_red = rbind(res_data_red, res_data_red3)

# combine results with results from original method
res_data_red$method = "permuted"
# read in original methods results
# res_data_red_og = readRDS("fogt_rebuild_results.RDS")
res_data_red_og = readRDS("fogt_rebuild_results_lindetrend.RDS")
res_data_red_og$method = "original"
# combine results
res_data_red = rbind(res_data_red, res_data_red_og)
res_data_red$method = factor(res_data_red$method, 
                             levels = c("original", "permuted"))
res_data_red$pretty_sector %>% table()
res_data_red$pretty_sector[res_data_red$pretty_sector == "Bellingshausen Amundsen Sea"] = "Bellingshausen\nAmundsen Sea"
res_data_red$pretty_sector[res_data_red$pretty_sector == "King Hakon VII"] = "King Haakon VII"

library("ggplot2")

# plot combined results
ggplot(res_data_red[res_data_red$sector != "total",]) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~pretty_sector, nrow = 1) +
  geom_boxplot(aes(x = season, y = val_cor, col = method)) +
  scale_color_manual(values = c("firebrick", "steelblue"), 
                     labels = c("original", "permuted")) +
  scale_x_discrete("Season") + 
  scale_y_continuous("Validation Correlation") +
  theme_bw() +
  coord_cartesian(ylim = c(-0.5, 1))
ggsave("plots/fogt_rebuild_valCor_permuted.png", height = 4, width = 12)

ggplot(res_data_red[res_data_red$sector != "total",]) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~pretty_sector, nrow = 1) +
  geom_boxplot(aes(x = season, y = val_ce, col = method)) +
  scale_color_manual(values = c("firebrick", "steelblue"),
                     labels = c("original", "permuted")) +
  scale_x_discrete("Sector") +
  scale_y_continuous("Validation CE") +
  theme_bw() +
  coord_cartesian(ylim = c(-0.5, 1))
ggsave("plots/fogt_rebuild_valCe_permuted.png", height = 4, width = 12)
