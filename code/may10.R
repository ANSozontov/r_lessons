library(tidyverse)
theme_set(theme_bw())
# data load ---------------------------------------------------------------
df <- readxl::read_xlsx("testdata.xlsx", sheet = "main_data")
tr <- gsheet::gsheet2tbl('https://docs.google.com/spreadsheets/d/1mhCBWmbA2gOCmiXJYfY6LHzE-kumEcK2cbeuA_kYVI4')
pl <- rio::import('https://github.com/ANSozontov/Mord_State_Reserve/raw/master/testdata.xlsx', 
                  sheet = 3, skip = 1)

# rotation ----------------------------------------------------------------
t(df[,7:64]) %>% View

# long & wide -------------------------------------------------------------
pivot_longer(df, names_to = "taxa", values_to = "num", -c("year", "season", "site", "plot", "trap", "duration"))
pivot_longer(df, names_to = "taxa", values_to = "num", -c(1:6))

# get.trap-days -----------------------------------------------------------
long <- 
  pivot_longer(df, names_to = "taxa", values_to = "num", -c(1:6)) %>% 
  group_by(year, season, site, plot, taxa, duration) %>% 
  summarise(num = sum(num), n = n(), .groups = "drop") %>% 
  mutate(lovsut = n*duration) %>% 
  mutate(num = num/lovsut*100) %>% 
  select(-n, -lovsut, -duration)

# class.div ---------------------------------------------------------------
div <- long %>% 
  group_by(year, season, site, plot) %>% 
  summarise(abu = sum(num), nsp = length(num[num>0]), .groups = "drop") %>% 
  mutate(mar = nsp/sqrt(abu), men = nsp/log(abu))

# join spec_traits --------------------------------------------------------
left_join(long, tr, by = "taxa") %>% 
  select(year, season, site, plot, num, wings.brach_02, wings.macro_02) %>% 
  # hum.xero_03, hum.meso_03, hum.hygro_03
  pivot_longer(names_to = "wings", values_to = "val", -c(1:5)) %>% 
  filter(num > 0, val > 0) %>% 
  group_by(season, wings) %>% 
  summarise(mean.w = mean(num), ci = sd(num), 
            .groups = "drop_last") #%>%  mutate(mean.w = mean.w/sum(mean.w)*100)

# join site chars ---------------------------------------------------------
div <- left_join(div, pl, by = c("site", "plot"))

# rows subset -------------------------------------------------------------
long %>% slice(1:10)
long %>% sample_n(10)
long %>% sample_frac(0.001)

# repeat ------------------------------------------------------------------
# betadiv (the same)
bet <- numeric()
while (length(bet)<1000) {
b <- df[,7:ncol(df)] %>% 
  sample_n(nrow(.)-1) %>% 
  select(sample(1:ncol(.), ncol(.) - 1))
b[b>0] <- 1
a <- mean(apply(b, 1, sum))
g <- apply(b, 2, sum)
g <- length(g[g>0])
bet <- c(bet, g/a)
}
quantile(bet, c(0.025, 0.975))

library(foreach)
bet <- numeric()
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
bet <- foreach(i = 1:5000, .combine = 'c') %dopar% {
  b <- dplyr::select(dplyr::sample_n(df[,7:ncol(df)], nrow(df)-1), 
                     sample(1:58, 57))
  b[b>0] <- 1
  a <- mean(apply(b, 1, sum))
  g <- apply(b, 2, sum)
  g <- length(g[g>0])
  g/a
}
parallel::stopCluster(cl)
quantile(bet, c(0.025, 0.975))

#lm (modify)
res <- expand.grid(otk = c("abu", "nsp", "mar", "men"), 
  pred = c("km", "nspec", "proec.area", "biomas", "Cu", "Pb", "Cd")) %>% 
  mutate(pred = as.character(pred), otk = as.character(otk), 
         pval = NA, r2 = NA, aic = NA, test = NA)
for(i in 1:nrow(res)) { 
  fit <- select(div, res$pred[i], res$otk[i])
  colnames(fit) <- c("pred", "otk")
  fit <- lm(otk ~ pred, data = fit)
  res$aic[i] <- round(AIC(fit))
  res$test[i] <- round(shapiro.test(fit$residuals)$p.value, 4)
  fit <- summary(fit)
  res$r2[i] <- round(fit$r.squared, 2)
  fit <- fit$fstatistic
  res$pval[i] <- round(pf(fit[1], fit[2], fit[3], lower.tail=F), 5) + 0.0001
}

# B.3. Advanced -----------------------------------------------------------
#p.adj
round(p.adjust(res$pval, method = "BY"), 4) # c("holm", "hommel", "BH", "fdr")
round(p.adjust(res$pval, method = "hochberg"), 4)
round(p.adjust(res$pval, method = "bonferroni"), 4)

res %>% 
  mutate(p.adj = round(p.adjust(res$pval, method = "BY"), 4))

car::qqPlot(div$nsp)

# hill
vegan::renyi(df[,7:64], scales = c(0,1, 2, 3, 10, Inf), hill = TRUE) # %>% as_tibble %>% filter(`0` > 0)
df %>% mutate(nsp = apply(df[,7:64], 1, function(a){length(a[a>0])})) %>% 
  filter(nsp > 0) %>% 
  select(7:64) %>% 
  vegan::renyi(., scales = c(0, 1, 2, 3, Inf), hill = TRUE)
vegan::renyi(df[,7:64], scales = c(0, 1, 2, 3, Inf), hill = TRUE) %>% 
  cbind(df[,1:6], .) %>% 
  as_tibble %>% 
  filter(`0` > 0, site %in% c("K27S", "K32N", "K04S", "K05N")) %>% 
  select(3, 7:11) %>% 
  pivot_longer(names_to = "ord", values_to = "hil", -site) %>% 
  group_by(site, ord) %>% 
  summarise( ci = sd(hil), hil = mean(hil), .groups = "drop") %>% 
  transmute(site, ord, h = paste0(round(hil, 1), "±", round(ci, 1))) %>% 
  pivot_wider(names_from = ord, values_from = h)

#FD
dis <- tr %>% 
  arrange(taxa) %>% 
  mutate(size.own_01 = log(tr$size.own_01)) %>% 
  column_to_rownames("taxa") %>% 
  mutate_all(function(a){as.vector(scale(a))}) %>% 
  vegan::vegdist(method = "gower")
df3 <- df %>% mutate(nsp = apply(df[,7:64], 1, function(a){length(a[a>0])})) %>% 
  filter(nsp > 0) %>% 
  mutate(id = paste0(year, site, 1:nrow(.))) %>% 
  select(7:64, 66) %>% 
  pivot_longer(names_to = "sp", values_to = "v", -id) %>% 
  arrange(sp) %>% 
  pivot_wider(names_from = "sp", values_from = "v") %>% 
  column_to_rownames("id") # %>% as.matrix()
fd <- FD::dbFD(dis, df3, calc.CWM = FALSE)
fd <- map_df(fd, cbind) %>% 
  cbind(id = names(fd$nbsp)) %>% 
  mutate(year = as.numeric(substr(id, 1, 4)), site = substr(id, 5, 8)) 
# ggplot(aes(x = site, color = site, y = FEve)) + geom_boxplot()

#raref
library(iNEXT)
df2 <- long %>% 
  transmute(id = paste0(year, site, season, plot), taxa, num) %>% 
  pivot_wider(names_from = taxa, values_from = num, values_fill = 0, values_fn = sum) %>% 
  column_to_rownames("id") %>% 
  mutate(nsp = apply(., 1, function(a){length(a[a>0])})) %>% 
  filter(nsp > 1) %>% 
  select(-nsp)
rar <- df2 %>% 
  t %>% 
  as.data.frame %>% 
  lapply(function(a){a[a>0]}) %>% 
  keep(function(x){length(x) > 1}) %>% 
  iNEXT(q = 0, datatype = "abundance", size = c(10, 20, 50, 100, 500, 2000))
#q == 0, 1, 2, reduce size vector in in case of q > 0

rar <- rar$iNextEst %>% 
  reduce(rbind) %>% 
  mutate(id = rep(names(rar$iNextEst), times= sapply(rar$iNextEst, nrow))) %>% 
  mutate(m = case_when(method == "observed" ~ 1, TRUE ~ m)) %>% 
  as_tibble() %>% 
  filter(m %in% c(10, 20, 50, 100, 500, 2000) | method == "observed") %>% 
  transmute(m = paste0("m", m/1000), order, qD, year = substr(id, 1, 4), 
            site = substr(id, 5, 8)) 

rar %>% # filter(m == "m1") %>% ##
  filter(qD > 0) %>% 
  ggplot(aes(x = site, y = qD, color = site)) + 
  geom_boxplot() + 
  labs(x = "", y = "n_spec") + facet_wrap(~m)

#permanova
fctrs <- df2 %>% 
  rownames_to_column("id") %>% 
  as_tibble %>% 
  transmute(year = substr(id, 1, 4), 
            site = substr(id, 5, 8), 
            season = substr(id, 9, nchar(id)-1), 
            plot = substr(id, nchar(id), nchar(id)))
dis <- vegan::vegdist(df2, method = "bray", binary = FALSE)

vegan::adonis(dis ~ year + site + season, data = fctrs, permutations = 499) # site:season


# save_all ----------------------------------------------------------------
rm(b, cl, a, bet, dis, fit, g, i)
save.image(file = "O://!_R10/export.RData")















