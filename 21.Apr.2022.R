library(tidyverse)
theme_set(theme_bw())
nsp <- function(a){length(a[a>0])}

# data load ---------------------------------------------------------------
df <- readxl::read_xlsx("data_example_22.03.2022.xlsx", sheet = "main")
tr <- rio::import('https://docs.google.com/spreadsheets/d/1nm2oT4Ps4gyN2Bzn3U-_1-0lkXKOYYny1JKtHn5MFrk') |> as_tibble()
pl <- rio::import('https://github.com/ANSozontov/R_lessons/raw/master/data_example_22.03.2022.xlsx', 
                  sheet = 2) |> as_tibble()

apply(df[,-1], 1, sum)
tibble(df[,1], N = apply(df[,-1], 1, sum)) %>% 
    left_join(pl, by = "id") %>% 
    filter(method != "manual collecting") %>% 
    mutate(biotop = paste0(biotop, ": ", year)) %>% 
    ggplot(aes(x = biotop, y = N)) + 
    geom_boxplot()

# long & wide + future ---------------------------------------------
df %>% 
    pivot_longer(names_to = "taxa", values_to = "num", -id) 

df %>% 
    pivot_longer(names_to = "taxa", values_to = "num", -id) %>% 
    left_join(pl, by = "id") %>% 
    left_join(tr, by = "taxa") %>% 
    group_by(biotop, storey) %>% 
    summarise(N = sum(num), .groups = "drop")

df %>% 
    pivot_longer(names_to = "taxa", values_to = "num", -id) %>% 
    left_join(pl, by = "id") %>% 
    left_join(tr, by = "taxa") %>% 
    filter(habitats != "water") %>% 
    mutate(biotop = paste0(biotop, ": ", year)) %>% 
    group_by(biotop, storey) %>% 
    summarise(N = sum(num), .groups = "drop") %>% 
    ggplot(aes(x = biotop, y = N, fill = storey)) + 
        geom_col(position = "fill")

# class.div ---------------------------------------------------------------
df %>% pivot_longer(names_to = "taxa", values_to = "num", -id) %>% 
  left_join(select(pl, id, biotop, year), by = "id") %>% 
  group_by(id, year, biotop) %>% 
  summarise(N = sum(num), S = nsp(num), .groups = "drop") %>% 
  mutate(Margalef = S/sqrt(N), Menhinnick = S/log(N), biotop = paste0(biotop, ": ", year)) %>% 
  pivot_longer(values_to = "val", names_to = "index", -1:-5) %>% 
  ggplot(aes(x = biotop, y = val, color = index)) + 
    geom_boxplot(color = "black", size = 1) + 
    geom_jitter(size = 4, alpha = 0.65)+
    facet_wrap(~index, scales = "free") + 
    theme(legend.position = "none") + 
    labs(x = NULL, y = NULL)

# repeat ------------------------------------------------------------------
# betadiv (the same)

B <- full_join(pl, df, by = "id") %>% 
  select(-id, -dmy, -line, -method, -avg_temp) %>% 
  rename(A = n.traps)
  
b <- B %>% mutate(A = apply(.[,4:ncol(.)], 1, nsp)) %>% 
  group_by(biotop, year) %>% 
  summarise_all(mean)  %>% 
  ungroup %>% 
  unite("biotop", biotop, year) %>% 
  transmute(biotop, A, G = apply(.[,3:ncol(.)], 1, nsp), W = G/A)
b

bet <-tibble(biotop = character(), W = numeric())
set.seed(1)
while (nrow(bet)<50*4) {
  b <- B[sample(1:nrow(B), nrow(B)-1),c(1:3, sample(4:ncol(B), ncol(B)-4) )] %>% 
    mutate(A = apply(.[,4:ncol(.)], 1, nsp)) %>% 
    group_by(biotop, year) %>% 
    summarise_all(mean)  %>% 
    ungroup %>% 
    unite("biotop", biotop, year) %>% 
    transmute(biotop, A, G = apply(.[,3:ncol(.)], 1, nsp), W = G/A)
  bet <- rbind(bet, b[,c(1,4)])
}

bet %>% 
  group_by(biotop) %>% 
  summarise(L = quantile(W, 0.025), R = quantile(W, 0.975)) %>% 
  left_join(b, by = "biotop") %>% 
  ggplot(aes(x = biotop, y = W, ymin = L, ymax = R, fill = biotop)) + 
  geom_col() + 
  geom_pointrange(size = 1) +
  labs(x = NULL, y = "Whittaker's W", title = "Beta-diversity, 250 iterations")

library(foreach)
cl <- parallel::makeCluster(parallel::detectCores()-1)
doParallel::registerDoParallel(cl)
bet.par <- foreach(i = 1:999, .combine = 'rbind') %dopar% {
  library(dplyr)
  B[sample(1:nrow(B), nrow(B)-1),c(1:3, sample(4:ncol(B), ncol(B)-4) )] %>% 
    mutate(A = apply(.[,4:ncol(.)], 1, nsp)) %>% 
    group_by(biotop, year) %>% 
    summarise_all(mean)  %>% 
    ungroup %>% 
    transmute(biotop = paste0(biotop, "_", year), A, G = apply(.[,3:ncol(.)], 1, nsp), W = G/A)
}
parallel::stopCluster(cl)
bet.par %>% 
  group_by(biotop) %>% 
  summarise(L = quantile(W, 0.025), R = quantile(W, 0.975)) %>% 
  left_join(b, by = "biotop") %>% 
  ggplot(aes(x = biotop, y = W, ymin = L, ymax = R, fill = biotop)) + 
  geom_col() + 
  geom_pointrange(size = 1) +
  labs(x = NULL, y = "Whittaker's W", title = "Beta-diversity, 250 iterations")

ggplot(bet.par, aes(x = W, fill = biotop, color = biotop)) + 
  geom_density(kernel = "rectangular") + 
  geom_vline(mapping = aes(xintercept = L, color = biotop), size = 1.1,
             data = summarise(group_by(bet.par, biotop), L = quantile(W, 0.025), R = quantile(W, 0.975)))+ 
  geom_vline(mapping = aes(xintercept = R, color = biotop), size = 1.2,
             data = summarise(group_by(bet.par, biotop), L = quantile(W, 0.025), R = quantile(W, 0.975)))


# DONT RUN ----------------------------------------------------------------
#lm (modify)
# res <- expand.grid(otk = c("abu", "nsp", "mar", "men"), 
#                    pred = c("km", "nspec", "proec.area", "biomas", "Cu", "Pb", "Cd")) %>% 
#     mutate(pred = as.character(pred), otk = as.character(otk), 
#            pval = NA, r2 = NA, aic = NA, test = NA)
# for(i in 1:nrow(res)) { 
#     fit <- select(div, res$pred[i], res$otk[i])
#     colnames(fit) <- c("pred", "otk")
#     fit <- lm(otk ~ pred, data = fit)
#     res$aic[i] <- round(AIC(fit))
#     res$test[i] <- round(shapiro.test(fit$residuals)$p.value, 4)
#     fit <- summary(fit)
#     res$r2[i] <- round(fit$r.squared, 2)
#     fit <- fit$fstatistic
#     res$pval[i] <- round(pf(fit[1], fit[2], fit[3], lower.tail=F), 5) + 0.0001
# }
# 
# #p.adj
# round(p.adjust(res$pval, method = "BY"), 4) # c("holm", "hommel", "BH", "fdr")
# round(p.adjust(res$pval, method = "hochberg"), 4)
# round(p.adjust(res$pval, method = "bonferroni"), 4)
# res %>% 
#     mutate(p.adj = round(p.adjust(res$pval, method = "BY"), 4))
# car::qqPlot(div$nsp)

# B.3. Advanced -----------------------------------------------------------

# hill
df[,-1] %>% 
  vegan::renyi(., scales = c(0, 1, 2, 3, Inf), hill = TRUE) %>% 
  as_tibble

vegan::renyi(df[,-1], scales = c(seq(0, 3, by = 0.1), 5, 10,
        # seq(3.5, 10, 0.5), 
        Inf), hill = TRUE) %>% 
  cbind(df[,1], .) %>% 
  as_tibble %>% 
  pivot_longer(names_to = "X", values_to = "H", -id) %>% 
  mutate(X = as.numeric(X)) %>% 
  left_join(select(pl, id, biotop, year), by = "id") %>% 
  group_by(year, biotop, X) %>% 
  summarise(SD = sd(H), H = mean(H), .groups = "drop") %>% 
  mutate(
    X = case_when(X <= 3 ~ X,
                  X == Inf ~ 6,
                  X == 5 ~ 4, 
                  X == 10 ~ 5),
         biotop = paste0(biotop, "_", year)) %>%
  # filter(biotop == "post-fire forest_2001") %>% View
  filter(X != Inf) %>% 
  ggplot(aes(x = X, y = H, color = biotop, 
             ymin = X - SD, ymax = X + SD)) + 
  geom_ribbon() +
  geom_line() + 
  scale_x_continuous(breaks = c(0:6), labels = c(0:3, 5, 10, Inf))



  filter(`0` > 0, site %in% c("K27S", "K32N", "K04S", "K05N")) %>% 
  select(3, 7:11) %>% 
  pivot_longer(names_to = "ord", values_to = "hil", -site) %>% 
  group_by(site, ord) %>% 
  summarise( ci = sd(hil), hil = mean(hil), .groups = "drop") %>% 
  transmute(site, ord, h = paste0(round(hil, 1), "?", round(ci, 1))) %>% 
  pivot_wider(names_from = ord, values_from = h)

#FD
dis <- tr %>% 
    arrange(taxa) %>% 
    mutate(size_01 = log(tr$size_01)) %>% 
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
ggplot(fd, aes(x = site, color = site, y = FEve)) + geom_boxplot()

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















