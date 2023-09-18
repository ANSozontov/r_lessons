# load --------------------------------------------------------------------
library(readxl)
library(tidyverse)
# wos1 <- read_excel("D:/wos_2.xlsx", sheet = 1)
# wos2 <- read_excel("D:/wos_2.xlsx", sheet = 2)
# within(wos1, Category <-data.frame(do.call('rbind', 
#   strsplit(as.character(Category), '; ', fixed=TRUE))))

db <- read_excel("O:/P_Monogrpah/DB 10.3 DarwinCore.xlsx", sheet = "data")
sp <- read_excel("O:/P_Monogrpah/DB 10.3 DarwinCore.xlsx", sheet = "species")

# columns -----------------------------------------------------------------


select(db, Fam, Gen, Spec, mmm, fff, sm,
       sf, jjj, Locality, Date_ok, Biotop) 
select(mutate(db, sp = paste0(Gen, " ", Spec)), 
       sp, mmm, fff, sm, sf, juv = jjj)
select(db, 9:12, 5:8)



# sum(c(1, 2,6))
# c(1, 2,6) %>% sum(.)

db %>% 
    mutate(sp = paste0(Gen, " ", Spec)) 
    select(sp, mmm, fff, sm, sf, jjj)

db %>% 
    transmute(
              mmm, fff, sm, sf, jjj, 
              sp = paste0(Gen, " ", Spec)) %>% 
    rename(speCC123 = sp)

# rows --------------------------------------------------------------------
d <- db %>% 
    select(Fam, Gen, Spec, mmm, fff, sm,
       sf, jjj, Locality, Date_ok, Biotop) 
d %>% 
    # sample_n(5)
    sample_frac(0.001)

d %>% slice(1001, 10, 1)

d$Fam == "Pholcidae"
which(d$Fam == "Pholcidae")

d %>% slice(which(d$Fam == "Pholcidae"))

d %>% filter(Fam == "Clubionidae", mmm >= 1, Spec == "sp.")

d %>% filter(Fam == "Clubionidae", mmm >= 1, Spec == "sp.") %>% 
    View()

d %>% filter(Fam == "Clubionidae", mmm >= 1, Spec == "sp.") %>% 
    pull(Date_ok) %>% 
    unique() %>% 
    sort()

# summarizing -------------------------------------------------------------

d$mmm[is.na(d$mmm)] <- 0
d$fff[is.na(d$fff)] <- 0

d %>% group_by(Fam) %>% 
    summarise(n_m = sum(mmm), n_f = sum(fff)) %>% 
    ungroup()

d %>% group_by(Fam) %>% 
    summarise(n_m = sum(mmm), n_f = sum(fff), .groups = "drop")

d %>% group_by(Fam, Biotop) %>% 
    summarise(n_m = sum(mmm), n_f = sum(fff),
              sdm = sd(mmm), sdf = sd(fff), .groups = "drop") %>%
    filter(n_m >= 1 | n_f >1) %>%
    transmute(Fam, Biotop, mm = paste0(n_m, " ± ", round(sdm, 2)), 
              ff = paste0(n_f, " ± ", round(sdf, 2))
    ) %>% 
    arrange(desc(Biotop))
    




