library(tidyverse)

# slice
# arrange
# sample_frac(d, 0.001)
# sample_n(d, 10)
#pivot_longer/wider

d <- select(db, Fam, Gen, Spec, mmm, fff, sm, sf, jjj, 
       N, E, Date_ok, Biotop)
select(db, -Et, -`Et prev`, -CT) 
mutate(select(d, -Fam, -Gen, -Spec, -Date_ok, -Biotop), 
       abundance = mmm + fff + sm + sf + jjj)
d[,4:8][is.na(d[,4:8])] <- 0

d %>% 
    select(-Fam, -Gen, -Spec, -Date_ok, -Biotop) %>% 
    mutate(abundance = mmm + fff + sm + sf + jjj)
d %>% 
    select(-Fam, -Gen, -Spec, -Date_ok, -Biotop) %>% 
    transmute(N, E, 
        `abundance 1` = mmm + fff + sm + sf + jjj) %>% 
    rename(abu = `abundance 1`) 
d %>% group_by(Fam) %>% 
    summarise(rec <- n(), .groups = "drop")

d %>% mutate(abu = mmm + fff + sm + sf + jjj) %>% 
    group_by(Fam) %>% 
    summarise(rec = n(), abu2 = sum(abu), .groups = "drop")
# "keep", "drop_last"
d %>% mutate(abu = mmm + fff + sm + sf + jjj) %>% 
    group_by(Fam) %>% 
    summarise(rec = n(), abu2 = sum(abu)) %>% 
    arrange(desc(abu2)) # arrange(abu2)

d %>% transmute(sp = paste0(Gen, "_", Spec), 
    abu = mmm + fff + sm + sf + jjj, 
    N, E, Date_ok, Biotop) %>% 
    group_by(sp, N, E, Biotop, Date_ok) %>% 
    summarise(abu = sum(abu), .groups = "drop") %>% 
    #ungroup() 
    group_by(N, E, Biotop, Date_ok) %>% 
    mutate(alpha = n()) %>% 
    ungroup() %>% 
    filter(alpha > 3, Biotop != "Faunistic record") %>%  
    group_by(Biotop) %>% 
    summarise(gamma = n(), mean.alpha = mean(alpha), 
              betha.whittaker = gamma / mean.alpha) %>% 
    filter(gamma != mean.alpha)

d %>% transmute(Species = paste0(Gen, " ", Spec), 
                abu = mmm + fff + sm + sf + jjj, 
                N, E, Date_ok, Biotop) %>% 
    left_join(traits, by = "Species") %>% 
    filter(!is.na(hunting1)) %>% 
    group_by(Biotop, hunting1) %>% 
    summarise(abu = sum(abu)) %>% 
    mutate(hunt = round(abu/sum(abu)*100))

    

    


