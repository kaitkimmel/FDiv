############# TRAIT NUMBER GRAPHS ##################

#Libraries
library(here)
library(ggplot2)
library(ggpubr)

#Pull data from stats
source(here("R/n_traitStats_simplified.R"))
source(here("R/n_traitStats_euc_simplified.R"))

#### GOWER ####
Frich_g <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.1, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.2, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.3, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.4, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.blue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.black, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = fr.1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FRic, color = community), data = cdr.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = fr.blue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FRic, color = community), data = sev.full, size = 3) + 
  #geom_point(aes(x = n_trait, y = fit, color = community), data = fr.blue, size = 3) + 
  #geom_point(aes(x = n_trait, y = fit, color = community), data = fr.black, size = 3) +
  scale_color_manual(values = c("#9D8F0F","#000000", "navyblue", "#00B7FF","#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "FRich") +
  theme_pubr()

FEve_g <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.1, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.2, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.3, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.4, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.blue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.black, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = fe.1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FEve, color = community), data = cdr.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = fe.blue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FEve, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#000000", "navyblue", "#00B7FF","#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "FEve") +
  theme_pubr()

FDis_g <-ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.1, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.2, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.3, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.4, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.blue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.black, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = fdis.1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FDis, color = community), data = cdr.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.blue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FDis, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#000000", "navyblue", "#00B7FF","#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "FDis") +
  theme_pubr()

FDiv_g <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.1, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.2, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.3, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.4, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.blue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.black, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = fdiv.1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FDiv, color = community), data = cdr.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.blue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FDiv, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#000000", "navyblue", "#00B7FF","#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "FDiv") +
  theme_pubr()

RaoQ_g <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.1, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.2, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.3, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.4, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.blue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.black, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = rq.1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = RaoQ, color = community), data = cdr.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = rq.blue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = RaoQ, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#000000", "navyblue", "#00B7FF","#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "Rao Q") +
  theme_pubr()

kderichness_g <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.1, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.2, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.3, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.4, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.blue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.black, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.alpha, color = community), data = cdr.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.blue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.alpha, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#000000", "navyblue", "#00B7FF","#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Richness") +
  theme_pubr()

kdeevenness_g <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.1, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.2, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.3, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.4, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.blue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.black, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.evenness, color = community), data = cdr.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.blue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.evenness, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#000000", "navyblue", "#00B7FF","#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Evenness") +
  theme_pubr()

kdedispersion_g <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.1, alpha = 0.1, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.2, alpha = 0.1, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.3, alpha = 0.1, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.4, alpha = 0.1, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.blue, alpha = 0.1, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.black, alpha = 0.1, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.dispersion, color = community), data = cdr.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.blue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.black, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.dispersion, color = community), data = sev.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#000000", "navyblue", "#00B7FF","#731279", "#075A13"), name = "Community") +
  labs(x = "Number of Traits", y = "KDE Dispersion") +
  theme_pubr()


#### EUCLIDEAN ####
Frich_e <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.e1, alpha = 0.3, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.e2, alpha = 0.3, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.e3, alpha = 0.3, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.e4, alpha = 0.3, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.eblue, alpha = 0.3, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fr.eblack, alpha = 0.3, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = fr.e1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.e2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.e3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FRic, color = community), data = cdre.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = fr.eblue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fr.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FRic, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#000000", "navyblue", "#00B7FF","#731279", "#075A13"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "FRich") +
  theme_pubr()

FEve_e <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.e1, alpha = 0.3, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.e2, alpha = 0.3, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.e3, alpha = 0.3, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.e4, alpha = 0.3, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.eblue, alpha = 0.3, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fe.eblack, alpha = 0.3, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = fe.e1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.e2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.e3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FEve, color = community), data = cdre.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = fe.eblue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fe.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FEve, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#000000", "navyblue", "#00B7FF","#731279", "#075A13"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "FEve") +
  theme_pubr()

FDis_e <-ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.e1, alpha = 0.3, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.e2, alpha = 0.3, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.e3, alpha = 0.3, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.e4, alpha = 0.3, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.eblue, alpha = 0.3, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdis.eblack, alpha = 0.3, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = fdis.e1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.e2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.e3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FDis, color = community), data = cdre.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.eblue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fdis.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FDis, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#000000", "navyblue", "#00B7FF","#731279", "#075A13"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "FDis") +
  theme_pubr()

FDiv_e <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.e1, alpha = 0.3, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.e2, alpha = 0.3, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.e3, alpha = 0.3, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.e4, alpha = 0.3, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.eblue, alpha = 0.3, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = fdiv.eblack, alpha = 0.3, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = fdiv.e1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.e2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.e3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = FDiv, color = community), data = cdre.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.eblue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = fdiv.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = FDiv, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#000000", "navyblue", "#00B7FF","#731279", "#075A13"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "FDiv") +
  theme_pubr()

RaoQ_e <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.e1, alpha = 0.3, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.e2, alpha = 0.3, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.e3, alpha = 0.3, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.e4, alpha = 0.3, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.eblue, alpha = 0.3, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = rq.eblack, alpha = 0.3, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = rq.e1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.e2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.e3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = RaoQ, color = community), data = cdre.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = rq.eblue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = rq.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = RaoQ, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#000000", "navyblue", "#00B7FF","#731279", "#075A13"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "Rao Q") +
  theme_pubr()

kderichness_e <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.e1, alpha = 0.3, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.e2, alpha = 0.3, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.e3, alpha = 0.3, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.e4, alpha = 0.3, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.eblue, alpha = 0.3, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.alpha.eblack, alpha = 0.3, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.e1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.e2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.e3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.alpha, color = community), data = cdre.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.eblue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.alpha.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.alpha, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#000000", "navyblue", "#00B7FF","#731279", "#075A13"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "KDE Richness") +
  theme_pubr()

kdeevenness_e <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.e1, alpha = 0.3, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.e2, alpha = 0.3, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.e3, alpha = 0.3, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.e4, alpha = 0.3, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.eblue, alpha = 0.3, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.evenness.eblack, alpha = 0.3, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.e1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.e2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.e3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.evenness, color = community), data = cdre.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.eblue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.evenness.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.evenness, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#000000", "navyblue", "#00B7FF","#731279", "#075A13"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "KDE Evenness") +
  theme_pubr()

kdedispersion_e <- ggplot() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.e1, alpha = 0.3, fill = "#9D8F0F") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.e2, alpha = 0.3, fill = "#731279") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.e3, alpha = 0.3, fill = "#00B7FF") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.e4, alpha = 0.3, fill = "#075A13") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.eblue, alpha = 0.3, fill = "navyblue") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, x = n_trait), data = kde.dispersion.eblack, alpha = 0.3, fill = "black") +
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.e1, lwd = 2, color = "#9D8F0F") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.e2, lwd = 2, color = "#731279") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.e3, lwd = 2, color = "#00B7FF") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.e4, lwd = 2, color = "#075A13") + 
  geom_point(aes(x = n_trait, y = kde.dispersion, color = community), data = cdre.full, size = 3) + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.eblue, lwd = 2, color = "navyblue") + 
  geom_line(aes(x= n_trait, y = fit), data = kde.dispersion.eblack, lwd = 2, color = "black") + 
  geom_point(aes(x = n_trait, y = kde.dispersion, color = community), data = seve.full, size = 3) + 
  scale_color_manual(values = c("#9D8F0F","#000000", "navyblue", "#00B7FF","#731279", "#075A13"), name = "Community") +
  labs(x = "Maximum Trait-Trait Correlation", y = "KDE Dispersion") +
  theme_pubr()

png(here('Figures/n_trait_full.png'), height = 9, width = 12, units = 'in', res = 150)
ggarrange(plotlist = list(Frich_g,kderichness_g, Frich_e,kderichness_e, FEve_g,  
                          kdeevenness_g, FEve_e,kdeevenness_e, FDis_g,  FDiv_g, FDis_e,FDiv_e, 
                          RaoQ_g, kdedispersion_g,RaoQ_e, kdedispersion_e), ncol = 4, nrow = 4, 
          labels = c("A", "B", 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P'), common.legend = TRUE)
dev.off()
