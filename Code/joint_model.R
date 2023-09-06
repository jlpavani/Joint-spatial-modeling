library(dplyr)
library(rgdal)
library(spdep)
library(leaflet)
library(leafsync)
library(INLA)
library(INLAMSM)

# **************************************************************************** #
# Load Data
# **************************************************************************** #
load("dta_deng_chik.RData")
data <- data[order(data$id),]

Output.Areas <- readOGR(".", "CE_Municipios_2019") 

# Adjancency Matrix
W.nb <- poly2nb(Output.Areas, row.names = 1:length(Output.Areas)) 
W.sp <- as(nb2mat(W.nb, style = "B"), "Matrix")

# **************************************************************************** #
# Joint model
# **************************************************************************** #
# Intercept
data$intercept <- as.factor(data$disease)

# Dummy indices for space
data$s.dummy <- NA

# Spacial indices for specific effects
data$s.1 <- NA
data$s.1[data$disease == "Dengue"] <- as.numeric(as.factor(data$id_area[data$disease == "Dengue"]))

data$s.2 <- NA
data$s.2[data$disease == "Chikungunya"] <- as.numeric(as.factor(data$id_area[data$disease == "Chikungunya"]))

# Indices for spatial disease-specific effects
data$id_area1 <- data$s.1; data$id_area2 <- data$s.2

# Log-Normal prior with zero mean and precision 1/5.9 on spatial and temporal weights
prior.beta.s <- list(prior = "normal", param = c(0, 1 / 5.9), fixed = FALSE, initial = 0.01)

# Flat prior on sigma:
prior.prec <- list(prior = "expression: logdens = -log_precision / 2; return(logdens)", initial = 0)

inla.scale <- FALSE

formula <- observed ~ -1 + temperature + rurality + intercept + 
           f(id_area1, model = "besag",  scale.model = inla.scale, graph = W.sp, hyper = list(prec = prior.prec)) + 
           f(id_area2, model = "besag",  scale.model = inla.scale, graph = W.sp, hyper = list(prec = prior.prec)) + 
           f(s.dummy, model = "besag",  scale.model = inla.scale, graph = W.sp, hyper = list(prec = prior.prec)) + 
           f(s.1, copy = "s.dummy",  range = c(0, Inf), hyper = list(beta = prior.beta.s)) + 
           f(s.2, copy = "s.dummy",  range = c(0, Inf), hyper = list(beta = prior.beta.s))

out_fit <- inla(formula, data = data, E = expected, family = "poisson", 
                verbose = FALSE, control.predictor = list(compute=TRUE), 
                control.compute=list(dic = TRUE, waic = TRUE, config = TRUE),
                control.inla(strategy = "laplace"))

out_fit_poi_final <- inla.rerun(out_fit)
summary(out_fit) 

# **************************************************************************** #
# Results: Spatial Effects
# **************************************************************************** #
mycolors <- c("#063c63", "#0f6aaa", "#82b9e0", "#a2d0f2", "#ed5353", "#e21616",
              "#a21010", "#610909")

# Shared
Output.Areas@data$shared <- out_fit$summary.random$s.dummy[, "mean"]
colours <- colorBin(palette = mycolors, domain = Output.Areas@data$shared, 
                    bins = seq(min(Output.Areas@data$shared), max(Output.Areas@data$shared), len=9))
map_shared <- leaflet(data=Output.Areas, options = leafletOptions(zoomControl = FALSE)) %>% 
              addTiles() %>%
              addPolygons(fillColor = ~colours(shared), color="white", weight=1, fillOpacity = 1) %>% 
              addLegend(pal = colours, labFormat = labelFormat(digits = 1), values = Output.Areas@data$shared,
                        opacity = 1, title="Shared", position = "bottomright") %>%
              addScaleBar(position="bottomleft")

# Specific - dengue
Output.Areas@data$dengue <- out_fit$summary.random$id_area1[, "mean"]
colours <- colorBin(palette = mycolors, domain = Output.Areas@data$dengue, 
                    bins = seq(min(Output.Areas@data$dengue), max(Output.Areas@data$dengue), len=9))
map_dengue_sp <- leaflet(data=Output.Areas, options = leafletOptions(zoomControl = FALSE)) %>% 
                 addTiles() %>%
                 addPolygons(fillColor = ~colours(dengue), color="white", weight=1, fillOpacity = 1) %>% 
                 addLegend(pal = colours, labFormat = labelFormat(digits = 1), values = Output.Areas@data$dengue,
                           opacity = 1, title="Dengue", position = "bottomright") %>%
                 addScaleBar(position="bottomleft")

# Specific - chikungunya
Output.Areas@data$chiku <- out_fit$summary.random$id_area2[, "mean"]
colours <- colorBin(palette = mycolors, domain = Output.Areas@data$chiku, 
                    bins = seq(min(Output.Areas@data$chiku), max(Output.Areas@data$chiku), len=9))
map_chiku_sp <- leaflet(data=Output.Areas, options = leafletOptions(zoomControl = FALSE)) %>% 
                addTiles() %>%
                addPolygons(fillColor = ~colours(chiku), color="white", weight=1, fillOpacity = 1) %>% 
                addLegend(pal = colours, labFormat = labelFormat(digits = 1), values = Output.Areas@data$chiku,
                          opacity = 1, title="Chikungunya", position = "bottomright") %>%
                addScaleBar(position="bottomleft")

# Total - dengue
Output.Areas@data$dengue_tot <- out_fit$summary.random$id_area1[, "mean"] + out_fit$summary.random$s.1[, "mean"]
colours <- colorBin(palette = mycolors, domain = Output.Areas@data$dengue_tot, 
                    bins = seq(min(Output.Areas@data$dengue_tot), max(Output.Areas@data$dengue_tot), len=9))
map_dengue_tot <- leaflet(data=Output.Areas, options = leafletOptions(zoomControl = FALSE)) %>% 
                  addTiles() %>%
                  addPolygons(fillColor = ~colours(dengue_tot), color="white", weight=1, fillOpacity = 1) %>% 
                  addLegend(pal = colours, labFormat = labelFormat(digits = 1), values = Output.Areas@data$dengue_tot,
                            opacity = 1, title="Dengue", position = "bottomright") %>%
                  addScaleBar(position="bottomleft")

# Total - chikungunya
Output.Areas@data$chiku_tot <- out_fit$summary.random$id_area2[, "mean"] + out_fit$summary.random$s.2[, "mean"]
colours <- colorBin(palette = mycolors, domain = Output.Areas@data$chiku_tot, 
                    bins = seq(min(Output.Areas@data$chiku_tot), max(Output.Areas@data$chiku_tot), len=9))
(map_chiku_tot <- leaflet(data=Output.Areas, options = leafletOptions(zoomControl = FALSE)) %>% addTiles() %>%
                  addPolygons(fillColor = ~colours(chiku_tot), color="white", weight=1, fillOpacity = 1) %>% 
                  addLegend(pal = colours, labFormat = labelFormat(digits = 1), values = Output.Areas@data$chiku_tot, 
                            opacity = 1, title="Chikungunya", position = "bottomright") %>%
                  addScaleBar(position="bottomleft")%>% 
                  setView(lng = -39, lat = -5.5, zoom = 7 ) %>%
                  setMaxBounds( lng1 = -39.4, lat1 = -5.5, lng2 = -39.2, lat2 = -5))

# **************************************************************************** #
# Results: Relative Risks
# **************************************************************************** #
Output.Areas@data$RR_d <- out_fit$summary.fitted.values[1:184,"mean"]
Output.Areas@data$LL_d <- out_fit$summary.fitted.values[1:184,"0.025quant"]
Output.Areas@data$UL_d <- out_fit$summary.fitted.values[1:184,"0.975quant"]

colours <- colorBin(palette = mycolors, domain = Output.Areas@data$RR_d, 
                    bins = c(0, 0.3, 0.6, 0.975, 1.025, 1.5, 2, max(Output.Areas@data$RR_d)))
(map_RR_d2 <- leaflet(data=Output.Areas) %>% addTiles() %>%
              addPolygons(fillColor = ~colours(RR_d), color="white", weight=1, fillOpacity = 1) %>%
              addLegend(pal = colours, values = Output.Areas@data$RR_d,opacity = 1,
                        labFormat = labelFormat(digits = 1), title="Dengue") %>%
              addScaleBar(position="bottomleft") %>% 
              setView(lng = -39, lat = -5.5, zoom = 7 ) %>%
              setMaxBounds( lng1 = -39.4, lat1 = -5.5, lng2 = -39.2, lat2 = -5))


Output.Areas@data$RR_c <- out_fit$summary.fitted.values[185:368,"mean"]
Output.Areas@data$LL_c <- out_fit$summary.fitted.values[185:368,"0.025quant"]
Output.Areas@data$UL_c <- out_fit$summary.fitted.values[185:368,"0.975quant"]

colours <- colorBin(palette = mycolors, domain = Output.Areas@data$RR_c, 
                    bins = c(0, 0.3, 0.6, 0.975, 1.025, 1.5, 2, max(Output.Areas@data$RR_c)))
(map_RR_c2 <- leaflet(data=Output.Areas) %>% addTiles() %>%
              addPolygons(fillColor = ~colours(RR_c), color="white", weight=1, fillOpacity = 1) %>%
              addLegend(pal = colours, values = Output.Areas@data$RR_c, opacity = 1,
                        labFormat = labelFormat(digits = 1), title="Chikungunya") %>%
              addScaleBar(position="bottomleft")%>% 
              setView(lng = -39, lat = -5.5, zoom = 7 ) %>%
              setMaxBounds( lng1 = -39.4, lat1 = -5.5, lng2 = -39.2, lat2 = -5))

# **************************************************************************** #
# Results: Exceedance Probabilities
# **************************************************************************** #
Output.Areas@data$exc_prob_d <- sapply(out_fit$marginals.fitted.values, 
                                       FUN = function(marg){1 - inla.pmarginal(q = 1, 
                                                                               marginal = marg)})[1:184]

colours <- colorBin(palette = mycolors, domain = Output.Areas@data$exc_prob_d, 
                    bins = seq(min(Output.Areas@data$exc_prob_d),
                               max(Output.Areas@data$exc_prob_d), len=6))
(map_exc_prob_d <- leaflet(data=Output.Areas) %>% addTiles() %>%
                   addPolygons(fillColor = ~colours(exc_prob_d), color="white", weight=1, fillOpacity = 1) %>% 
                   addLegend(pal = colours, labFormat = labelFormat(digits = 1), 
                             values = Output.Areas@data$exc_prob_d, opacity = 1, title="Dengue") %>%
                   addScaleBar(position="bottomleft") %>% 
                   setView(lng = -39, lat = -5.5, zoom = 7 ) %>%
                   setMaxBounds( lng1 = -39.4, lat1 = -5.5, lng2 = -39.2, lat2 = -5))

Output.Areas@data$exc_prob_c <- sapply(out_fit$marginals.fitted.values, 
                                       FUN = function(marg){1 - inla.pmarginal(q = 1, 
                                                                               marginal = marg)})[185:368]

colours <- colorBin(palette = mycolors, domain = Output.Areas@data$exc_prob_c, 
                    bins = seq(min(Output.Areas@data$exc_prob_c),
                               max(Output.Areas@data$exc_prob_c), len=6))
(map_exc_prob_c <- leaflet(data=Output.Areas) %>% addTiles() %>%
                   addPolygons(fillColor = ~colours(exc_prob_c), color="white", weight=1, fillOpacity = 1) %>% 
                   addLegend(pal = colours, labFormat = labelFormat(digits = 1), values = Output.Areas@data$exc_prob_c,
                             opacity = 1, title="Chikungunya") %>%
                   addScaleBar(position="bottomleft") %>%
                   setView(lng = -39, lat = -5.5, zoom = 7 ) %>%
                   setMaxBounds( lng1 = -39.4, lat1 = -5.5, lng2 = -39.2, lat2 = -5))