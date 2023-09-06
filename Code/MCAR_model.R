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
# MCAR
# **************************************************************************** #
data$s <- NA
data$s[data$disease == "Dengue"] <- as.numeric(as.factor(data$id_area[data$disease == "Dengue"]))
data$s[data$disease == "Chikungunya"] <- as.numeric(as.factor(data$id_area[data$disease == "Chikungunya"]))

# Number of diseases
k <- length(unique(data$disease))

alpha.min <- 0.99; alpha.max <- 1 
model <- inla.MCAR.model(k = k, W = W.sp, alpha.min = alpha.min, 
                         alpha.max = alpha.max)
formula <- observed ~ 0 + disease + f(s, model = model) + temperature + rurality
IMCAR <- inla(formula, data = data, E = expected, family = "poisson",
              control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
              control.predictor = list(compute = TRUE),
              control.inla(strategy = "laplace"))
IMCAR <- inla.rerun(IMCAR)
summary(IMCAR) 

hyper.imcar <- inla.MCAR.transform(IMCAR, k, model = "IMCAR",
                                   alpha.min = alpha.min, alpha.max = alpha.max)

# **************************************************************************** #
# Results: Spatial Effects
# **************************************************************************** #
n <- nrow(W.sp)
mycolors <- c("#063c63", "#0f6aaa", "#82b9e0", "#a2d0f2", "#ed5353", "#e21616",
              "#a21010", "#610909")

Output.Areas@data$IMCAR_dengue <- IMCAR$summary.random$s[1:n, "mean"]
colours <- colorBin(palette = mycolors, domain = Output.Areas@data$IMCAR_dengue,
                    bins = seq(min(Output.Areas@data$IMCAR_dengue), max(Output.Areas@data$IMCAR_dengue), len=9))
(map_IMCAR_dengue_sp <- leaflet(data=Output.Areas, options = leafletOptions(zoomControl = FALSE)) %>% 
    addTiles() %>%
    addPolygons(fillColor = ~colours(IMCAR_dengue), color="white", weight=1,
                fillOpacity = 1) %>% addLegend(pal = colours, labFormat = labelFormat(digits = 1),
                                               values = Output.Areas@data$IMCAR_dengue,
                                               opacity = 1, title="MCAR dengue",
                                               position = "bottomright") %>%
    addScaleBar(position="bottomleft"))


Output.Areas@data$IMCAR_chiku <- IMCAR$summary.random$s[(n+1):(2*n), "mean"]
colours <- colorBin(palette = mycolors, domain = Output.Areas@data$IMCAR_chiku,
                    bins = seq(min(Output.Areas@data$IMCAR_chiku), max(Output.Areas@data$IMCAR_chiku), len=9))
(map_IMCAR_chiku_sp <- leaflet(data=Output.Areas, options = leafletOptions(zoomControl = FALSE)) %>% 
    addTiles() %>%
    addPolygons(fillColor = ~colours(IMCAR_chiku), color="white", weight=1,
                fillOpacity = 1) %>% addLegend(pal = colours, labFormat = labelFormat(digits = 2),
                                               values = Output.Areas@data$IMCAR_chiku,
                                               opacity = 1, title="MCAR chikungunya",
                                               position = "bottomright") %>%
    addScaleBar(position="bottomleft"))

# **************************************************************************** #
# Results: Relative Risk
# **************************************************************************** #
Output.Areas@data$IMCAR_dengue <- IMCAR$summary.fitted[1:n, "mean"]
colours <- colorBin(palette = mycolors, domain = Output.Areas@data$RR_d, 
                    bins = c(0, 0.3, 0.6, 0.975, 1.025, 1.5, 2,
                             max(Output.Areas@data$RR_d)))

(map_IMCAR_dengue <- leaflet(data=Output.Areas, options = leafletOptions(zoomControl = FALSE)) %>% 
    addTiles() %>%
    addPolygons(fillColor = ~colours(IMCAR_dengue), color="white", weight=1,
                fillOpacity = 1) %>% 
    addLegend(pal = colours, values = Output.Areas@data$IMCAR_dengue, opacity = 1, 
              labFormat = labelFormat(digits = 1), title = "Dengue", 
              position = "bottomright") %>%
    addScaleBar(position="bottomleft") %>% 
    setView(lng = -39, lat = -5.5, zoom = 7 ) %>%
    setMaxBounds( lng1 = -39.4, lat1 = -5.5, lng2 = -39.2, lat2 = -5))

Output.Areas@data$IMCAR_chiku <- IMCAR$summary.fitted[(n+1):(2*n), "mean"]
colours <- colorBin(palette = mycolors, domain = Output.Areas@data$RR_c, 
                    bins = c(0, 0.3, 0.6, 0.975, 1.025, 1.5, 2, 
                             max(Output.Areas@data$RR_c)))
(map_IMCAR_chiku <- leaflet(data=Output.Areas, options = leafletOptions(zoomControl = FALSE)) %>% 
    addTiles() %>%
    addPolygons(fillColor = ~colours(IMCAR_chiku), color="white", weight=1,
                fillOpacity = 1) %>% 
    addLegend(pal = colours, values = Output.Areas@data$IMCAR_chiku, opacity = 1, 
              labFormat = labelFormat(digits = 1), title = "Chikungunya",
              position = "bottomright") %>%
    addScaleBar(position="bottomleft") %>% 
    setView(lng = -39, lat = -5.5, zoom = 7 ) %>%
    setMaxBounds( lng1 = -39.4, lat1 = -5.5, lng2 = -39.2, lat2 = -5))