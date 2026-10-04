# Calcule les résultats lus par l'application et les enregistre dans data/processed/
# À lancer depuis la racine du projet : source("R/resultats.R")

source(here::here("R", "preparation.R"))
source(here::here("R", "analyse.R"))

base = construire_base()
carte = construire_carte(base)

# Clusters LISA : contiguïté, puis 5 plus proches voisins
listw = poids_queen(carte)
carte$lisa = clusters_lisa(carte, listw)
carte$lisa2 = clusters_lisa(carte, poids_knn(carte, k = 5))

# Coefficients locaux de la GWR
gwr = estimer_gwr(carte)
gwr_sf = st_as_sf(gwr$SDF)
carte$coef_gini = gwr_sf$gini
carte$coef_taux_pauvrete = gwr_sf$taux_pauvrete
carte$coef_densite = gwr_sf[["log(densite_pop)"]]
carte$coef_part_seniors = gwr_sf$part_seniors
carte$r2_local = gwr_sf$Local_R2

dir.create(here::here("data", "processed"), showWarnings = FALSE)
saveRDS(carte, here::here("data", "processed", "carte_resultats.rds"))
