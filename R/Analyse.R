# Poids spatiaux, clusters LISA et régression géographiquement pondérée
# Utilisé par le rapport et par le script de résultats.

pacman::p_load(sf, spdep, GWmodel)

# Modèle étudié
f = taux_creation ~ taux_pauvrete + gini + log(densite_pop) + part_seniors

# Matrice de poids par contiguïté (Queen)
poids_queen = function(carte) {
  mapC.nb = poly2nb(carte, queen = TRUE)
  nb2listw(mapC.nb, style = "W", zero.policy = TRUE)
}

# Matrice de poids par k plus proches voisins
poids_knn = function(carte, k = 5) {
  coords = st_coordinates(st_centroid(st_geometry(carte)))
  nb2listw(knn2nb(knearneigh(coords, k = k)), style = "W")
}

# Clusters LISA du taux de création pour une matrice de poids donnée
clusters_lisa = function(carte, listw) {
  # LISA : indice de Moran local
  lmoran = localmoran(carte$taux_creation, listw, zero.policy = TRUE)
  pval = lmoran[, ncol(lmoran)]

  moy = mean(carte$taux_creation)
  voisins = lag.listw(listw, carte$taux_creation, zero.policy = TRUE)

  # seuil à 0.10 ; les p-value manquantes sont traitées comme non significatives
  s = !is.na(pval) & pval <= 0.10

  lisa = rep("Non significatif", nrow(carte))
  lisa[s & carte$taux_creation >= moy & voisins >= moy] = "Haut-Haut"
  lisa[s & carte$taux_creation <= moy & voisins <= moy] = "Bas-Bas"
  lisa[s & carte$taux_creation >= moy & voisins <= moy] = "Haut-Bas"
  lisa[s & carte$taux_creation <= moy & voisins >= moy] = "Bas-Haut"
  factor(lisa, levels = c("Haut-Haut","Bas-Bas","Haut-Bas","Bas-Haut","Non significatif"))
}

# Régression géographiquement pondérée
estimer_gwr = function(carte) {
  carte_sp = as(carte, "Spatial")
  # Choix du nombre de voisins autour de chaque point
  bw = bw.gwr(f, data = carte_sp, approach = "AICc",kernel = "bisquare", adaptive = TRUE)
  # Estimation de GWR
  gwr = gwr.basic(f, data = carte_sp, bw = bw,kernel = "bisquare", adaptive = TRUE)
  gwr
}
