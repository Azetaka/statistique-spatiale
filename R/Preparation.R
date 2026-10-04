# Construction de la base départementale et de la carte
# Utilisé par le rapport et par l'application.

pacman::p_load(readxl, dplyr, readr, stringr, sf)

# Base : une ligne par département de France métropolitaine (96 lignes)
construire_base = function() {
  F_CREA = here::here("data", "raw", "Creation_Entreprise.xlsx")
  F_POP = here::here("data", "raw", "Population.xlsx")
  F_REV = here::here("data", "raw", "DISP_DEP.csv")
  F_PAUV = here::here("data", "raw", "FILO2021_DISP_PAUVRES_DEP.csv")

  # 1. Créations d'entreprises du commerce en 2024 (feuille DEP_UL)
  crea = read_excel(F_CREA, sheet = "DEP_UL", skip = 5, col_names = FALSE)
  names(crea) = c("code_dep","nom_dep","activite","forme_legale", paste0("an_", 2012:2024))
  crea = crea |>
    filter(str_starts(activite, "Commerce"), forme_legale == "Total") |>
    transmute(code_dep = as.character(code_dep), nom_dep, creations_commerce = an_2024)

  # Population 2024 : total + classes d'âge (bloc "Ensemble")
  pop = read_excel(F_POP, sheet = "2024", skip = 5, col_names = FALSE) |>
    transmute(code_dep = as.character(...1),
              pop_0_19 = as.numeric(...3), pop_75p = as.numeric(...7),
              population = as.numeric(...8))

  # Filosofi : revenu médian (Q221) + Gini (GI21)
  rev = read_delim(F_REV, delim = ";", locale = locale(decimal_mark = ","),
                   col_select = c(CODGEO, Q221, GI21),
                   col_types = cols(CODGEO = col_character(), Q221 = col_double(), GI21 = col_double())) |>
    transmute(code_dep = CODGEO, revenu_median = Q221, gini = GI21)

  # --- 4. Filosofi : taux de pauvreté 60 % (TP6021) ---
  pauv = read_delim(F_PAUV, delim = ";", locale = locale(decimal_mark = ","),
                    col_select = c(CODGEO, TP6021),
                    col_types = cols(CODGEO = col_character(), TP6021 = col_double())) |>
    transmute(code_dep = CODGEO, taux_pauvrete = TP6021)

  # base finale
  base = crea |>
    full_join(pop, by = "code_dep") |> full_join(rev, by = "code_dep") |>
    full_join(pauv, by = "code_dep") |>
    filter(str_detect(code_dep, "^(\\d{2}|2A|2B)$")) |>
    mutate(taux_creation = 1000 * creations_commerce / population,
           part_jeunes  = 100 * pop_0_19 / population,
           part_seniors = 100 * pop_75p  / population) |>
    select(code_dep, nom_dep, creations_commerce, population, taux_creation,
           revenu_median, taux_pauvrete, gini, part_jeunes, part_seniors) |>
    arrange(code_dep)

  stopifnot(nrow(base) == 96, all(!is.na(base$taux_creation)))
  base
}

# Carte : la base jointe au fond de carte, avec la densité et le type de territoire
construire_carte = function(base) {
  URL_GEO = here::here("data", "raw", "departements.geojson")
  dep_geo = st_read(URL_GEO, quiet = TRUE) |>
    st_transform(2154) |>
    filter(!str_starts(code, "97")) |>
    rename(code_dep = code)

  carte = dep_geo |>
    left_join(base, by = "code_dep")

  carte$densite_pop = carte$population / (as.numeric(st_area(carte)) / 1e6)
  carte$type_terr = ifelse(carte$densite_pop >= median(carte$densite_pop, na.rm = TRUE),"Urbain (dense)", "Rural (peu dense)")
  carte
}
