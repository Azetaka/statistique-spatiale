# Application Shiny : créations d'entreprises de commerce par département
# Elle lit data/processed/carte_resultats.rds, produit par R/resultats.R,
# et réutilise la formule et les fonctions de R/analyse.R.

library(shiny)
library(bslib)
library(leaflet)
library(sf)
library(ggplot2)
library(spdep)
library(spatialreg)
library(lmtest)
library(DT)

# spdep définit aussi une fonction card() : on appelle donc celle de bslib par son nom complet.

source("R/analyse.R")

carte = readRDS("data/processed/carte_resultats.rds")   # Lambert-93 : sert aux calculs
carte_wgs = st_transform(carte, 4326)                    # WGS84 : sert à l'affichage


# ── Libellés ─────────────────────────────────────────────────────────────────

DONNEES = c(
  "Taux de création (pour 1 000 hab.)" = "taux_creation",
  "Revenu médian (€)" = "revenu_median",
  "Taux de pauvreté (%)" = "taux_pauvrete",
  "Indice de Gini" = "gini",
  "Part des 75 ans et plus (%)" = "part_seniors",
  "Densité (hab./km²)" = "densite_pop"
)

INDICATEURS = list(
  "Données" = DONNEES,
  "Clusters LISA du taux de création" = c(
    "Voisinage par contiguïté" = "lisa",
    "Voisinage des 5 plus proches" = "lisa2"
  ),
  "Effets locaux (GWR)" = c(
    "Effet du Gini" = "coef_gini",
    "Effet du taux de pauvreté" = "coef_taux_pauvrete",
    "Effet de la densité" = "coef_densite",
    "Effet de la part de seniors" = "coef_part_seniors",
    "R² local" = "r2_local"
  )
)

AIDE = c(
  "Données" = "Les départements sont répartis en cinq classes de même effectif.",
  "Clusters LISA du taux de création" = "Haut-Haut : un département où l'on crée beaucoup, entouré de voisins où l'on crée beaucoup. Bas-Bas : l'inverse. Seuil de significativité : 10 %.",
  "Effets locaux (GWR)" = "Coefficient estimé autour de chaque département : l'effet de la variable sur le taux de création n'est pas le même partout."
)

TERMES_MODELE = c(
  "(Intercept)" = "Constante",
  "taux_pauvrete" = "Taux de pauvreté",
  "gini" = "Indice de Gini",
  "log(densite_pop)" = "Densité (logarithme)",
  "part_seniors" = "Part des 75 ans et plus"
)

COULEURS_LISA = c("Haut-Haut" = "#bd0026", "Bas-Bas" = "#2c7bb6", "Haut-Bas" = "#feb24c",
                  "Bas-Haut" = "#addd8e", "Non significatif" = "#f0f0f0")
COULEURS_TERRITOIRE = c("Rural (peu dense)" = "#78c679", "Urbain (dense)" = "#f03b20")
ETIQUETTES_AXE = scales::label_number(decimal.mark = ",", big.mark = " ")


# ── Fonctions d'affichage ────────────────────────────────────────────────────

famille = function(var) names(INDICATEURS)[sapply(INDICATEURS, function(x) var %in% x)]
libelle = function(var) {
  tous = unlist(unname(INDICATEURS))
  names(tous)[tous == var]
}

# Nombre à la française, avec des décimales adaptées à l'ordre de grandeur
fmt = function(x) {
  dec = ifelse(abs(x) >= 100, 0, ifelse(abs(x) >= 10, 1, ifelse(abs(x) >= 1, 2, 3)))
  mapply(function(v, d) format(round(v, d), nsmall = d, big.mark = " ", decimal.mark = ","), x, dec)
}
fmt_p = function(p) ifelse(p < 0.001, "<\u00a00,001", format(round(p, 3), nsmall = 3, decimal.mark = ","))
etoiles = function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", "")))

# Carte choroplèthe : clusters si `valeurs` est un facteur, cinq classes sinon
carte_choroplethe = function(valeurs, titre, palette = "OrRd", inverser = FALSE) {
  if (is.factor(valeurs)) {
    niveaux = levels(droplevels(valeurs))
    couleurs = unname(COULEURS_LISA[as.character(valeurs)])
    etiquettes = sprintf("<b>%s</b><br>%s", carte_wgs$nom_dep, valeurs)
    legende = function(m) addLegend(m, "bottomleft", colors = unname(COULEURS_LISA[niveaux]),
                                    labels = niveaux, title = titre, opacity = 1)
  } else {
    bornes = unique(quantile(valeurs, probs = seq(0, 1, 0.2), na.rm = TRUE))
    pal = colorBin(palette, domain = valeurs, bins = bornes, reverse = inverser)
    couleurs = pal(valeurs)
    etiquettes = sprintf("<b>%s</b><br>%s : %s", carte_wgs$nom_dep, titre, fmt(valeurs))
    k = length(bornes) - 1
    legende = function(m) addLegend(m, "bottomleft", colors = pal((bornes[-1] + bornes[-(k + 1)]) / 2),
                                    labels = sprintf("%s à %s", fmt(bornes[-(k + 1)]), fmt(bornes[-1])),
                                    title = titre, opacity = 1)
  }
  leaflet(carte_wgs, options = leafletOptions(zoomSnap = 0.25, attributionControl = FALSE)) |>
    addPolygons(layerId = ~code_dep, fillColor = couleurs, fillOpacity = 1,
                color = "#9a9a9a", weight = 0.6, label = lapply(etiquettes, HTML),
                highlightOptions = highlightOptions(color = "#333333", weight = 2, bringToFront = TRUE)) |>
    legende()
}


# ── Modèles (estimés une fois, avec le voisinage par contiguïté du rapport) ──

listw_ref = poids_queen(carte)
ols = lm(f, data = carte)
slx = lmSLX(f, data = carte, listw = listw_ref, zero.policy = TRUE)
sar = lagsarlm(f, data = carte, listw = listw_ref, zero.policy = TRUE)
sem = errorsarlm(f, data = carte, listw = listw_ref, zero.policy = TRUE)
sdm = lagsarlm(f, data = carte, listw = listw_ref, Durbin = TRUE, zero.policy = TRUE)

bp = bptest(ols)
moran_residus = lm.morantest(ols, listw_ref, zero.policy = TRUE)
# Selon la version de spdep, les tests LM s'appellent lm.RStests ou lm.LMtests
tests_lm = if (exists("lm.RStests")) {
  lm.RStests(ols, listw_ref, test = "all", zero.policy = TRUE)
} else {
  lm.LMtests(ols, listw_ref, test = "all", zero.policy = TRUE)
}
# Versions robustes : RLMerr / RLMlag (anciens noms) ou adjRSerr / adjRSlag (nouveaux noms)
lm_err = tests_lm[[grep("^(RLM|adjRS)err$", names(tests_lm))]]
lm_lag = tests_lm[[grep("^(RLM|adjRS)lag$", names(tests_lm))]]

coefs = summary(ols)$coefficients
TABLE_COEFS = data.frame(
  Variable = unname(TERMES_MODELE[rownames(coefs)]),
  Coefficient = paste0(fmt(coefs[, 1]), " ", etoiles(coefs[, 4])),
  `Écart-type` = fmt(coefs[, 2]),
  `p-value` = fmt_p(coefs[, 4]),
  check.names = FALSE, row.names = NULL
)
names(TABLE_COEFS)[4] = "p\u2011value"   # trait d'union insécable

lecture = function(p, oui, non) ifelse(p < 0.05, oui, non)
TABLE_TESTS = data.frame(
  Test = c("Breusch-Pagan", "Moran sur les résidus",
           "LM robuste, variable décalée (ρ)", "LM robuste, erreur (λ)"),
  Statistique = fmt(c(bp$statistic, moran_residus$statistic, lm_lag$statistic, lm_err$statistic)),
  `p-value` = fmt_p(c(bp$p.value, moran_residus$p.value, lm_lag$p.value, lm_err$p.value)),
  Lecture = c(
    lecture(bp$p.value, "Variance des résidus non constante", "Variance des résidus constante"),
    lecture(moran_residus$p.value, "Résidus autocorrélés dans l'espace", "Pas d'autocorrélation spatiale des résidus"),
    lecture(lm_lag$p.value, "Dépendance spatiale détectée", "Pas de dépendance spatiale"),
    lecture(lm_err$p.value, "Dépendance spatiale détectée", "Pas de dépendance spatiale")
  ),
  check.names = FALSE, row.names = NULL
)
names(TABLE_TESTS)[3] = "p\u2011value"

aic = c("MCO" = AIC(ols), "SLX" = AIC(slx), "SAR" = AIC(sar), "SEM" = AIC(sem), "SDM" = AIC(sdm))
aic = sort(aic)
TABLE_AIC = data.frame(Modèle = names(aic), AIC = fmt(aic), check.names = FALSE, row.names = NULL)

spatial_justifie = moran_residus$p.value < 0.05 | lm_lag$p.value < 0.05 | lm_err$p.value < 0.05
CONCLUSION = if (spatial_justifie) {
  "Au moins un test signale une dépendance spatiale dans les résidus : un modèle spatial est à envisager."
} else {
  "Les résidus des MCO ne présentent pas d'autocorrélation spatiale et les tests LM robustes ne sont pas significatifs. Aucun modèle spatial n'est justifié : la ressemblance entre départements voisins s'explique par leurs caractéristiques communes, pas par un effet de voisinage."
}

residus = residuals(ols)


# ── Interface ────────────────────────────────────────────────────────────────

ui = page_navbar(
  title = "Créations d'entreprises de commerce par département (2024)",
  fillable = "Carte",
  header = tags$style(".leaflet-container { background: #ffffff; }"),

  nav_panel(
    "Carte",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        selectInput("indicateur", "Indicateur à cartographier", choices = INDICATEURS),
        textOutput("aide"),
        hr(),
        h6("Fiche du département"),
        uiOutput("fiche")
      ),
      bslib::card(full_screen = TRUE, leafletOutput("carte", height = "100%"))
    )
  ),

  nav_panel(
    "Autocorrélation spatiale",
    layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        width = 340,
        radioButtons("type_voisinage", "Définition du voisinage",
                     choices = c("Contiguïté (frontière commune)" = "queen",
                                 "k plus proches voisins" = "knn")),
        conditionalPanel(
          "input.type_voisinage == 'knn'",
          sliderInput("k", "Nombre de voisins (k)", min = 2, max = 10, value = 5, step = 1)
        ),
        p("L'indice de Moran mesure la ressemblance entre départements voisins. Positif et significatif : les départements où l'on crée beaucoup sont entourés de départements où l'on crée beaucoup."),
        p("Change le voisinage pour vérifier que le résultat ne dépend pas de ce choix.")
      ),
      layout_columns(
        value_box("Indice de Moran", textOutput("moran_i"), height = "120px"),
        value_box("p-value", textOutput("moran_p"), height = "120px"),
        value_box("Clusters Haut-Haut / Bas-Bas", textOutput("nb_clusters"), height = "120px"),
        fill = FALSE
      ),
      layout_columns(
        bslib::card(card_header("Nuage de Moran"), plotOutput("nuage_moran", height = "440px")),
        bslib::card(card_header("Clusters LISA (seuil de 10 %)"), leafletOutput("carte_lisa", height = "440px")),
        fill = FALSE
      )
    )
  ),

  nav_panel(
    "Modèle",
    layout_columns(
      bslib::card(
        card_header("Régression par les moindres carrés ordinaires"),
        p("Variable expliquée : taux de création pour 1 000 habitants."),
        tableOutput("table_coefs"),
        textOutput("qualite_modele"),
        p(class = "text-muted small", "Seuils de significativité : *** 1 %, ** 5 %, * 10 %.")
      ),
      bslib::card(
        card_header("Tests sur les résidus"),
        tableOutput("table_tests"),
        p(strong("Conclusion."), CONCLUSION)
      ),
      col_widths = c(6, 6), fill = FALSE
    ),
    layout_columns(
      bslib::card(
        card_header("Comparaison avec les modèles spatiaux"),
        p("Plus l'AIC est faible, meilleur est le compromis entre ajustement et complexité."),
        tableOutput("table_aic")
      ),
      bslib::card(
        card_header("Résidus des MCO"),
        p("En rouge, les départements où l'on crée plus que ce que prévoit le modèle. En bleu, moins."),
        leafletOutput("carte_residus", height = "380px")
      ),
      col_widths = c(4, 8), fill = FALSE
    )
  ),

  nav_panel(
    "Explorer",
    layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        width = 340,
        selectInput("var_x", "Axe horizontal", choices = DONNEES, selected = "gini"),
        selectInput("var_y", "Axe vertical", choices = DONNEES, selected = "taux_creation"),
        textOutput("correlation"),
        hr(),
        textOutput("dep_survol")
      ),
      bslib::card(
        fill = FALSE,
        card_header("Relation entre deux variables"),
        plotOutput("nuage", height = "420px", hover = hoverOpts("survol", delay = 80))
      ),
      bslib::card(
        fill = FALSE,
        card_header("Les 96 départements"),
        DTOutput("table_dep", fill = FALSE)
      )
    )
  )
)


# ── Serveur ──────────────────────────────────────────────────────────────────

server = function(input, output, session) {

  # Onglet Carte
  output$aide = renderText(AIDE[[famille(input$indicateur)]])

  output$carte = renderLeaflet({
    var = input$indicateur
    gwr_local = famille(var) == "Effets locaux (GWR)"
    carte_choroplethe(carte_wgs[[var]], libelle(var),
                      palette = if (gwr_local) "RdBu" else "OrRd", inverser = gwr_local)
  })

  output$fiche = renderUI({
    clic = input$carte_shape_click
    if (is.null(clic)) {
      return(p("Clique sur un département pour afficher ses chiffres."))
    }
    d = st_drop_geometry(carte[carte$code_dep == clic$id, ])
    rang = rank(-carte$taux_creation)[carte$code_dep == clic$id]
    tagList(
      strong(sprintf("%s (%s)", d$nom_dep, d$code_dep)),
      tags$table(class = "table table-sm",
                 tags$tr(tags$td("Taux de création"), tags$td(sprintf("%s pour 1 000 hab.", fmt(d$taux_creation)))),
                 tags$tr(tags$td("Rang"), tags$td(sprintf("%d sur %d", rang, nrow(carte)))),
                 tags$tr(tags$td("Créations en 2024"), tags$td(fmt(d$creations_commerce))),
                 tags$tr(tags$td("Population"), tags$td(fmt(d$population))),
                 tags$tr(tags$td("Revenu médian"), tags$td(sprintf("%s €", fmt(d$revenu_median)))),
                 tags$tr(tags$td("Taux de pauvreté"), tags$td(sprintf("%s %%", fmt(d$taux_pauvrete)))),
                 tags$tr(tags$td("Indice de Gini"), tags$td(fmt(d$gini))),
                 tags$tr(tags$td("Densité"), tags$td(sprintf("%s hab./km²", fmt(d$densite_pop)))),
                 tags$tr(tags$td("Cluster LISA"), tags$td(as.character(d$lisa)))
      )
    )
  })

  # Onglet Autocorrélation spatiale
  voisinage = reactive({
    if (input$type_voisinage == "queen") poids_queen(carte) else poids_knn(carte, k = input$k)
  })
  moran = reactive(moran.test(carte$taux_creation, voisinage(), zero.policy = TRUE))
  lisa = reactive(clusters_lisa(carte, voisinage()))

  output$moran_i = renderText(fmt(moran()$estimate[1]))
  output$moran_p = renderText(fmt_p(moran()$p.value))
  output$nb_clusters = renderText({
    n = table(lisa())
    sprintf("%d / %d", n[["Haut-Haut"]], n[["Bas-Bas"]])
  })

  output$nuage_moran = renderPlot({
    d = data.frame(nom_dep = carte$nom_dep, taux_creation = carte$taux_creation,
                   voisins = lag.listw(voisinage(), carte$taux_creation, zero.policy = TRUE))
    mx = mean(d$taux_creation); my = mean(d$voisins)
    d$influent = abs(d$taux_creation - mx) > 1.5 * sd(d$taux_creation) | abs(d$voisins - my) > 1.5 * sd(d$voisins)
    ggplot(d, aes(taux_creation, voisins)) +
      geom_vline(xintercept = mx, color = "grey70", linetype = "dashed") +
      geom_hline(yintercept = my, color = "grey70", linetype = "dashed") +
      geom_point(color = "#f03b20", alpha = 0.7, size = 2.5) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.8) +
      ggrepel::geom_text_repel(data = subset(d, influent), aes(label = nom_dep), size = 4,
                               color = "grey20", max.overlaps = 20) +
      scale_x_continuous(labels = ETIQUETTES_AXE) + scale_y_continuous(labels = ETIQUETTES_AXE) +
      labs(x = "Taux de création du département", y = "Taux de création moyen des voisins") +
      theme_minimal(base_size = 15)
  })

  output$carte_lisa = renderLeaflet(carte_choroplethe(lisa(), "Clusters LISA"))

  # Onglet Modèle
  output$table_coefs = renderTable(TABLE_COEFS, align = "lrrr")
  output$qualite_modele = renderText(sprintf(
    "R² = %s, R² ajusté = %s, %d départements.",
    fmt(summary(ols)$r.squared), fmt(summary(ols)$adj.r.squared), nrow(carte)
  ))
  output$table_tests = renderTable(TABLE_TESTS, align = "lrrl")
  output$table_aic = renderTable(TABLE_AIC, align = "lr")
  output$carte_residus = renderLeaflet(
    carte_choroplethe(residus, "Résidu", palette = "RdBu", inverser = TRUE)
  )

  # Onglet Explorer
  donnees = st_drop_geometry(carte)

  output$nuage = renderPlot({
    ggplot(donnees, aes(.data[[input$var_x]], .data[[input$var_y]], color = type_terr)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(values = COULEURS_TERRITOIRE, name = NULL) +
      scale_x_continuous(labels = ETIQUETTES_AXE) + scale_y_continuous(labels = ETIQUETTES_AXE) +
      labs(x = libelle(input$var_x), y = libelle(input$var_y)) +
      theme_minimal(base_size = 15) +
      theme(legend.position = "top")
  })

  output$correlation = renderText(sprintf(
    "Corrélation entre les deux variables : %s",
    format(round(cor(donnees[[input$var_x]], donnees[[input$var_y]]), 2), nsmall = 2, decimal.mark = ",")
  ))

  output$dep_survol = renderText({
    proche = nearPoints(donnees, input$survol, xvar = input$var_x, yvar = input$var_y,
                        maxpoints = 1, threshold = 12)
    if (nrow(proche) == 0) {
      return("Survole un point pour lire le nom du département.")
    }
    sprintf("%s : %s en abscisse, %s en ordonnée.", proche$nom_dep,
            fmt(proche[[input$var_x]]), fmt(proche[[input$var_y]]))
  })

  output$table_dep = renderDT({
    t = donnees[, c("nom_dep", "taux_creation", "revenu_median", "taux_pauvrete", "gini",
                    "part_seniors", "densite_pop", "lisa")]
    datatable(
      t, rownames = FALSE, fillContainer = FALSE,
      colnames = c("Département", "Taux de création", "Revenu médian", "Taux de pauvreté",
                   "Gini", "Part des 75 ans et plus", "Densité", "Cluster LISA"),
      options = list(
        pageLength = 10, order = list(list(1, "desc")),
        columnDefs = list(list(className = "dt-right text-end", targets = 1:6)),
        language = list(search = "Rechercher :", lengthMenu = "Afficher _MENU_ lignes",
                        info = "_START_ à _END_ sur _TOTAL_ départements",
                        paginate = list(previous = "Précédent", `next` = "Suivant"))
      )
    ) |>
      formatRound(c("taux_creation", "taux_pauvrete", "part_seniors"), 2, dec.mark = ",", mark = " ") |>
      formatRound("gini", 3, dec.mark = ",", mark = " ") |>
      formatRound(c("revenu_median", "densite_pop"), 0, dec.mark = ",", mark = " ")
  })
}

shinyApp(ui, server)
