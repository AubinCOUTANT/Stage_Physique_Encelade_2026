# Modélisation des propriétés physiques de la surface glacée d'Encelade

Ce dépôt contient les ressources numériques développées dans le cadre de mon stage de Physique (2026) au **Laboratoire de Planétologie et Géosciences (LPG)** de Nantes.

## Présentation du Projet
L'objectif de ce travail est de modéliser l'évolution microstructurale du régolithe d'Encelade en simulant la compétition thermodynamique entre le **frittage** (*sintering*) et la **sublimation**. Ces simulations visent à contraindre la stabilité mécanique de la surface pour préparer l'atterrissage de la future mission **ESA-L4**.

## Index des Ressources

Les scripts ci-dessous sont indexés selon les mentions "Code source" figurant dans les légendes des figures du rapport de stage.


* **[Code 1] : Thermodynamique de la glace**
    * **Fichier** : `andreas_thermo.py`
    * **Description** : Implémentation de la loi d'Andreas (2007) pour le calcul de la pression de vapeur saturante ($P_{sat}$) sur la gamme de températures d'Encelade. Sert de dépendance aux autres modules.
    
* **[Code 2] : Modèle de Sublimation (Validation)**
    * **Fichier** : `sublimation_model.py`
    * **Description** : Résolution numérique (RK4) de l'érosion des grains selon la loi de Hertz-Knudsen et le modèle géométrique de Gundlach et al. (2018).
    * **Utilisation** : Reproduit les données expérimentales de référence pour la perte de masse (Figure 2 du rapport).

* **[Code 3] : Modèle de Sintering (Interaction dynamique)**
    * **Fichier** : `sintering_model.py`
    * **Description** : Algorithme utilisant un solveur **Runge-Kutta 4 vectoriel** pour résoudre le système d'équations différentielles couplées. Il simule l'interaction simultanée entre la croissance du pont de glace (*neck*) par diffusion et son érosion par sublimation.
    * **Utilisation** : Génère la courbe de validation montrant l'évolution temporelle du neck ainsi que les phases de "Neck Evolution" et "Solidification" (Figure 3 du rapport).

* **[Code 4] : Cartographie de la stabilité thermique (Zonage)**
    * **Fichier** : `enceladus_extrapolation.py`
    * **Description** : Application du modèle de frittage aux environnements thermiques spécifiques d'Encelade. Le code simule l'évolution sur des temps géologiques pour trois zones distinctes : les plaines inertes (80 K), les marges des *Tiger Stripes* (100-160 K) et les points chauds actifs (180-220 K).
    * **Utilisation** : Produit le graphique de zonage thermique évaluant la consolidation de la surface pour chaque région (Figure 4 du rapport).

* **[Code 5] : Sensibilité Granulométrique (Plaines - 80 K)**
    * **Fichier** : `enceladus_grain_size_effect_80K.py`
    * **Description** : Étude de l'influence de la granulométrie sur la cinétique de frittage à température fixée (80 K). Le code intègre la théorie **JKR** (Johnson-Kendall-Roberts) pour calculer l'adhésion initiale et simule l'évolution pour des grains de 0,5 à 100 $\mu m$.
    * **Utilisation** : Produit la Figure 6 du rapport, illustrant l'absence de consolidation par frittage dans les plaines.

* **[Code 6] : Sensibilité Granulométrique (Zones Actives - 120-160 K)**
    * **Fichier** : `enceladus_grain_size_effect_120-160K.py`
    * **Description** : Simulation de la cinétique de frittage dans les conditions thermiques des zones actives (120 K, 140 K, 160 K). Le code génère une figure multi-panneaux (disposition triangulaire) mettant en évidence la divergence des temps de consolidation selon la taille des grains et la température.
    * **Utilisation** : Produit la Figure 7 du rapport.
    * 
* **[Code 7] : Résistance à la Traction Statique**
    * **Fichier** : `enceladus_tensile_strength_synthesis.py`
    * **Description** : Calcul de la *Tensile Strength* ($Y$) théorique d'un empilement granulaire en fonction du rayon de particule et de la température, sans prise en compte du frittage dynamique. Compare deux scénarios de porosité ($\phi=0.1$ et $\phi=0.3$).
    * **Utilisation** : Produit la Figure 8 du rapport.

* **[Code 8] : Évolution Dynamique de la Résistance (Sintering)**
    * **Fichier** : `enceladus_tensile_strength_neck_evolution.py`
    * **Description** : Modélisation dynamique de l'augmentation de la cohésion induite par la croissance du pont de glace (*neck*) dans les zones actives "chaudes" (120-160 K). Utilise une visualisation en triangle pour comparer les régimes.
    * **Utilisation** : Produit la Figure 9 du rapport, montrant le renforcement mécanique rapide du sol.

* **[Code 9] : Approximation du Facteur Géométrique ($\Xi$)**
    * **Fichier** : `enceladus_xi_factor_approximation.py`
    * **Description** : Routine d'approximation numérique du facteur géométrique $\Xi$ (Xi) nécessaire au calcul de la conductivité thermique. Extrapole les données de référence de Chan & Tien (1973) aux milieux à haute porosité ($\phi=0,1$) via une loi de puissance.
    * **Utilisation** : Produit la Figure 10 du rapport.

* **[Code 10] : Conductivité Thermique Dynamique**
    * **Fichier** : `enceladus_thermal_conductivity_dynamic.py`
    * **Description** : Modélisation de l'évolution de la conductivité thermique effective du régolithe lors du processus de frittage à 120 K. Intègre le facteur $\Xi$ calculé précédemment.
    * **Utilisation** : Produit la Figure 11 du rapport.

---
**Auteur** : Aubin COUTANT (L3 Physique)  
**Encadrants** : Gabriel TOBIE, Riccardo ARTONI  
**Laboratoire** : LPG - Nantes Université / CNRS UMR 6112
