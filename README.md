# Modélisation des propriétés physiques de la surface glacée d'Encelade

Ce dépôt contient les ressources numériques développées dans le cadre de mon stage de Physique (2026) au **Laboratoire de Planétologie et Géosciences (LPG)** de Nantes.

## Présentation du Projet
L'objectif de ce travail est de modéliser l'évolution microstructurale du régolithe d'Encelade en simulant la compétition thermodynamique entre le **frittage** (*sintering*) et la **sublimation**. Ces simulations visent à contraindre la stabilité mécanique de la surface pour préparer l'atterrissage de la future mission **ESA-E4**.

## Index des Ressources (référencées dans le rapport)

Les scripts ci-dessous sont indexés selon les mentions "Code source" figurant dans les légendes des figures de mon rapport.

### [Code 1] : Thermodynamique de la glace
* **Fichier** : `andreas_thermo.py`
* **Description** : Implémentation de la loi d'Andreas (2007) pour le calcul de la pression de vapeur saturante ($P_{sat}$).
* **Utilisation** : Ce code a servi à générer la Figure 2 du rapport (Évolution de $P_{sat}$ sur la gamme 60-180 K).

### [Code 2] : Modèle de Sublimation (Validation)
* **Fichier** : `sublimation_model.py`
* **Description** : Résolution numérique (RK4) de l'érosion des grains selon la loi de Hertz-Knudsen et le modèle de Gundlach et al. (2018).
* **Utilisation** : Ce code permet de reproduire les données expérimentales de référence (Validation de la dynamique de sublimation présentée en Section 2.3).

---
**Auteur** : Aubin COUTANT (L3 Physique)
**Encadrants** : Gabriel TOBIE, Riccardo ARTONI
**Laboratoire** : LPG - Nantes Université / CNRS UMR 6112
