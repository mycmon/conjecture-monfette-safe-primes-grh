## Reliability

### 1. What the program measures robustly

- **Segmented vectorized sieve.**  
  Safe prime counting is based on a double segmented sieve (for p and (p−1)/2), tested up to x = 2·10⁹ and beyond.
- **Theoretical constant C₂.**  
  \( C_2 \) is computed as the Euler product  
  \( \displaystyle C_2 = \prod_{p\ge 3} \frac{p(p-2)}{(p-1)^2} \)  
  using a 1,000,000‑prime base, sufficient for the displayed precision.
- **Reproducible reports.**  
  Each run generates a Markdown report including:
  - parameters (x_max, ε, time),
  - numerical values (π_SG, li₂, C₂_emp, errors),
  - exported figures,
  - automatic analysis of the three plots.

### 2. What the program does *not* claim to prove

- It does **not** prove the Hardy–Littlewood conjecture.
- It does **not** prove the Generalized Riemann Hypothesis (GRH).
- It provides **numerical data compatible** with these conjectures, in particular:
  - the error stays far below the envelope \( \sqrt{x}\log^2 x \),
  - the normalized ratio |err| / (√x·log²x) is in the 10⁻³–10⁻² range for tested values.

### 3. Limits and transparency

- Results depend on double‑precision arithmetic and the libraries used (`numpy`, `scipy`).

- GRH-type envelopes are interpreted as **compatibility tests**, not proofs.

- The code is structured to be readable, modifiable, and verifiable by other researchers.

  -----

  ***Français***

  -----

  

## Fiabilité

### 1. Ce que le programme mesure de façon robuste

- **Crible segmenté vectorisé.**  
  Le comptage des premiers sûrs repose sur un double crible segmenté (p et (p−1)/2), testé jusqu’à x = 2·10⁹ et au‑delà.
- **Constante C₂ théorique.**  
  \( C_2 \) est calculée comme produit eulérien  
  \( \displaystyle C_2 = \prod_{p\ge 3} \frac{p(p-2)}{(p-1)^2} \)  
  avec une base de 1 000 000, suffisante pour la précision affichée.
- **Rapports reproductibles.**  
  Chaque exécution génère un rapport Markdown incluant :
  - les paramètres (x_max, ε, temps),
  - les valeurs numériques (π_SG, li₂, C₂_emp, erreurs),
  - les graphiques exportés,
  - une analyse automatique des trois figures.

### 2. Ce que le programme ne prétend pas démontrer

- Il **ne prouve pas** la conjecture de Hardy–Littlewood.
- Il **ne prouve pas** l’hypothèse de Riemann généralisée (GRH).
- Il fournit des **données numériques compatibles** avec ces conjectures, en particulier :
  - l’erreur reste très largement sous l’enveloppe \( \sqrt{x}\log^2 x \),
  - le ratio normalisé |err| / (√x·log²x) est de l’ordre de 10⁻³ à 10⁻² sur les plages testées.

### 3. Limites et transparence

- Les résultats dépendent de la précision numérique de `double` et des bibliothèques utilisées (`numpy`, `scipy`).
- Les bornes de type GRH sont interprétées comme **tests de compatibilité**, pas comme preuves.
- Le code est structuré pour être lisible, modifiable et vérifiable par d’autres chercheurs.

---

---

- 
