## 1. Cryptography (practical impact)

- **Safe primes and cryptographic groups.**  
  The program studies safe primes p = 2q + 1, used in Diffie–Hellman (RFC 7919), SRP, ElGamal, and other protocols.
- **Generation speed-up.**  
  The primorial sieve based on the p–e law reduces the number of tested candidates to ~6.2% of integers, giving a **16× speed-up** when generating 2048‑bit safe primes.
- **Standard sizing.**  
  The density of n‑bit safe primes is  
  \( \displaystyle \frac{C_2}{(n \log 2)^2} \) with \( C_2 = 0.6601619 \),  
  allowing precise estimation of generation cost for standard sizes (2048, 3072, 4096 bits).

## 2. Number theory (conceptual contribution)

- **Combinatorial justification of local factors.**  
  The p–e law gives an exact identity for admissible residues modulo primorials, providing a rigorous combinatorial basis for the Hardy–Littlewood local factor \((1 - k/p)\).
  
- **Logarithmic correction.**  
  Numerical data suggests a correction of the form  
  \( \pi_{SG}(x) \approx C_2 \,\mathrm{li}_2(x)\,(1 + 3.2/\log x) \),  
  consistent with expected asymptotic expansions.
  
- **Visualization tool.**  
  The interface allows numerical exploration of the convergence of \( C_2^{emp}(x) \) toward \( C_2 \) and the behavior of the normalized error.
  
  ----------------
  
  ***Français***
  
  ----
  
  

### 1. Cryptographie (impact pratique)

- **Premiers sûrs et groupes cryptographiques.**  
  Le programme étudie les premiers sûrs p = 2q+1, utilisés dans Diffie–Hellman (RFC 7919), SRP, ElGamal, etc.
- **Accélération de la génération.**  
  Le crible primorial basé sur la loi p–e réduit le nombre de candidats testés à ~6.2 % des entiers, soit un gain d’environ **16×** pour la génération de premiers sûrs de 2048 bits.
- **Dimensionnement des standards.**  
  La densité des premiers sûrs de n bits est donnée par  
  \( \displaystyle \frac{C_2}{(n \log 2)^2} \) avec \( C_2 = 0.6601619 \), ce qui permet d’estimer précisément le coût de génération pour les tailles standard (2048, 3072, 4096 bits).

### 2. Théorie des nombres (apport conceptuel)

- **Justification combinatoire des facteurs locaux.**  
  La loi p–e donne une identité exacte sur les résidus admissibles modulo les primoriaux, ce qui fournit une base combinatoire rigoureuse au facteur local \((1 - k/p)\) de Hardy–Littlewood.
- **Correction logarithmique.**  
  Les données numériques suggèrent une correction de type  
  \( \pi_{SG}(x) \approx C_2 \,\mathrm{li}_2(x)\,(1 + 3.2/\log x) \),  
  cohérente avec les développements asymptotiques attendus.
- **Outil de visualisation.**  
  L’interface permet d’explorer numériquement la convergence de \( C_2^{emp}(x) \) vers \( C_2 \) et le comportement de l’erreur normalisée.