# Plan sur 6 mois

* Collaborations
  * Papier Déconvolution ?
  * Papier Lucille ?
  * Papier Margherita/RNA splicing ?


# Projet de thèse
* **Supervision** : Faire des réunions toutes les deux semaines en gardant la flexibilité de communiquer en direct ou par email
* **Objectif** : créer un multiplex network de la réponse aux IDHi
* **Plan pour les 6 prochains mois**

## I. Gene regulatory network  
* [ ] Travailler sur les algorithmes de gene regulatory network  
* [ ] Les comparer en utilisant notamment la prédiction de l'activité des TFs.
  * [ ] + Analyses de réseau (Eigenvalue, Pagerank etc.)

## II. Réseau 3D de chromatine
* [ ] Générer le réseau de chromatine des données CD34+ avec l'algo chicago. (Les data sont en cours de téléchargement).
* [ ] Continuer mes analyses d'assortativité de la chromatine pour définir quel réseau de chromatine utiliser pour mes analyses (un réseau par groupe de patient).
  * [ ] Selon les résultats, faire des analyses non supervisées en formant des groupes de patients en fonction des résultats d'assortativité/score leukemia stemness.
* [ ] Faire des réseaux avec la méthylation en feature et des réseaux avec l'expression des gènes + activité TF en feature.
  * [ ] + Analyses de réseau (Eigenvalue, Pagerank etc.)

## III. Réseau réactions métaboliques
* [ ] Prédire l'état des réactions sur le jeu de données d'expressions génique (en cours)
* [ ] Calculer les réactions significativement différentes entre les différents groupes de patients.
  * [ ] Selon les résultats, faire des analyses non supervisées en formant des groupes de patients en fonction des profiles de réaction (les comparer avec des analyses de clustering d'expression de gène brute)
* [ ] Faire un réseau de ces réactions en ne gardant que celles différemment actives
  * [ ] + Analyses de réseau (Eigenvalue, Pagerank etc.)

## IV. Open chromatin access
* [ ] Analyser les ATAC-seq
* [ ] Tester l'assortativité et voir si ça se rapproche avec les résultats d'assortativité de l'expression des gènes et de la méthylation des données de lignées cellulaires
* [ ] Affiner le Gene regulatory network au niveau des TFs

## V. Multiplex Network Analysis
* [ ] Faire la biblio des analyses de multiplex possible
* [ ] Analyser le / les différents multiplex

# In vitro

## Mai/juin 2022
* [ ] Summary all biological data
* [ ] define additional in vitro experiments (ATACseq, proteomics, CYTOF with IDH cell lines +/-IDHi treatment)

* [ ] Poster pour Jobim

* [ ] Identifier tous les IDH mutants avec les kmer dans beataml, Verhaak, TCGA : approches nathaniel et GSEA sur les stem cells ?
* [x] Etudier le papier de John Dick sur les stem cell et voir comment construire cet aspects dans les IDHmut versus wt et boucler avec le papier de Koichi Nature Comm
* [ ] MYC-STAT1-3-RELA : Correlation avec LSC17 score ?
* [ ] Lister tous les databases DNA methylome de patients LAM IDHmt et mutant : TCGA, Koichi, Figueroa
* [ ] deconvolution: blastes versus normales ?

## Juillet/Aout 2022
* [ ] Travailler sur nos data RNAseq de nos 200 patients avec données IDH
* [ ] Les TF en commun et différents entres les différents algo/models testés par Alexis
* [ ] Analysis cell type trajectory/psuedotime en methylation with prediction to IDHi response?

## Septembre/Octobre 2022
* [ ] Faire les figures rassemblant toutes les data (in silico et biological) en format article
* [ ] puis les légendes
* [ ] puis travailler sur les take-home messages à mettre en avant !
* [ ] définir le journal

## Novembre 2022 
* [ ] Comité de thèse à planifier

# Question/partie bio pour le papier à réfléchir
* [ ] validation epitranscriptomics coupled to CYTOF avec Lucille et Montpellier
* [ ] cell line regulatory network; est ce que nous faisons un ATACseq et methylome ?
* [ ] Definition de la leukemic stem cell des IDH mutants (lien avec le papier de Koichi et Dick)
* [ ] TET2/vitC: stem cell (lien avec le papier de Koichi et possible approche therapeutique et validation in vivo du papier d’alexis)
* [ ] epigenetic drives metabolic/catabolic flexibility: challenges galactose, glutamine starvation (lien avec nouveau papier Lucille sur la catabolic flexibility)
