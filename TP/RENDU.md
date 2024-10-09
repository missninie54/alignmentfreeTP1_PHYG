Rapport TME_2 PHYG
Eugénie GENESTANT
M2 BIM BMC

L'objectif de ce TP est de calculer des indices de Jaccard entre les paires d'organismes procaryotes pour retrouver les familles présentes.

Pour cela, différents fichiers sont fournis au sein du répértoire TP.

____

__init__.py : fichier d'initialisation du module.

loading.py : comprend la fonction de chargement des fichiers

kmers.py : contient les fonctions pour manipuler et comparer les k-mers.

matrice_jaccard_distance.ipynb : comprend le code pour obtenir la matrice de distance et le dendrogramme construit à partir des indices de jaccard (modules numpy, matplotlib, scipy et seaborn ont été utilisés)

__main__.py : comprend l'exécution du chargement du dossier "data" contenant les fichiers au format fasta téléchargés au préalable sur la base de données European Nucleotide Archive (ENA). Les fichiers fasta peuvent être compréssés ou non. Pour l'analyses, les 5 fichiers fasta.gz téléchargés manuellement sont les suivants :

GCA_000013265.1.fasta.gz
GCA_000008865.2.fasta.gz
GCA_030271835.1.fasta.gz
GCA_000005845.2.fasta.gz
GCA_000069965.1.fasta.gz

Note : Les fichiers fasta peuvent être chargés qu'ils soient compréssés ou non.

__

Méthodes implémentées dans le fichier kmers.py :

- 'encode_nucl' : convertit un nucléotide en un nombre entier.
Chaque nucléotide est associé à un entier :
'A' -> 0 ; 'C' -> 1 ; 'T' -> 2 ; 'G' -> 3
Cela permet de représenter les séquences d'ADN sous une forme numérique, facilitant les manipulations informatiques.

- 'encode_kmer' : convertit une séquence de k nucléotides (k-mer) en un entier.
Fonctionnement :
    La fonction initialise un entier kmer à 0.
    Pour chaque nucléotide dans la séquence de longueur k :
        Elle décale les bits de kmer de 2 positions vers la gauche (pour faire de la place pour le nouveau nucléotide).
        Elle ajoute la valeur numérique du nucléotide (obtenue via encode_nucl) aux bits de droite de kmer.
    Le résultat est un entier unique qui représente le k-mer d'origine.


- 'enumerate_kmers' : générer tous les k-mers d'une séquence donnée.
Fonctionnement :
    La fonction commence par calculer le k-mer initial (les k premiers nucléotides de la séquence) et le retourne (yield).
    Elle utilise un masque (1 << (2 * k)) - 1 pour s'assurer que seules les 2k bits les plus significatifs sont conservés lors du décalage des bits.
    Pour chaque position suivante dans la séquence :
        Elle décale le k-mer précédent de 2 bits à gauche.
        Elle applique le masque pour conserver les 2k bits les plus significatifs.
        Elle ajoute la valeur numérique du nouveau nucléotide (obtenue via encode_nucl) à la fin.
    Elle retourne chaque k-mer généré à partir de la séquence initiale.


- 'jaccard' : calcule l'indice de Jaccard entre deux ensembles de séquences.
Fonctionnement :
    La fonction prend deux listes de séquences (sequencesA et sequencesB) et la longueur des k-mers (k).
    Pour chaque séquence dans sequencesA et sequencesB :
        Elle génère tous les k-mers (à l'aide de enumerate_kmers) et les ajoute aux ensembles kmersA et kmersB respectivement.
    Elle calcule l'intersection (intersect) des deux ensembles de k-mers (les k-mers communs).
    Elle calcule l'union (union) des deux ensembles de k-mers (tous les k-mers uniques).
    Elle imprime la taille des k-mers en commun et la taille totale des k-mers uniques.
    L'indice de Jaccard est calculé comme le ratio de la taille de l'intersection sur la taille de l'union.
    Formule : J(A,B) = ∣A∩B∣ / ∣AUB∣​
        ∣A∩B∣ est le nombre de k-mers communs entre les deux sets. (notée intersection)
        ∣A∪B∣ est le nombre total de k-mers uniques entre les deux sets. (notée union)
    L'indice de Jaccard mesure la similarité entre les deux ensembles de k-mers, variant entre 0 (aucune similarité) et 1 (identiques).
  

Instructions
Pour exécuter le script principal et calculer les indices de Jaccard, il faut utiliser la commande suivante :
```sh
python3 -m TP
```

Dans le jupyter notebook : 

Production d'une matrice de distance à partir des indices de Jaccard obtenus.
La matrice de Jaccard fournie contient les indices de similarité entre les paires d'organismes suivants.


GCA_000013265.1 (Escherichia coli UTI89)
GCA_000008865.2 (Escherichia coli O157)
GCA_030271835.1 (Proteus appendicitidis)
GCA_000005845.2 (Escherichia coli K-12 MG1655)
GCA_000069965.1 (Proteus mirabilis HI4320)

Matrice de distance :

    GCA_000013265.1	GCA_000008865.2	GCA_030271835.1	GCA_000005845.2	GCA_000069965.1
GCA_000013265.1	0.0	0.3048	0.0008	0.3361	0.0008
GCA_000008865.2	0.3048	0.0	0.0008	0.446	0.0008
GCA_030271835.1	0.0008	0.0008	0.0	0.0009	0.0262
GCA_000005845.2	0.3361	0.446	0.0009	0.0	0.0009
GCA_000069965.1	0.0008	0.0008	0.0262	0.0009	0.0

![alt text](image.png) (lien vers la matrice lue avec Excel)

Chaque valeur dans la matrice représente l'indice de Jaccard entre deux organismes, c'est-à-dire la proportion de gènes partagés par rapport au total de gènes uniques et partagés.


Comparaisons Entre les Organismes : Résumé des Relations

    Escherichia coli UTI89 (GCA_000013265.1) et Escherichia coli O157
    (GCA_000008865.2) :
        Similarité modérée (0.3048), confirmant qu'ils appartiennent à la même espèce mais sont des souches différentes.

    Escherichia coli K-12 MG1655 (GCA_000005845.2) et les autres souches d'Escherichia coli :
        Similarité élevée avec E. coli O157
        (0.4460) et modérée avec E. coli UTI89 (0.3361), ce qui reflète la variabilité génétique intra-spécifique.

    Proteus appendicitidis (GCA_030271835.1) et Proteus mirabilis (GCA_000069965.1) :
        Faible similarité (0.0262), indiquant des différences génétiques importantes mais avec quelques gènes communs.

    Comparaisons entre Escherichia coli et Proteus spp. :
        Très faible similarité (0.0008 à 0.0009), comme attendu pour des espèces différentes.



Production d'un dendrogramme 'matrice_jaccard_distance.ipynb' : 

![alt text](image-1.png)

Aide pour comprendre les relations de similarité ou de distance entre plusieurs organismes. En utilisant les indices de Jaccard, qui mesurent la similarité entre les ensembles de gènes présents dans chaque organisme, nous pouvons visualiser comment ces organismes sont liés les uns aux autres.

    Clusters Principaux :
        Escherichia coli UTI89 (GCA_000013265.1) et Escherichia coli O157
        (GCA_000008865.2) sont proches l'un de l'autre avec un indice de Jaccard de 0.3048, indiquant une similarité génomique significative.
        Proteus appendicitidis (GCA_030271835.1) et Proteus mirabilis HI4320 (GCA_000069965.1) montrent également une certaine proximité, avec des indices de Jaccard relativement faibles (0.0262 entre eux), indiquant des différences génomiques importantes mais avec des éléments communs.


Conclusion :

1. Identification des Souches d'Escherichia coli
    Les deux souches d'Escherichia coli (GCA_000013265.1 et GCA_000008865.2) ainsi que la souche K-12 (GCA_000005845.2) montrent des indices de Jaccard modérés à élevés, indiquant qu'elles partagent une proportion de gènes. Ces souches appartiennent toutes à la même famille, Enterobacteriaceae, qui est connue pour regrouper divers membres du genre Escherichia.

2. Comparaison avec Proteus spp.
    Les indices de Jaccard entre les souches d'Escherichia coli et les organismes du genre Proteus (GCA_030271835.1 et GCA_000069965.1) sont très faibles (proches de 0). Les organismes Proteus sont également de la famille Enterobacteriaceae, mais leurs différences génétiques significatives suggèrent qu'ils appartiennent à des sous-genres ou des groupes phylogénétiquement distincts à l'intérieur de cette famille.

3. Clusters Phylogénétiques
    Dans le dendrogramme, les souches d'Escherichia coli se regroupent ensemble, tandis que les souches de Proteus se regroupent séparément. Cela confirme la proximité génétique entre les souches d'Escherichia coli et la divergence avec Proteus. Les souches d'Escherichia coli forment un cluster distinct, indiquant leur similarité, alors que les Proteus sont regroupés en un autre cluster, soulignant leur divergence.


