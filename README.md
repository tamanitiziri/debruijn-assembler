# De Bruijn Graph Assembler

Un script Python pour assembler des séquences à partir de fichiers FASTQ en utilisant des **graphes de De Bruijn** et l'algorithme de **Hierholzer** pour trouver des chemins eulériens.

---

## Description

Ce projet permet de :

1. Générer les **k-mers** à partir de séquences FASTQ (supporte fichiers `.fastq` et `.fastq.gz`).
2. Construire un **graphe de De Bruijn** avec des arcs étiquetés par les k-mers.
3. Trouver un **chemin ou cycle eulérien** pour reconstruire la séquence originale.
4. Visualiser le graphe de De Bruijn via **PyGraphviz**.
5. Générer des fichiers de sortie :
   - `kmers_raw.txt` : tous les k-mers extraits.
   - `kmers_counted.txt` : k-mers comptés (`sort | uniq -c`).
   - `edges_labeled.txt` : arcs du graphe avec k-mers comme étiquettes.
   - `assembled_sequence.txt` : séquence reconstruite.
   - `debruijn_graph.pdf` : visualisation du graphe.

---

## Prérequis

- Python ≥ 3.8
- Modules Python :
  ```bash
  pip install biopython matplotlib networkx pygraphviz
