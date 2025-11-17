#!/usr/bin/env python3

from Bio import SeqIO
import gzip
from concurrent.futures import ProcessPoolExecutor
import subprocess
import os
import matplotlib.pyplot as plt
import networkx as nx

def open_fastq(path):
    """Ouvre un fichier FASTQ, qu'il soit gzippé ou non."""
    with open(path, 'rb') as test_f:
        magic_number = test_f.read(2)
    if magic_number == b'\x1f\x8b':
        return gzip.open(path, 'rt')
    else:
        return open(path, 'r')

def generate_kmers(sequence: str, k: int):
    """Génère tous les k-mers de taille k à partir d'une séquence donnée."""
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

#def process_record(record):
#    """Génère les k-mers pour un record FASTQ donné."""
#    return generate_kmers(str(record.seq), 3)


###############################################
# FONCTIONS SCALABLES POUR GROS FASTQ
###############################################

def build_kmer_database_to_file(fastq_path: str, k: int, output_file: str):
    """
    Génère un fichier contenant TOUS les k-mers extraits du FASTQ.
    Version scalable : aucun stockage en RAM.
    
    Args:
        fastq_path: Chemin du fichier FASTQ (.fastq ou .fastq.gz)
        k: Taille des k-mers
        output_file: Fichier de sortie où écrire les k-mers
    """
    print(f"[INFO] Génération des k-mers {k}-mer → {output_file}")
    
    with open(output_file, "w") as out:
        with open_fastq(fastq_path) as handle:
            for record in SeqIO.parse(handle, "fastq"):
                seq = str(record.seq).upper()

                # Génération en streaming
                for kmer in generate_kmers(seq, k):
                    out.write(kmer + "\n")

    print("[OK] Fichier brut de k-mers généré.")


def count_kmers_with_sort(kmer_file: str, output_count_file: str):
    """
    Utilise sort | uniq -c pour compter les k-mers sans utiliser Python.
    Version scalable et très rapide, même avec des milliards de k-mers.
    
    Args:
        kmer_file: Fichier contenant 1 k-mer par ligne
        output_count_file: Fichier résultant contenant "<count> <kmer>"
    """
    print("[INFO] Comptage des k-mers via sort | uniq -c ...")

    # sort et uniq -c permettent un comptage externe, scalable
    cmd = f"sort {kmer_file} | uniq -c > {output_count_file}"
    subprocess.run(cmd, shell=True, check=True)

    print(f"[OK] Comptage terminé → {output_count_file}")


def build_debruijn_edges_to_file_with_labels(kmer_file: str, output_edges_file: str):
    """
    Construit les arcs du graphe de De Bruijn avec étiquette = kmer.
    Ne supprime PAS les doublons : plusieurs occurrences = plusieurs arcs.

    Format : prefix \t suffix \t kmer
    """
    print(f"[INFO] Génération des arcs étiquetés → {output_edges_file}")
    
    with open(kmer_file) as km_in, open(output_edges_file, "w") as edges_out:
        for line in km_in:
            kmer = line.strip()
            if len(kmer) < 2:
                continue
            
            prefix = kmer[:-1]
            suffix = kmer[1:]

            # Format complet avec étiquette
            edges_out.write(f"{prefix}\t{suffix}\t{kmer}\n")
    
    print("[OK] Arcs étiquetés générés.")



def sort_debruijn_edges(edge_file: str, output_sorted_file: str):
    """
    Trie les arcs du graphe de De Bruijn par préfixe.
    Permet de regrouper les suffixes ensemble sans tout mettre en RAM.
    
    Args:
        edge_file: Fichier "prefix suffix" non trié
        output_sorted_file: Fichier trié par clé (prefix)
    """
    print("[INFO] Tri des arcs du graphe...")

    cmd = f"sort {edge_file} > {output_sorted_file}"
    subprocess.run(cmd, shell=True, check=True)

    print(f"[OK] Graphe trié généré → {output_sorted_file}")

###############################################
# VISUALISATION DU GRAPHE DE DE BRUIJN
###############################################


def visualize_graph(edge_file: str, max_edges=300):
    print("[INFO] Construction du graphe pour visualisation...")

    G = nx.MultiDiGraph()

    with open(edge_file) as f:
        for i, line in enumerate(f):
            if i > max_edges:
                break

            prefix, suffix, kmer = line.strip().split("\t")
            G.add_edge(prefix, suffix, label=kmer)

    pos = nx.spring_layout(G, seed=42)

    plt.figure(figsize=(12, 10))

    nx.draw(G, pos, with_labels=True,
            node_size=500, font_size=8,
            arrowsize=10, node_color="#66c2a5")

    # 🔥 Labels corrects pour MultiDiGraph
    edge_labels = {(u, v, k): d["label"]
                   for u, v, k, d in G.edges(keys=True, data=True)}

    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=6)

    plt.title("Graphe de De Bruijn (k-mers = étiquettes d'arêtes)")
    plt.tight_layout()
    plt.show()




#######################################
path = "exemple.fastq"

k = 3

# Étape A — Génération des k-mers (streaming)
build_kmer_database_to_file(path, k, "kmers_raw.txt")

# Étape B — Comptage (sort | uniq -c)
count_kmers_with_sort("kmers_raw.txt", "kmers_counted.txt")

# Étape C — Construction des arcs du graphe
build_debruijn_edges_to_file_with_labels("kmers_raw.txt", "edges_labeled.txt")


# Étape D — Tri final du graphe
#sort_debruijn_edges("edges_labeled.txt","edges_labeled_sorted.txt")

#visualize_graph("edges_labeled.txt")


#################################################################
# Charger les arcs depuis ton fichier
import pygraphviz as pgv
edges_file = "edges_labeled.txt"
G = pgv.AGraph(directed=True, strict=False)  # strict=False pour autoriser multi-arêtes

with open(edges_file) as f:
    for line in f:
        prefix, suffix, kmer = line.strip().split("\t")
        # Ajouter l'arête avec étiquette
        G.add_edge(prefix, suffix, label=kmer)

# Mise en page automatique avec Graphviz (dot, neato, fdp, etc.)
G.layout(prog='dot')  # 'dot' pour DAG, 'neato' ou 'fdp' pour plus libre

# Sauvegarder en PNG ou PDF
G.draw("debruijn_graph.png")
G.draw("debruijn_graph.pdf")

print("Graphe généré : debruijn_graph.png")