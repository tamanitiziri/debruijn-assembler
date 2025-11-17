#!/usr/bin/env python3

from Bio import SeqIO
import gzip
from concurrent.futures import ProcessPoolExecutor
import subprocess
import os
import matplotlib.pyplot as plt
import networkx as nx
from collections import defaultdict, deque
import pygraphviz as pgv

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


class EulerianPathFinder:
    def __init__(self, edges_file: str):
            """
            Initialise avec le fichier d'arcs étiquetés.
            """
            from collections import defaultdict, deque

            self.graph = defaultdict(list)  # sommet -> liste des arcs (suffix, kmer)
            self.arc_count = 0
            
            # Pour suivre les degrés
            self.in_degree = defaultdict(int)
            self.out_degree = defaultdict(int)

            # Charger le graphe
            self.load_graph(edges_file)
            
    def load_graph(self, edges_file: str):
        """Charge le graphe depuis le fichier d'arcs."""
        print("[INFO] Chargement du graphe...")
        with open(edges_file) as f:
            for line in f:
                prefix, suffix, kmer = line.strip().split("\t")
                self.graph[prefix].append((suffix, kmer))
                self.out_degree[prefix] += 1
                self.in_degree[suffix] += 1
                self.arc_count += 1
        print(f"[OK] Graphe chargé: {self.arc_count} arcs, {len(self.graph)} sommets")
    
    def find_start_node(self):
        """
        Trouve le nœud de départ pour le cycle eulérien.
        S'il y a un sommet avec out_degree > in_degree, c'est le début.
        Sinon, on prend un sommet arbitraire avec des arcs sortants.
        """
        start_node = None
        
        # Chercher un sommet avec plus d'arcs sortants qu'entrants
        all_nodes = set(self.out_degree.keys()) | set(self.in_degree.keys())
        
        for node in all_nodes:
            if self.out_degree[node] > self.in_degree[node]:
                return node
        
        # Si aucun trouvé, prendre un sommet avec des arcs sortants
        for node in self.graph:
            if self.graph[node]:  # a des arcs sortants
                return node
        
        # Fallback: premier sommet du graphe
        return next(iter(self.graph)) if self.graph else None
    
    def hierholzer_algorithm(self):
        """
        Implémentation de l'algorithme de Hierholzer pour trouver un cycle eulérien.
        Complexité: O(n) où n est le nombre d'arcs.
        
        Returns:
            Liste des k-mers dans l'ordre du cycle eulérien
        """
        if self.arc_count == 0:
            return []
        
        # Faire une copie du graphe pour pouvoir le modifier
        graph_copy = defaultdict(deque)
        for node in self.graph:
            graph_copy[node] = deque(self.graph[node])
        
        # Étape 1: Choisir un sommet de départ
        start_node = self.find_start_node()
        if start_node is None:
            raise ValueError("Aucun sommet de départ trouvé")
        
        print(f"[INFO] Début du cycle eulérien au sommet: {start_node}")
        
        # Stack pour le cycle principal
        stack = [start_node]
        cycle = []  # Stockera les k-mers dans l'ordre
        
        # Étape 2: Construire le cycle
        while stack:
            current_node = stack[-1]
            
            # Si le nœud courant a encore des arcs sortants
            if graph_copy[current_node]:
                # Prendre le premier arc disponible
                next_node, kmer_label = graph_copy[current_node].popleft()
                
                # Ajouter le nœud suivant à la pile
                stack.append(next_node)
                
                # Stocker le k-mer (optionnel, selon ce qu'on veut en sortie)
                # cycle.append(kmer_label)
                
            else:
                # Plus d'arcs sortants, ajouter au cycle et backtrack
                if len(stack) > 1:
                    # L'arc qu'on vient de parcourir en backtrack
                    from_node = stack[-2]
                    to_node = stack[-1]
                    
                    # Trouver le k-mer correspondant (pourrait être optimisé)
                    for i, (n, kmer) in enumerate(self.graph[from_node]):
                        if n == to_node:
                            cycle.append(kmer)
                            break
                
                stack.pop()
        
        # Inverser le cycle car on a construit à l'envers lors du backtrack
        cycle.reverse()
        
        # Vérifier qu'on a utilisé tous les arcs
        if len(cycle) != self.arc_count:
            print(f"[ATTENTION] Cycle incomplet: {len(cycle)}/{self.arc_count} arcs utilisés")
        else:
            print(f"[SUCCÈS] Cycle eulérien trouvé: {len(cycle)} arcs")
        
        return cycle
    
    def find_eulerian_path(self):
        """
        Trouve un chemin eulérien (si le graphe n'a pas de cycle eulérien).
        Ajoute un arc fictif si nécessaire.
        """
        # Identifier les nœuds avec déséquilibre
        start_candidates = []
        end_candidates = []
        
        all_nodes = set(self.out_degree.keys()) | set(self.in_degree.keys())
        
        for node in all_nodes:
            diff = self.out_degree[node] - self.in_degree[node]
            if diff == 1:
                start_candidates.append(node)
            elif diff == -1:
                end_candidates.append(node)
            elif diff != 0:
                print(f"[ATTENTION] Sommet {node} a un déséquilibre important: {diff}")
        
        # Si parfaitement équilibré, utiliser Hierholzer normal
        if not start_candidates and not end_candidates:
            print("[INFO] Graphe équilibré - recherche de cycle eulérien")
            return self.hierholzer_algorithm()
        
        # Sinon, ajouter un arc fictif pour créer un cycle
        if len(start_candidates) == 1 and len(end_candidates) == 1:
            start_node = start_candidates[0]
            end_node = end_candidates[0]
            fictive_kmer = f"FICTIVE_{start_node}_{end_node}"
            
            print(f"[INFO] Ajout d'un arc fictif: {end_node} -> {start_node}")
            
            # Ajouter l'arc fictif
            self.graph[end_node].append((start_node, fictive_kmer))
            self.out_degree[end_node] += 1
            self.in_degree[start_node] += 1
            self.arc_count += 1
            
            # Trouver le cycle eulérien
            cycle = self.hierholzer_algorithm()
            
            # Retirer l'arc fictif et réorganiser le cycle
            return self._remove_fictive_edge(cycle, fictive_kmer, start_node, end_node)
        else:
            raise ValueError("Graphe ne peut pas avoir de chemin eulérien")
    
    def _remove_fictive_edge(self, cycle, fictive_kmer, start_node, end_node):
        """
        Retire l'arc fictif et réorganise le cycle en chemin.
        """
        if fictive_kmer not in cycle:
            print("[ATTENTION] Arc fictif non trouvé dans le cycle")
            return cycle
        
        # Trouver la position de l'arc fictif
        fictive_index = cycle.index(fictive_kmer)
        
        # Réorganiser le cycle pour commencer au vrai début
        if fictive_index + 1 < len(cycle):
            new_cycle = cycle[fictive_index + 1:] + cycle[:fictive_index]
        else:
            new_cycle = cycle[:fictive_index]
        
        print(f"[INFO] Chemin eulérien réorganisé: début={start_node}, fin={end_node}")
        return new_cycle
    
    def reconstruct_sequence(self, kmer_path: list, k: int):
        """
        Reconstruit la séquence ADN à partir du chemin de k-mers.
        
        Args:
            kmer_path: Liste des k-mers dans l'ordre du chemin
            k: Taille des k-mers
            
        Returns:
            Séquence ADN reconstruite
        """
        if not kmer_path:
            return ""
        
        # Commencer avec le premier k-mer
        sequence = kmer_path[0]
        
        # Pour chaque k-mer suivant, ajouter seulement la dernière base
        for next_kmer in kmer_path[1:]:
            sequence += next_kmer[-1]
        
        return sequence


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

###############################################################
# Étape E — Calcul du chemin/cycle eulérien


euler_finder = EulerianPathFinder("edges_labeled.txt")

# Chercher un cycle eulérien si possible
kmer_path = euler_finder.find_eulerian_path()

# Reconstruire la séquence complète à partir du chemin
sequence_reconstructed = euler_finder.reconstruct_sequence(kmer_path, k)

print("Séquence reconstruite")
#garder la sequen construite dans un fichier assemblied_sequence.txt
#print(sequence_reconstructed)
with open("assembled_sequence.txt", "w") as f:
    f.write(sequence_reconstructed + "\n")



#################################################################visualisation du graphe de De Bruijn avec pygraphviz
# Charger les arcs depuis ton fichier

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
#G.draw("debruijn_graph.png")
G.draw("debruijn_graph.pdf")

print("Graphe généré : debruijn_graph.png")