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

class EulerianPathFinder:
    def __init__(self, edges_file: str, k: int):
        """
        Initialise avec le fichier d'arcs étiquetés.
        """
        self.graph = defaultdict(list)  # sommet -> liste des arcs (suffix, kmer)
        self.arc_count = 0
        self.k = k
        
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
                parts = line.strip().split("\t")
                if len(parts) == 3:
                    prefix, suffix, kmer = parts
                    self.graph[prefix].append((suffix, kmer))
                    self.out_degree[prefix] += 1
                    self.in_degree[suffix] += 1
                    self.arc_count += 1
        print(f"[OK] Graphe chargé: {self.arc_count} arcs, {len(self.graph)} sommets")
    
    def find_start_node(self):
        """
        Trouve le nœud de départ pour le chemin eulérien.
        Dans un vrai assemblage, on prend souvent le nœud avec le plus grand excédent sortant.
        """
        start_node = None
        max_diff = -float('inf')
        
        all_nodes = set(self.out_degree.keys()) | set(self.in_degree.keys())
        
        for node in all_nodes:
            diff = self.out_degree[node] - self.in_degree[node]
            if diff > max_diff:
                max_diff = diff
                start_node = node
        
        # Si aucun nœud avec excédent sortant, prendre celui avec le plus d'arcs sortants
        if start_node is None:
            for node in all_nodes:
                if self.out_degree[node] > 0 and (start_node is None or self.out_degree[node] > self.out_degree[start_node]):
                    start_node = node
        
        return start_node
    
    def find_longest_greedy_path(self):
        """
        Cherche un chemin eulérien approximativement le plus long possible.
        """
        if self.arc_count == 0:
            return []

        graph_copy = defaultdict(deque)
        for node in self.graph:
            # Trier les arcs pour favoriser les sommets avec plus d’options
            graph_copy[node] = deque(sorted(self.graph[node], key=lambda x: -self.out_degree.get(x[0], 0)))

        start_node = self.find_start_node()
        if start_node is None:
            for node in graph_copy:
                if graph_copy[node]:
                    start_node = node
                    break
        if start_node is None:
            return []

        print(f"[INFO] Début du chemin long au sommet: {start_node}")

        stack = [(start_node, [], set())]  # node, path, used_edges
        longest_path = []

        while stack:
            current_node, path_so_far, used_edges = stack.pop()

            extended = False
            for next_node, kmer in graph_copy[current_node]:
                edge_id = (current_node, next_node, kmer)
                if edge_id not in used_edges:
                    new_path = path_so_far + [kmer]
                    new_used = used_edges | {edge_id}
                    stack.append((next_node, new_path, new_used))
                    extended = True

            if not extended:
                # Aucun arc non visité : comparer la longueur
                if len(path_so_far) > len(longest_path):
                    longest_path = path_so_far

        print(f"[INFO] Longueur du chemin trouvé: {len(longest_path)} arcs")
        return longest_path

        
    
    def reconstruct_sequence(self, kmer_path: list):
        """
        Reconstruit la séquence ADN à partir du chemin de k-mers.
        """
        if not kmer_path:
            return ""
        
        # Commencer avec le premier k-mer
        sequence = kmer_path[0]
        
        # Pour chaque k-mer suivant, ajouter seulement la dernière base
        for next_kmer in kmer_path[1:]:
            sequence += next_kmer[-1]
        
        return sequence

    def find_connected_components(self):
        """
        Trouve les composantes connexes du graphe pour un assemblage par contigs.
        """
        visited = set()
        components = []
        
        def dfs(node, component):
            stack = [node]
            while stack:
                current = stack.pop()
                if current not in visited:
                    visited.add(current)
                    component.append(current)
                    for neighbor, _ in self.graph.get(current, []):
                        if neighbor not in visited:
                            stack.append(neighbor)
                    # Also check incoming neighbors
                    for potential_source in self.graph:
                        for neighbor, _ in self.graph[potential_source]:
                            if neighbor == current and potential_source not in visited:
                                stack.append(potential_source)
        
        for node in self.graph:
            if node not in visited:
                component = []
                dfs(node, component)
                components.append(component)
        
        return components

    def assemble_contigs(self):
        """
        Assemble des contigs à partir des composantes connexes.
        """
        components = self.find_connected_components()
        print(f"[INFO] {len(components)} composantes connexes trouvées")
        
        contigs = []
        for i, component in enumerate(components):
            if len(component) > 1:  # Ignorer les composantes trop petites
                # Créer un sous-graphe pour cette composante
                subgraph = {}
                for node in component:
                    subgraph[node] = [edge for edge in self.graph[node] if edge[0] in component]
                
                # Assembler cette composante
                contig = self._assemble_component(subgraph, component)
                if contig:
                    contigs.append(contig)
                    print(f"[INFO] Contig {i+1}: {len(contig)} bp")
        
        return contigs
    
    def _assemble_component(self, subgraph, component):
        """
        Assemble un contig à partir d'une composante connexe.
        """
        if not subgraph:
            return ""
        
        # Trouver un bon point de départ (nœud avec peu d'entrées)
        start_node = None
        for node in component:
            in_deg = sum(1 for n in subgraph for edge in subgraph[n] if edge[0] == node)
            if in_deg == 0 or (start_node is None and subgraph.get(node)):
                start_node = node
                break
        
        if start_node is None:
            start_node = next(iter(subgraph))
        
        # Parcourir le graphe
        current = start_node
        sequence = current
        visited_edges = set()
        max_steps = len(component) * 2  # Éviter les boucles infinies
        
        for _ in range(max_steps):
            if current not in subgraph or not subgraph[current]:
                break
            
            # Prendre le premier arc non visité
            found_next = False
            for i, (next_node, kmer) in enumerate(subgraph[current]):
                edge_id = (current, next_node, kmer)
                if edge_id not in visited_edges:
                    visited_edges.add(edge_id)
                    sequence += kmer[-1] if len(kmer) == self.k else kmer[-1]
                    current = next_node
                    found_next = True
                    break
            
            if not found_next:
                break
        
        return sequence

###############################################
# FONCTIONS SCALABLES POUR GROS FASTQ
###############################################

def build_kmer_database_to_file(fastq_path: str, k: int, output_file: str):
    """
    Génère un fichier contenant TOUS les k-mers extraits du FASTQ.
    """
    print(f"[INFO] Génération des k-mers {k}-mer → {output_file}")
    
    with open(output_file, "w") as out:
        with open_fastq(fastq_path) as handle:
            for record in SeqIO.parse(handle, "fastq"):
                seq = str(record.seq).upper()
                
                # Filtrer les séquences trop courtes
                if len(seq) >= k:
                    # Génération en streaming
                    for i in range(len(seq) - k + 1):
                        kmer = seq[i:i+k]
                        # Filtrer les k-mers avec des caractères non-ADN
                        if all(base in 'ACGT' for base in kmer):
                            out.write(kmer + "\n")

    print("[OK] Fichier brut de k-mers généré.")

def count_kmers_with_sort(kmer_file: str, output_count_file: str):
    """
    Utilise sort | uniq -c pour compter les k-mers.
    """
    print("[INFO] Comptage des k-mers via sort | uniq -c ...")
    cmd = f"sort {kmer_file} | uniq -c > {output_count_file}"
    subprocess.run(cmd, shell=True, check=True)
    print(f"[OK] Comptage terminé → {output_count_file}")

def build_debruijn_edges_to_file_with_labels(kmer_file: str, output_edges_file: str):
    """
    Construit les arcs du graphe de De Bruijn avec étiquette = kmer.
    """
    print(f"[INFO] Génération des arcs étiquetés → {output_edges_file}")
    
    with open(kmer_file) as km_in, open(output_edges_file, "w") as edges_out:
        for line in km_in:
            kmer = line.strip()
            if len(kmer) < 2:
                continue
            
            prefix = kmer[:-1]
            suffix = kmer[1:]
            edges_out.write(f"{prefix}\t{suffix}\t{kmer}\n")
    
    print("[OK] Arcs étiquetés générés.")

def filter_low_coverage_kmers(kmer_count_file: str, output_filtered: str, min_coverage: int = 2):
    """
    Filtre les k-mers avec une couverture trop faible.
    """
    print(f"[INFO] Filtrage des k-mers avec couverture < {min_coverage}")
    
    with open(kmer_count_file) as f_in, open(output_filtered, "w") as f_out:
        for line in f_in:
            parts = line.strip().split()
            if len(parts) >= 2:
                count = int(parts[0])
                kmer = parts[1]
                if count >= min_coverage and all(base in 'ACGT' for base in kmer):
                    f_out.write(kmer + "\n")
    
    print("[OK] Filtrage terminé")

#######################################
# SCRIPT PRINCIPAL
#######################################
# SCRIPT PRINCIPAL — VERSION CONTIGS COMPLETS
#######################################

path = "reads.fastq.fq"
k = 31

try:
    # Étape A — Génération des k-mers (streaming)
    build_kmer_database_to_file(path, k, "kmers_raw.txt")

    # Étape B — Comptage (sort | uniq -c)
    count_kmers_with_sort("kmers_raw.txt", "kmers_counted.txt")

    # Étape B2 — Filtrage des k-mers à faible couverture
    filter_low_coverage_kmers("kmers_counted.txt", "kmers_filtered.txt", min_coverage=2)

    # Étape C — Construction des arcs du graphe
    build_debruijn_edges_to_file_with_labels("kmers_filtered.txt", "edges_labeled.txt")

    # Étape D — Assemblage
    euler_finder = EulerianPathFinder("edges_labeled.txt", k)
    
    # On récupère tous les contigs possibles
    print("[INFO] Assemblage de tous les contigs...")
    contigs = []

    # Essayer d'abord un chemin eulérien principal
    kmer_path = euler_finder.find_longest_greedy_path()

    if kmer_path:
        sequence_reconstructed = euler_finder.reconstruct_sequence(kmer_path)
        contigs.append(sequence_reconstructed)
        print(f"[INFO] Chemin eulérien principal: {len(sequence_reconstructed)} bp")

    # Puis compléter avec toutes les autres composantes
    contigs_from_components = euler_finder.assemble_contigs()
    contigs.extend(contigs_from_components)

    # Écriture finale dans un fichier FASTA
    with open("assembled_contigs.fasta", "w") as f:
        for i, contig in enumerate(contigs):
            if len(contig) >= k:  # éviter les contigs trop courts
                f.write(f">contig_{i+1}_length_{len(contig)}\n")
                f.write(contig + "\n")

    print(f"[SUCCÈS] {len(contigs)} contigs écrits dans assembled_contigs.fasta")

except Exception as e:
    print(f"[ERREUR] {e}")
    import traceback
    traceback.print_exc()
