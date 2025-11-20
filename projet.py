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
    def __init__(self, edges_file: str):
        """
        Initialise avec le fichier d'arcs étiquetés.
        """
        self.graph = defaultdict(list)  # sommet -> liste des arcs (suffix, kmer, is_fictive)
        self.arc_count = 0
        self.fictive_edges = []  # Pour suivre les arêtes fictives ajoutées
        
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
                self.graph[prefix].append((suffix, kmer, False))  # False = arête réelle
                self.out_degree[prefix] += 1
                self.in_degree[suffix] += 1
                self.arc_count += 1
        print(f"[OK] Graphe chargé: {self.arc_count} arcs, {len(self.graph)} sommets")
    
    def balance_graph(self):
        """
        Équilibre le graphe en ajoutant des arêtes fictives si nécessaire.
        Retourne le nombre d'arêtes fictives ajoutées.
        """
        print("[INFO] Analyse de l'équilibre du graphe...")
        
        start_nodes = []  # nœuds avec out_degree > in_degree
        end_nodes = []    # nœuds avec in_degree > out_degree
        
        all_nodes = set(self.out_degree.keys()) | set(self.in_degree.keys())
        
        for node in all_nodes:
            diff = self.out_degree[node] - self.in_degree[node]
            if diff > 0:
                start_nodes.append((node, diff))
            elif diff < 0:
                end_nodes.append((node, -diff))
        
        print(f"[INFO] Déséquilibre: {len(start_nodes)} start nodes, {len(end_nodes)} end nodes")
        
        # Ajouter des arêtes fictives pour équilibrer
        fictive_count = 0
        for (start_node, start_count), (end_node, end_count) in zip(start_nodes, end_nodes):
            # Ajouter une arête fictive de end_node vers start_node
            fictive_kmer = f"FICTIVE_{fictive_count}"
            self.graph[end_node].append((start_node, fictive_kmer, True))
            self.out_degree[end_node] += 1
            self.in_degree[start_node] += 1
            self.arc_count += 1
            fictive_count += 1
            
            self.fictive_edges.append((end_node, start_node, fictive_kmer))
        
        print(f"[INFO] {fictive_count} arêtes fictives ajoutées")
        return fictive_count
    
    def find_start_node(self):
        """
        Trouve le nœud de départ pour le cycle eulérien.
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
                next_node, kmer_label, is_fictive = graph_copy[current_node].popleft()
                
                # Ajouter le nœud suivant à la pile
                stack.append(next_node)
                
            else:
                # Plus d'arcs sortants, ajouter au cycle et backtrack
                if len(stack) > 1:
                    from_node = stack[-2]
                    to_node = stack[-1]
                    
                    # Trouver le k-mer correspondant
                    for i, (n, kmer, fictive) in enumerate(self.graph[from_node]):
                        if n == to_node:
                            cycle.append((kmer, fictive))  # Stocker avec l'info fictive
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
    
    def find_eulerian_path_with_contigs(self):
        """
        Trouve un chemin eulérien même pour les graphes non équilibrés.
        Retourne une liste de contigs (fragments).
        """
        print("[INFO] Recherche de chemin eulérien avec gestion des contigs...")
        
        # Équilibrer le graphe si nécessaire
        fictive_count = self.balance_graph()
        
        if fictive_count == 0:
            print("[INFO] Graphe déjà équilibré - recherche de cycle simple")
            cycle = self.hierholzer_algorithm()
            return [self._extract_sequence_from_path(cycle)]
        
        # Trouver le cycle eulérien dans le graphe équilibré
        cycle = self.hierholzer_algorithm()
        
        # Fragmenter en contigs aux arêtes fictives
        contigs = self._split_into_contigs(cycle)
        
        print(f"[SUCCÈS] Génération de {len(contigs)} contigs")
        return contigs
    
    def _split_into_contigs(self, cycle):
        """
        Fragmente le cycle en contigs aux points des arêtes fictives.
        """
        contigs = []
        current_contig = []
        
        for kmer, is_fictive in cycle:
            if is_fictive:
                # Arête fictive trouvée - terminer le contig courant
                if current_contig:
                    sequence = self.reconstruct_sequence(current_contig, 3)  # k=3 par défaut
                    if sequence:  # Ne pas ajouter de contigs vides
                        contigs.append(sequence)
                    current_contig = []
            else:
                # Arête réelle - ajouter au contig courant
                current_contig.append(kmer)
        
        # Ajouter le dernier contig
        if current_contig:
            sequence = self.reconstruct_sequence(current_contig, 3)
            if sequence:
                contigs.append(sequence)
        
        return contigs
    
    def _extract_sequence_from_path(self, path):
        """Extrait la séquence d'un chemin de k-mers."""
        kmers_only = [kmer for kmer, fictive in path if not fictive]
        return self.reconstruct_sequence(kmers_only, 3)
    
    def reconstruct_sequence(self, kmer_path: list, k: int):
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

# Les fonctions restantes sont identiques à votre version originale
def build_kmer_database_to_file(fastq_path: str, k: int, output_file: str):
    """Génère un fichier contenant TOUS les k-mers extraits du FASTQ."""
    print(f"[INFO] Génération des k-mers {k}-mer → {output_file}")
    
    with open(output_file, "w") as out:
        with open_fastq(fastq_path) as handle:
            for record in SeqIO.parse(handle, "fastq"):
                seq = str(record.seq).upper()
                for kmer in generate_kmers(seq, k):
                    out.write(kmer + "\n")
    print("[OK] Fichier brut de k-mers généré.")

def count_kmers_with_sort(kmer_file: str, output_count_file: str):
    """Compte les k-mers avec sort | uniq -c"""
    print("[INFO] Comptage des k-mers via sort | uniq -c ...")
    cmd = f"sort {kmer_file} | uniq -c > {output_count_file}"
    subprocess.run(cmd, shell=True, check=True)
    print(f"[OK] Comptage terminé → {output_count_file}")

def build_debruijn_edges_to_file_with_labels(kmer_file: str, output_edges_file: str):
    """Construit les arcs du graphe de De Bruijn avec étiquette = kmer."""
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

# Exemple d'utilisation
if __name__ == "__main__":
    path = "reads.fastq.fq"  # Chemin vers le fichier FASTQ
    k = 31

    # Étape A — Génération des k-mers
    build_kmer_database_to_file(path, k, "kmers_raw.txt")

    # Étape B — Comptage
    count_kmers_with_sort("kmers_raw.txt", "kmers_counted.txt")

    # Étape C — Construction des arcs du graphe
    build_debruijn_edges_to_file_with_labels("kmers_raw.txt", "edges_labeled.txt")

    # Étape D — Recherche du chemin eulérien avec gestion des contigs
    euler_finder = EulerianPathFinder("edges_labeled.txt")
    
    # Cette méthode gère automatiquement les graphes non eulériens
    contigs = euler_finder.find_eulerian_path_with_contigs()

    # Sauvegarde des contigs
    with open("assembled_contigs.txt", "w") as f:
        for i, contig in enumerate(contigs):
            f.write(f">contig_{i}\n{contig}\n")
    
    print(f"[TERMINÉ] {len(contigs)} contigs générés dans 'assembled_contigs.txt'")