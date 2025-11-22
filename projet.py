#!/usr/bin/env python3

from Bio import SeqIO
import gzip
import subprocess
import os
import sys
import argparse
from collections import defaultdict, deque

def open_fastq(path):
    """Ouvre un fichier FASTQ (généralement non compressé dans ce contexte)."""
    return open(path, 'r')

def generate_kmers(sequence: str, k: int):
    """Génère tous les k-mers de taille k à partir d'une séquence donnée."""
    if len(sequence) < k:
        return []
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

class EulerianPathFinder:
    def __init__(self, edges_file: str, k: int):
        """
        Initialise avec le fichier d'arcs étiquetés.
        k: Taille du k-mer utilisé dans le graphe.
        """
        self.k = k
        self.graph = defaultdict(list)  # sommet -> liste des arcs (suffix, kmer, is_fictive)
        self.arc_count = 0
        self.fictive_edges = []
        
        self.in_degree = defaultdict(int)
        self.out_degree = defaultdict(int)

        self.load_graph(edges_file)
        
    def load_graph(self, edges_file: str):
        """Charge le graphe depuis le fichier d'arcs."""
        print("[INFO] Chargement du graphe...")
        try:
            with open(edges_file) as f:
                for line in f:
                    prefix, suffix, kmer = line.strip().split("\t")
                    self.graph[prefix].append((suffix, kmer, False))
                    self.out_degree[prefix] += 1
                    self.in_degree[suffix] += 1
                    self.arc_count += 1
        except FileNotFoundError:
            print(f"[ERREUR] Fichier d'arcs non trouvé: {edges_file}")
            exit(1)
        print(f"[OK] Graphe chargé: {self.arc_count} arcs, {len(self.graph)} sommets")
        
    def balance_graph(self):
        """
        Équilibre le graphe en ajoutant des arêtes fictives si nécessaire.
        """
        print("[INFO] Analyse de l'équilibre du graphe...")
        
        start_nodes = []
        end_nodes = []
        
        all_nodes = set(self.out_degree.keys()) | set(self.in_degree.keys())
        
        for node in all_nodes:
            diff = self.out_degree[node] - self.in_degree[node]
            if diff > 0:
                start_nodes.append((node, diff))
            elif diff < 0:
                end_nodes.append((node, -diff))
        
        fictive_count = 0
        # Connexion simplifiée des nœuds déséquilibrés
        for i in range(min(len(start_nodes), len(end_nodes))):
            start_node, _ = start_nodes[i]
            end_node, _ = end_nodes[i]
            
            # Arête fictive de end_node vers start_node
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
        """Trouve le noeud de départ pour le cycle eulérien."""
        
        all_nodes = set(self.out_degree.keys()) | set(self.in_degree.keys())
        # Noeud naturel (out > in)
        for node in all_nodes:
            if self.out_degree[node] > self.in_degree[node]:
                return node
        
        # Graphe équilibré, prendre le premier nœud avec des arcs sortants
        for node in self.graph:
            if self.graph[node]:
                return node
        
        return next(iter(self.graph)) if self.graph else None
        
    def hierholzer_algorithm(self):
        """Implémentation de l'algorithme de Hierholzer."""
        if self.arc_count == 0:
            return []
        
        graph_copy = defaultdict(deque)
        for node in self.graph:
            graph_copy[node] = deque(self.graph[node])
        
        start_node = self.find_start_node()
        if start_node is None:
            print("[ATTENTION] Aucun sommet de départ trouvé.")
            return []
        
        print(f"[INFO] Début du cycle eulérien au sommet: {start_node}")
        
        stack = [start_node]
        cycle = []
        
        while stack:
            current_node = stack[-1]
            
            if graph_copy[current_node]:
                next_node, kmer_label, is_fictive = graph_copy[current_node].popleft()
                
                stack.append(next_node)
                cycle.append((kmer_label, is_fictive))
            else:
                stack.pop()

        if len(cycle) != self.arc_count:
            print(f"[ATTENTION] Cycle incomplet: {len(cycle)}/{self.arc_count} arcs utilisés.")
        else:
            print(f"[SUCCÈS] Cycle eulérien trouvé: {len(cycle)} arcs")
        
        return cycle
        
    def find_eulerian_path_with_contigs(self):
        """Trouve le chemin eulérien et le fragmente en contigs."""
        print("[INFO] Recherche de chemin eulérien avec gestion des contigs...")
        
        self.balance_graph()
        cycle = self.hierholzer_algorithm()
        contigs = self._split_into_contigs(cycle)
        
        print(f"[SUCCÈS] Génération de {len(contigs)} contigs")
        return contigs
        
    def _split_into_contigs(self, cycle):
        """Fragmente le cycle en contigs aux arêtes fictives."""
        contigs = []
        current_contig = []
        
        for kmer, is_fictive in cycle:
            if is_fictive:
                if current_contig:
                    sequence = self.reconstruct_sequence(current_contig, self.k)
                    if sequence:
                        contigs.append(sequence)
                    current_contig = []
            else:
                current_contig.append(kmer)
        
        # Ajouter le dernier contig
        if current_contig:
            sequence = self.reconstruct_sequence(current_contig, self.k)
            if sequence:
                contigs.append(sequence)
        
        return contigs
        
    def reconstruct_sequence(self, kmer_path: list, k: int):
        """Reconstruit la séquence ADN à partir du chemin de k-mers."""
        if not kmer_path:
            return ""
        
        sequence = kmer_path[0]
        
        for next_kmer in kmer_path[1:]:
            sequence += next_kmer[-1]
        
        return sequence

def build_kmer_database_to_file(fastq_path: str, k: int, output_file: str):
    """Génère un fichier contenant TOUS les k-mers extraits du FASTQ."""
    
    if not os.path.exists(fastq_path):
        raise FileNotFoundError(f"Fichier FASTQ non trouvé à: {fastq_path}")

    try:
        with open(output_file, "w") as out:
            with open_fastq(fastq_path) as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    seq = str(record.seq).upper()
                    for kmer in generate_kmers(seq, k):
                        out.write(kmer + "\n")
        print("[OK] Fichier brut de k-mers généré.")
    except Exception as e:
        print(f"[ERREUR] Échec de la génération des k-mers: {e}")
        exit(1)


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
        
    

def filter_contigs_by_length(input_file: str, output_file: str, min_length: int):
    """
    Filtre les contigs d'un fichier FASTA/FASTQ en fonction de leur longueur.
    Génère un nouveau fichier FASTA contenant uniquement les contigs >= min_length.
    """
    print(f"[INFO] Filtrage des contigs... (longueur minimale: {min_length})")
    
    try:
        # Lecture des enregistrements depuis le fichier FASTA
        records = list(SeqIO.parse(input_file, "fasta"))
        
        # Filtrage des enregistrements
        filtered_records = [
            record for record in records if len(record.seq) >= min_length
        ]
        
        # Écriture des contigs filtrés dans le nouveau fichier
        SeqIO.write(filtered_records, output_file, "fasta")
        
        print(f"[OK] {len(filtered_records)} contigs conservés sur {len(records)} → {output_file}")
        return len(filtered_records)
        
    except FileNotFoundError:
        print(f"[ATTENTION] Fichier d'entrée pour le filtrage non trouvé: {input_file}")
        return 0
    except Exception as e:
        print(f"[ERREUR] Échec du filtrage des contigs: {e}")
        return 0


if __name__ == "__main__":
    
    #  Définition des valeurs par défaut 
    DEFAULT_PATH = 'reads.fastq.fq'
    DEFAULT_K = 31
    DEFAULT_MIN_LEN = 1000 


    """G = pgv.AGraph(directed=True, strict=False)  # strict=False pour autoriser multi-arêtes

    with open(edges_file) as f:
        for line in f:
            prefix, suffix, kmer = line.strip().split("\t")
            # Ajouter l'arête avec étiquette
            G.add_edge(prefix, suffix, label=kmer)

    # Mise en page automatique avec Graphviz (dot, neato, fdp, etc.)
    G.layout(prog='dot')  # 'dot' pour DAG, 'neato' ou 'fdp' pour plus libre

    # Sauvegarder en PNG ou PDF
    #G.draw("debruijn_graph.png")
    G.draw("debruijn_graph.pdf")"""
    

    # Définiton du format de la ligne de commande par l'utilisateur  
    parser = argparse.ArgumentParser(
        description="Assembleur De Bruijn à Chemin Eulérien.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=f"Exemple: python3 projet.py --filename reads.fastq.fq -k 40 --filter 100\n"
               f"Par défaut: --filename '{DEFAULT_PATH}', -k {DEFAULT_K}, --filter {DEFAULT_MIN_LEN}"
    )
    
    # ARGUMENT MODIFIÉ : -f remplacé par --filename
    parser.add_argument(
        '--filename',
        type=str,
        default=DEFAULT_PATH,
        dest='path', # Utilisation de 'dest' pour maintenir le nom de variable 'path'
        help="Chemin vers le fichier FASTQ d'entrée (reads)"
    )
    
    parser.add_argument(
        '-k', '--kmer',
        type=int,
        default=DEFAULT_K,
        help="Taille du k-mer (k) pour l'assemblage"
    )
    
    # ARGUMENT MODIFIÉ : -l remplacé par --filter
    parser.add_argument(
        '--filter',
        type=int,
        default=DEFAULT_MIN_LEN,
        dest='min_length', # Utilisation de 'dest' pour maintenir le nom de variable 'min_length'
        help="Longueur minimale des contigs à conserver après l'assemblage"
    )

    args = parser.parse_args()
    

    # Récupération des arguments par leurs noms de 'dest'
    path = args.path
    k = args.kmer
    min_length = args.min_length 
    

    print("--- Assembleur De Bruijn à Chemin Eulérien ---")
    print(f"Fichier d'entrée: {path}")
    print(f"Taille du k-mer (k): {k}")
    print(f"Longueur minimale pour le filtrage (--filter): {min_length}")

    # Définition des noms de fichiers intermédiaires et de sortie
    kmer_raw_file = "kmers_raw.txt"
    kmer_counted_file = "kmers_counted.txt"
    edges_labeled_file = "edges_labeled.txt"
    output_contig_file = "assembled_contigs.fasta"
    output_contig_filtered_file = "assembled_contigs_filtered.fasta" # Nouveau fichier de sortie filtré

    # --- ÉTAPES D'ASSEMBLAGE ---
    
        # Étape A — Génération des k-mers
    build_kmer_database_to_file(path, k, kmer_raw_file)

        # Étape B — Comptage (facultatif)
    #count_kmers_with_sort(kmer_raw_file, kmer_counted_file)

        # Étape C — Construction des arcs du graphe de De Bruijn
    build_debruijn_edges_to_file_with_labels(kmer_raw_file, edges_labeled_file)
        
        # Étape D — Recherche du chemin eulérien et reconstruction des contigs
    euler_finder = EulerianPathFinder(edges_labeled_file, k) 
    contigs = euler_finder.find_eulerian_path_with_contigs()

        # Étape E — Sauvegarde des contigs (bruts)
    print(f"[INFO] Sauvegarde des contigs bruts → {output_contig_file}")
    with open(output_contig_file, "w") as f:
        for i, contig in enumerate(contigs):
            f.write(f">contig_{i} Length={len(contig)}\n{contig}\n")
    print(f"[OK] {len(contigs)} contigs enregistrés.")

        # Étape F — Filtrage des contigs
    filter_contigs_by_length(output_contig_file, output_contig_filtered_file, min_length)
        
        # Nettoyage des fichiers intermédiaires
    print("\n[INFO] Nettoyage des fichiers intermédiaires...")
    files_to_clean = [kmer_raw_file, kmer_counted_file]
    for file in files_to_clean:
        try:
            if os.path.exists(file):
                os.remove(file)
        except OSError as e:
            print(f"[ATTENTION] Impossible de supprimer {file}: {e}")

    print(f"\n[TERMINÉ] L'assemblage et le filtrage sont terminés.")
    print(f"Contigs bruts: '{output_contig_file}'")
    print(f"Contigs filtrés (>= {min_length} bp): '{output_contig_filtered_file}'")
    print("-" * 40)