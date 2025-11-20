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

def calculate_overlap(seq1, seq2, min_overlap_ratio=0.8):
    """Calcule le chevauchement maximum entre deux séquences."""
    max_overlap = 0
    min_len = min(len(seq1), len(seq2))
    required_overlap = int(min_len * min_overlap_ratio)
    
    for overlap_size in range(min_len, required_overlap - 1, -1):
        if seq1.endswith(seq2[:overlap_size]):
            max_overlap = max(max_overlap, overlap_size)
        if seq2.endswith(seq1[:overlap_size]):
            max_overlap = max(max_overlap, overlap_size)
    
    return max_overlap

def deduplicate_contigs(contigs, min_overlap_ratio=0.9):
    """
    Supprime les contigs qui sont des sous-séquences ou ont de forts chevauchements.
    """
    if not contigs:
        return []
    
    # Trier par longueur décroissante
    contigs.sort(key=len, reverse=True)
    unique_contigs = []
    
    for contig in contigs:
        is_duplicate = False
        
        # Vérifier si ce contig est contenu dans un contig plus long
        for unique in unique_contigs:
            if contig in unique:
                is_duplicate = True
                break
        
        # Vérifier les chevauchements importants
        if not is_duplicate:
            for unique in unique_contigs:
                overlap = calculate_overlap(contig, unique, min_overlap_ratio)
                min_len = min(len(contig), len(unique))
                if overlap >= min_len * min_overlap_ratio:
                    is_duplicate = True
                    break
        
        if not is_duplicate:
            unique_contigs.append(contig)
    
    return unique_contigs

class RobustDeBruijnAssembler:
    def __init__(self, edges_file: str, k: int):
        self.graph = defaultdict(list)
        self.arc_count = 0
        self.k = k
        self.in_degree = defaultdict(int)
        self.out_degree = defaultdict(int)
        self.edge_usage = defaultdict(int)  # Suivi de l'utilisation des arêtes
        self.load_graph(edges_file)
        
    def load_graph(self, edges_file: str):
        """Charge le graphe depuis le fichier d'arcs."""
        print("[INFO] Chargement du graphe...")
        with open(edges_file) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    prefix, suffix, kmer = parts[:3]
                    self.graph[prefix].append((suffix, kmer))
                    self.out_degree[prefix] += 1
                    self.in_degree[suffix] += 1
                    self.arc_count += 1
                    
                    # Initialiser le compteur d'usage pour cette arête
                    edge_id = (prefix, suffix, kmer)
                    self.edge_usage[edge_id] = 0
                    
        print(f"[OK] Graphe chargé: {self.arc_count} arcs, {len(self.graph)} sommets")
    
    def is_eulerian(self):
        """Vérifie si le graphe a un chemin/cycle eulérien."""
        start_nodes = 0
        end_nodes = 0
        
        all_nodes = set(self.out_degree.keys()) | set(self.in_degree.keys())
        
        for node in all_nodes:
            diff = self.out_degree[node] - self.in_degree[node]
            if diff == 1:
                start_nodes += 1
            elif diff == -1:
                end_nodes += 1
            elif diff != 0:
                return False
        
        return (start_nodes == 0 and end_nodes == 0) or (start_nodes == 1 and end_nodes == 1)
    
    def find_eulerian_path_or_cycle(self):
        """
        Implémentation robuste de Hierholzer avec tracking strict des arêtes.
        """
        if self.arc_count == 0:
            return []
        
        # Réinitialiser l'usage des arêtes
        for edge_id in self.edge_usage:
            self.edge_usage[edge_id] = 0
        
        # Faire une copie du graphe avec comptage des multi-arêtes
        graph_copy = defaultdict(deque)
        edge_counter = defaultdict(int)
        
        for node in self.graph:
            for suffix, kmer in self.graph[node]:
                edge_id = (node, suffix, kmer)
                edge_counter[edge_id] += 1
                graph_copy[node].append((suffix, kmer))
        
        # Trier les arcs pour favoriser ceux avec moins d'usage
        for node in graph_copy:
            graph_copy[node] = deque(sorted(graph_copy[node], 
                                          key=lambda x: self.edge_usage.get((node, x[0], x[1]), 0)))
        
        # Trouver le point de départ
        start_node = self._find_optimal_start_node()
        if start_node is None:
            return []
        
        print(f"[INFO] Début du parcours eulérien au sommet: {start_node}")
        
        stack = [start_node]
        path = []
        local_used_edges = defaultdict(int)  # Usage local pour ce parcours
        
        while stack:
            current_node = stack[-1]
            
            if graph_copy[current_node]:
                next_node, kmer_label = graph_copy[current_node].popleft()
                edge_id = (current_node, next_node, kmer_label)
                
                # Vérifier si on peut utiliser cette arête
                if local_used_edges[edge_id] < edge_counter[edge_id]:
                    local_used_edges[edge_id] += 1
                    self.edge_usage[edge_id] += 1  # Usage global
                    stack.append(next_node)
            else:
                if len(stack) > 1:
                    from_node = stack[-2]
                    to_node = stack[-1]
                    
                    # Trouver un k-mer correspondant non surutilisé
                    for n, kmer in self.graph[from_node]:
                        if n == to_node:
                            edge_id = (from_node, to_node, kmer)
                            if local_used_edges[edge_id] > 0:
                                path.append(kmer)
                                break
                stack.pop()
        
        path.reverse()
        
        # Vérifier le résultat
        used_edges = sum(local_used_edges.values())
        if used_edges == self.arc_count:
            print(f"[SUCCÈS] Chemin eulérien complet trouvé: {used_edges} arcs")
        else:
            print(f"[INFO] Chemin partiel: {used_edges}/{self.arc_count} arcs")
        
        return path
    
    def _find_optimal_start_node(self):
        """Trouve le meilleur nœud de départ pour le parcours eulérien."""
        all_nodes = set(self.out_degree.keys()) | set(self.in_degree.keys())
        start_candidates = []
        
        # Priorité 1: Nœuds avec excédent sortant (pour chemins eulériens)
        for node in all_nodes:
            if self.out_degree[node] - self.in_degree[node] == 1:
                return node
        
        # Priorité 2: Nœuds avec arcs sortants mais sans entrées (début naturel)
        for node in all_nodes:
            if self.out_degree[node] > 0 and self.in_degree[node] == 0:
                return node
        
        # Priorité 3: Nœuds avec le meilleur score combiné
        for node in all_nodes:
            if self.out_degree[node] > 0:
                # Score basé sur le déséquilibre et la "fraîcheur" des arêtes
                imbalance = self.out_degree[node] - self.in_degree[node]
                freshness = sum(1 for edge in self.graph[node] 
                              if self.edge_usage.get((node, edge[0], edge[1]), 0) == 0)
                score = imbalance + (freshness / max(1, len(self.graph[node])))
                start_candidates.append((node, score))
        
        if start_candidates:
            start_candidates.sort(key=lambda x: -x[1])
            return start_candidates[0][0]
        
        # Fallback: premier nœud avec arcs sortants
        for node in self.graph:
            if self.graph[node]:
                return node
        
        return None
    
    def find_connected_components(self):
        """Trouve les composantes connexes du graphe."""
        visited = set()
        components = []
        
        def dfs(node, component):
            stack = [node]
            while stack:
                current = stack.pop()
                if current not in visited:
                    visited.add(current)
                    component.append(current)
                    # Arcs sortants
                    for neighbor, _ in self.graph.get(current, []):
                        if neighbor not in visited:
                            stack.append(neighbor)
                    # Arcs entrants
                    for source in self.graph:
                        for neighbor, _ in self.graph[source]:
                            if neighbor == current and source not in visited:
                                stack.append(source)
        
        for node in self.graph:
            if node not in visited:
                component = []
                dfs(node, component)
                components.append(component)
        
        return components
    
    def assemble_contigs_from_components(self):
        """Assemble des contigs à partir des composantes connexes."""
        components = self.find_connected_components()
        print(f"[INFO] {len(components)} composantes connexes trouvées")
        
        contigs = []
        for i, component in enumerate(components):
            if len(component) >= 2:  # Ignorer les composantes triviales
                contig = self._assemble_single_component(component)
                if contig and len(contig) >= self.k:
                    contigs.append(contig)
                    print(f"[INFO] Contig {i+1}: {len(contig)} bp")
        
        return contigs
    
    def _assemble_single_component(self, component):
        """Assemble un contig à partir d'une composante connexe."""
        if not component:
            return ""
        
        # Créer un sous-graphe pour la composante
        subgraph = {}
        
        for node in component:
            # Filtrer les arêtes déjà utilisées et celles dans la composante
            available_edges = []
            for suffix, kmer in self.graph[node]:
                if suffix in component:
                    edge_id = (node, suffix, kmer)
                    # Préférer les arêtes peu utilisées
                    usage = self.edge_usage.get(edge_id, 0)
                    available_edges.append((suffix, kmer, usage))
            
            # Trier par usage croissant (préférer les arêtes fraîches)
            available_edges.sort(key=lambda x: x[2])
            subgraph[node] = [(suffix, kmer) for suffix, kmer, _ in available_edges]
        
        # Trouver le meilleur point de départ
        start_node = None
        best_score = -float('inf')
        
        for node in component:
            if subgraph.get(node):
                in_deg = sum(1 for n in subgraph for edge in subgraph[n] if edge[0] == node)
                out_deg = len(subgraph[node])
                # Score favorisant les nœuds avec peu d'entrées et beaucoup de sorties
                score = out_deg - (in_deg * 0.5)
                if score > best_score:
                    best_score = score
                    start_node = node
        
        if start_node is None:
            for node in component:
                if subgraph.get(node):
                    start_node = node
                    break
        
        if start_node is None:
            return ""
        
        # Parcourir la composante
        current = start_node
        sequence = current
        visited_edges = set()
        max_steps = len(component) * 3
        
        for step in range(max_steps):
            if current not in subgraph or not subgraph[current]:
                break
            
            # Choisir la meilleure arête suivante
            best_next = None
            best_kmer = None
            min_usage = float('inf')
            
            for next_node, kmer in subgraph[current]:
                edge_id = (current, next_node, kmer)
                if edge_id not in visited_edges:
                    usage = self.edge_usage.get(edge_id, 0)
                    if usage < min_usage:
                        min_usage = usage
                        best_next = next_node
                        best_kmer = kmer
            
            if best_next is None:
                break
                
            edge_id = (current, best_next, best_kmer)
            visited_edges.add(edge_id)
            self.edge_usage[edge_id] += 1  # Marquer comme utilisé
            sequence += best_kmer[-1]
            current = best_next
        
        return sequence
    
    def smart_assemble(self):
        """
        Stratégie d'assemblage intelligente avec dé-duplication.
        """
        print("[INFO] Début de l'assemblage intelligent...")
        
        # Essai 1: Chemin eulérien avec tracking amélioré
        euler_path = self.find_eulerian_path_or_cycle()
        
        all_contigs = []
        
        if euler_path:
            main_sequence = self.reconstruct_sequence(euler_path)
            coverage = len(euler_path) / self.arc_count
            
            if coverage > 0.95:  # 95% des arcs utilisés
                print(f"[SUCCÈS] Assemblage eulérien complet: {len(main_sequence)} bp")
                all_contigs = [main_sequence]
            else:
                print(f"[INFO] Assemblage eulérien partiel ({coverage:.1%}), recherche de contigs supplémentaires...")
                contigs = self.assemble_contigs_from_components()
                all_contigs = [main_sequence] + [c for c in contigs if c != main_sequence]
        else:
            print("[INFO] Aucun chemin eulérien trouvé, assemblage par contigs...")
            all_contigs = self.assemble_contigs_from_components()
        
        # Dé-duplication agressive
        print("[INFO] Dé-duplication des contigs...")
        initial_count = len(all_contigs)
        unique_contigs = deduplicate_contigs(all_contigs, min_overlap_ratio=0.85)
        
        print(f"[INFO] Dé-duplication: {initial_count} → {len(unique_contigs)} contigs")
        return unique_contigs
    
    def reconstruct_sequence(self, kmer_path: list):
        """Reconstruit la séquence ADN à partir du chemin de k-mers."""
        if not kmer_path:
            return ""
        
        sequence = kmer_path[0]
        for next_kmer in kmer_path[1:]:
            sequence += next_kmer[-1]
        
        return sequence

###############################################
# FONCTIONS DE PRÉ-TRAITEMENT AMÉLIORÉES
###############################################

def build_kmer_database_to_file(fastq_path: str, k: int, output_file: str):
    """Génère un fichier contenant TOUS les k-mers extraits du FASTQ."""
    print(f"[INFO] Génération des k-mers {k}-mer → {output_file}")
    
    total_sequences = 0
    total_kmers = 0
    
    with open(output_file, "w") as out:
        with open_fastq(fastq_path) as handle:
            for record in SeqIO.parse(handle, "fastq"):
                total_sequences += 1
                seq = str(record.seq).upper()
                
                if len(seq) >= k:
                    for i in range(len(seq) - k + 1):
                        kmer = seq[i:i+k]
                        # Filtrage strict des k-mers
                        if all(base in 'ACGT' for base in kmer):
                            out.write(kmer + "\n")
                            total_kmers += 1

    print(f"[OK] Fichier brut de k-mers généré: {total_sequences} séquences, {total_kmers} k-mers")

def count_kmers_with_sort(kmer_file: str, output_count_file: str):
    """Utilise sort | uniq -c pour compter les k-mers."""
    print("[INFO] Comptage des k-mers via sort | uniq -c ...")
    cmd = f"sort {kmer_file} | uniq -c > {output_count_file}"
    result = subprocess.run(cmd, shell=True, check=True)
    print(f"[OK] Comptage terminé → {output_count_file}")

def build_debruijn_edges_to_file_with_labels(kmer_file: str, output_edges_file: str):
    """Construit les arcs du graphe de De Bruijn avec étiquette = kmer."""
    print(f"[INFO] Génération des arcs étiquetés → {output_edges_file}")
    
    kmer_count = 0
    with open(kmer_file) as km_in, open(output_edges_file, "w") as edges_out:
        for line in km_in:
            kmer = line.strip()
            if len(kmer) >= 2:
                prefix = kmer[:-1]
                suffix = kmer[1:]
                edges_out.write(f"{prefix}\t{suffix}\t{kmer}\n")
                kmer_count += 1
    
    print(f"[OK] Arcs étiquetés générés: {kmer_count} k-mers")

def filter_low_coverage_kmers(kmer_count_file: str, output_filtered: str, min_coverage: int = 3):
    """
    Filtre les k-mers avec une couverture trop faible.
    Version améliorée avec statistiques.
    """
    print(f"[INFO] Filtrage des k-mers avec couverture < {min_coverage}")
    
    total_kmers = 0
    filtered_kmers = 0
    
    with open(kmer_count_file) as f_in, open(output_filtered, "w") as f_out:
        for line in f_in:
            total_kmers += 1
            parts = line.strip().split()
            if len(parts) >= 2:
                count = int(parts[0])
                kmer = parts[1]
                if count >= min_coverage and all(base in 'ACGT' for base in kmer):
                    f_out.write(kmer + "\n")
                    filtered_kmers += 1
    
    print(f"[OK] Filtrage terminé: {filtered_kmers}/{total_kmers} k-mers conservés ({filtered_kmers/total_kmers*100:.1f}%)")

#######################################
# SCRIPT PRINCIPAL
#######################################

def main():
    # Configuration - Ajuster selon les données
    path = "reads.fastq.fq"  # Chemin vers votre fichier FASTQ
    k = 31  # Taille des k-mers (31 pour standard, 51 pour meilleure précision)
    min_coverage = 3  # Filtrage des erreurs (2 pour sensibilité, 3 pour précision)
    
    try:
        print("=" * 60)
        print("ASSEMBLEUR DE BRUJIN ROBUSTE")
        print("=" * 60)
        
        # Étape A — Génération des k-mers
        print("\n[ÉTAPE 1] Génération des k-mers...")
        build_kmer_database_to_file(path, k, "kmers_raw.txt")

        # Étape B — Comptage
        print("\n[ÉTAPE 2] Comptage des k-mers...")
        count_kmers_with_sort("kmers_raw.txt", "kmers_counted.txt")

        # Étape C — Filtrage
        print(f"\n[ÉTAPE 3] Filtrage (couverture ≥ {min_coverage})...")
        filter_low_coverage_kmers("kmers_counted.txt", "kmers_filtered.txt", min_coverage)

        # Étape D — Construction du graphe
        print("\n[ÉTAPE 4] Construction du graphe de De Bruijn...")
        build_debruijn_edges_to_file_with_labels("kmers_filtered.txt", "edges_labeled.txt")

        # Étape E — Assemblage intelligent
        print("\n[ÉTAPE 5] Assemblage...")
        assembler = RobustDeBruijnAssembler("edges_labeled.txt", k)
        contigs = assembler.smart_assemble()

        # Étape F — Sauvegarde des résultats
        print("\n[ÉTAPE 6] Sauvegarde des résultats...")
        with open("final_assembly.fasta", "w") as f:
            total_length = 0
            for i, contig in enumerate(contigs):
                if len(contig) >= k:  # Filtrer les contigs trop courts
                    f.write(f">contig_{i+1}_length_{len(contig)}\n")
                    f.write(contig + "\n")
                    total_length += len(contig)
                    print(f"  Contig {i+1}: {len(contig)} bp")
            
            print(f"\n[SUCCÈS] Assemblage terminé!")
            print(f"  • {len(contigs)} contigs")
            print(f"  • {total_length} bp total")
            print(f"  • Taille moyenne: {total_length/max(1, len(contigs)):.0f} bp")
            print(f"  • Fichier de sortie: final_assembly.fasta")

    except Exception as e:
        print(f"\n[ERREUR] {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()