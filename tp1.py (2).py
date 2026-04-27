def adn_en_proteine(sequence_adn):
    """
    Transforme une séquence d'ADN en séquence de protéine.
    
    Args:
        sequence_adn (str): Séquence d'ADN (A, T, G, C)
    
    Returns:
        str: Séquence d'acides aminés ou message d'erreur
    """
    
    # Dictionnaire du code génétique (codon -> acide aminé)
    code_genetique = {
        # Alanine
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        # Arginine
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
        # Asparagine
        'AAT': 'N', 'AAC': 'N',
        # Acide aspartique
        'GAT': 'D', 'GAC': 'D',
        # Cystéine
        'TGT': 'C', 'TGC': 'C',
        # Glutamine
        'CAA': 'Q', 'CAG': 'Q',
        # Acide glutamique
        'GAA': 'E', 'GAG': 'E',
        # Glycine
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
        # Histidine
        'CAT': 'H', 'CAC': 'H',
        # Isoleucine
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
        # Leucine
        'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        # Lysine
        'AAA': 'K', 'AAG': 'K',
        # Méthionine (démarrage)
        'ATG': 'M',
        # Phénylalanine
        'TTT': 'F', 'TTC': 'F',
        # Proline
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        # Sérine
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
        # Thréonine
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        # Tryptophane
        'TGG': 'W',
        # Tyrosine
        'TAT': 'Y', 'TAC': 'Y',
        # Valine
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        # Codons stop
        'TAA': '*', 'TAG': '*', 'TGA': '*'
    }
    
    # Validation de la séquence d'entrée
    sequence_adn = sequence_adn.upper().strip()
    nucleotides_valides = set('ATGC')
    
    if not all(nucl in nucleotides_valides for nucl in sequence_adn):
        return "Erreur: La séquence contient des caractères non valides (A, T, G, C uniquement)"
    
    # Recherche du codon de démarrage (ATG)
    start_index = sequence_adn.find('ATG')
    if start_index == -1:
        return "Erreur: Aucun codon de démarrage (ATG) trouvé"
    
    # Extraction de la séquence à partir du codon de démarrage
    sequence_codante = sequence_adn[start_index:]
    
    # Traduction en protéine
    proteine = []
    
    for i in range(0, len(sequence_codante) - 2, 3):
        codon = sequence_codante[i:i+3]
        
        if len(codon) != 3:
            break
            
        if codon in code_genetique:
            acide_amine = code_genetique[codon]
            
            # Arrêt si codon stop
            if acide_amine == '*':
                break
                
            proteine.append(acide_amine)
        else:
            return f"Erreur: Codon non reconnu à la position {i} : {codon}"
    
    return ''.join(proteine)


def lire_fichier_adn(nom_fichier):
    """Lit une séquence ADN depuis un fichier"""
    try:
        with open(nom_fichier, 'r') as fichier:
            lignes = fichier.readlines()
            # Ignorer les lignes de commentaire et en-têtes FASTA
            sequence = ''
            for ligne in lignes:
                if not ligne.startswith('>') and not ligne.startswith('#'):
                    sequence += ligne.strip()
            return sequence
    except FileNotFoundError:
        return None


def sauvegarder_proteine(nom_fichier, proteine):
    """Sauvegarde la séquence protéique dans un fichier"""
    with open(nom_fichier, 'w') as fichier:
        fichier.write(proteine)
        # Écrire la séquence avec des sauts de ligne tous les 60 caractères
        for i in range(0, len(proteine), 60):
            fichier.write(proteine[i:i+60] + '\n')


# Exemple d'utilisation
if __name__ == "__main__":
    # Exemple avec une séquence simple
    sequence_test = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    
    print("=" * 50)
    print("TRADUCTION ADN → PROTÉINE")
    print("=" * 50)
    
    print(f"\nSéquence ADN: {sequence_test}")
    
    proteine = adn_en_proteine(sequence_test)
    print(f"Séquence protéique: {proteine}")
    
    # Exemple avec plusieurs cadres de lecture
    print("\n" + "=" * 50)
    print("CADRES DE LECTURE POSSIBLES")
    print("=" * 50)
    
    sequence_longue = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGGGCCGCTGAA"
    
    for frame in range(3):
        seq_cadre = sequence_longue[frame:]
        proteine_frame = adn_en_proteine(seq_cadre)
        print(f"Cadre {frame + 1}: {proteine_frame}")
    
    # Exemple avec gestion d'erreur
    print("\n" + "=" * 50)
    print("GESTION DES ERREURS")
    print("=" * 50)
    
    sequence_erreur = "ATGGCCXTTGTAATG"  # Contient 'X' invalide
    resultat = adn_en_proteine(sequence_erreur)
    print(f"Test avec séquence invalide: {resultat}")
    
    # Fonctionnalités supplémentaires
    print("\n" + "=" * 50)
    print("INFORMATIONS SUPPLÉMENTAIRES")
    print("=" * 50)
    
    # Statistiques de la protéine
    proteine_stat = adn_en_proteine(sequence_test)
    if not proteine_stat.startswith("Erreur"):
        acides_amines = set(proteine_stat)
        print(f"Longueur de la protéine: {len(proteine_stat)} acides aminés")
        print(f"Acides aminés présents: {', '.join(sorted(acides_amines))}")
        
        # Composition en pourcentage
        print("\nComposition détaillée:")
        for aa in sorted(acides_amines):
            count = proteine_stat.count(aa)
            pourcentage = (count / len(proteine_stat)) * 100
            print(f"  {aa}: {count} ({pourcentage:.1f}%)")