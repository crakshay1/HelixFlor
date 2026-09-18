# Inclure des librairies GBOT
import os
import sys
import argparse
import re

base_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(base_dir, '..', 'flagdbplante.cfg')
os.environ['DATABASE_FILE'] = config_path
os.environ['DEBUG'] = "False"

# Ajout du chemin de recherche des librairies (Elle se trouve dans un repertoire au dessus)
p = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, p)

from FLAGdbServer import setting                                    # module d'initialisation des variables
from FLAGdbServer.tools import logger                               # Module d'initialisation des logs
from FLAGdbServer.database import FLAGdb                            # Module de connection a la base
from FLAGdbServer.apis.species import speciesDAO, SpeciesNotFound   # Module pour interroger des especes
from FLAGdbServer.apis.feature import featureDAO                    # Module pour interroger des features
from FLAGdbServer.apis.sequence import sequenceDAO, rev_compl       # Module pour interroger des séquences
from FLAGdbServer.apis.qualifier import qualifierDAO                # Module pour interroger des qualifiers
from ToolsServer.extractSequence import calcProtein                 # Module pour prédire les protéines

# Méthodes de l'outil

def normalize_chr(chr_name):
    # Scaffold : on ne touche pas
    if re.match(r"^[Ss]caffold", chr_name):
        return chr_name

    # Mitochondrie
    if chr_name.endswith(("M", "m")):
        return "chrM"

    # Chloroplaste
    if chr_name.endswith(("C", "c")):
        return "chrC"

    # Chromosomes numériques
    number = re.search(r"\d+", chr_name)
    if number:
        return "chr" + str(int(number.group()))

    return chr_name


def are_overlapping(start1: int, stop1: int, start2: int, stop2: int, seuil: int) -> bool:
    """
        Prend en paramètre les starts et stops de deux éléments.
        Renvoie True si les deux éléments se chevauchent, sinon False.
    """
    return min(stop1, stop2) - max(start1, start2) > seuil # Si c'est négatif -> Pas de chevauchement


def find_missingAnnotationV1(id : int, feature : str) -> list :
    # Cette fonction servait de V1 et posait problème en termes de compléxité... D'où sa version améliorée comparing_annotation() !
    """
        Prend en paramètre l'id de l'espèce et le feat d'intérêt.
        Renvoie un dict avec pour clés Ratés et Ajoutés dont les valeurs correspondent à une liste
        de ce que Helixer a raté par rapport à l'annotation experte, et une liste pour ce qu'il a ajouté
        par rapport à l'annotation experte.
    """
    compteur_expert = 0
    compteur_Helixer = 0
    liste_missed = []
    liste_added = []
    h_feat = f"helixer_{feature}"

    liste_seq = sequenceDAO.getSequencesBySpeciesId(id) # Je choppe toutes les séquences pour une espèce
    for sequence in liste_seq : # Pour chacune des séquences, je vais chercher les deux annotations pour un même feat
        liste_feat = featureDAO.getAllFeatsForSequence(sequence['id'], feature) # Ouuuf j'accède à tous les feats
        liste_feat_Helixer = featureDAO.getAllFeatsForSequence(sequence['id'], h_feat) # Pareil pour tous les feat Helixer
        compteur_expert += len(liste_feat)
        compteur_Helixer += len(liste_feat_Helixer)
        # Les feat des experts non pris en compte par Helixer
        for feat in liste_feat :
            chevauchement = False
            for feat_h in liste_feat_Helixer:
                if are_overlapping(feat['start'], feat['stop'], feat_h['start'], feat_h['stop'], 0): 
                    # De cette manière là j'envisage toutes les possibilités de chevauchement
                    chevauchement = True
                    break
            if not chevauchement : # Si les feat se chevauchent pas, alors Helixer a raté le feat annoté par les experts
                liste_missed.append(feat)

        # Les feat Helixer en plus par rapport à l'annotation experte
        for feat_h in liste_feat_Helixer:
            chevauchement = False
            for feat in liste_feat:
                if are_overlapping(feat['start'], feat['stop'], feat_h['start'], feat_h['stop'], 0):
                    chevauchement = True
                    break
            if not chevauchement:
                liste_added.append(feat_h)
        #print(sequence['accession'])

    return {
        'Missed' : liste_missed,
        'Added' : liste_added,
        'Expert' : compteur_expert,
        'Helixer' : compteur_Helixer
    }


def obtain_both_feats(id : int, feature: str) -> list:
    """
        Prend en compte l'id de l'espèce et le feature d'intérêt.
        Renvoie une liste des features qui sont présents chez Helixer et chez les experts
    """
    liste_compare= []
    h_feat = f"helixer_{feature}"
    liste_seq = sequenceDAO.getSequencesBySpeciesId(id) # Je récupère toutes les séquences d'une espèce
    for sequence in liste_seq : # Pour chacune des séquences, je vais chercher les feats annotées de la même manière
        liste_feat = featureDAO.getAllFeatsForSequence(sequence['id'], feature) # Ouuuf j'accède à tous les feats
        liste_feat_Helixer = featureDAO.getAllFeatsForSequence(sequence['id'], h_feat) # Pareil pour tous les feat Helixer
        # Le feat des experts chevauchant celui de Helixer
        for feat in liste_feat :
            for feat_h in liste_feat_Helixer:
                if are_overlapping(feat['start'], feat['stop'], feat_h['start'], feat_h['stop'], 0):
                    # De cette manière là j'envisage toutes les possibilités de chevauchement
                    liste_compare.append((feat,feat_h)) # Si les feat se chevauchent, alors on peut alors vérifier si ils sont pareils
                    break # Seulement un feat expert chevauche un feat Helixer (pas de transcrit alternatif par Helixer)
    return liste_compare


def getFeatwithIdfAndType(id_feat, feattype, species_id):
    """
        Petite fonction pour éviter le crash du programme 
        lorsqu'un CDS n'a pas de mRNA correspondant (souvent isoforme).
    """
    try:
        return featureDAO.getFeatByIdfAndType(id_feat, feattype, species_id)
    except Exception:
        return None


def comparing_annotations(features1: list, features2: list, threshold: int = 0) -> dict:
    """
        Prend en deux listes de features.
        Renvoie un dictionnaire avec pour clés :
        - Added : Ce que features2 a de plus que features1
        - Missed : Ce que features2 n'a pas de features1
        - Common : Ce que features2 et features1 ont en commun     

        Le but ici est de pouvoir chercher des overlaps entre deux différentes annotations.  
        De cette manière plusieurs recherches peuvent être faites :
        - Helixer (CDS) VS Experts (CDS)
        - Helixer (mRNA) VS Experts (mRNA)
        - Helixer ADDED VS Experts (TEs)
        - Helixer (CDS) VS PFAM
    """
    # features1 = featureDAO.getAllFeatsForSpecies(34, "CDS")
    # features2 = featureDAO.getAllFeatsForSpecies(34, "helixer_CDS")
    liste_missed = []
    liste_added = []
    liste_common = []

    # On trie les features pour permettre une recherche plus efficace des overlaps
    features1 = sorted(features1, key=lambda x: x['start'])
    features2 = sorted(features2, key=lambda x: x['start'])

    # Qu'est-ce que features2 n'a pas de features1
    # Sweep line du coup... Parce que je veux optimiser mes fonctions
    j = 0
    for feat in features1:
        chevauchement = False
        while j < len(features2) and features2[j]['stop'] < feat['start']:
            j += 1 # Je saute jusqu'à trouver une vraie zone potentielle d'overlap
        k = j
        while k < len(features2) and features2[k]['start'] <= feat['stop']: # Je vais chercher tous les overlaps tel un gourmand
            if are_overlapping(feat['start'], feat['stop'], features2[k]['start'], features2[k]['stop'], threshold) and feat['name'] == features2[k]['name'] :
                liste_common.append((feat, features2[k]))
                chevauchement = True
                break
            k += 1
        if not chevauchement:
            liste_missed.append(feat)

    # Qu'est-ce que features2 a de plus par rapport à features1
    i = 0
    for feat_h in features2:
        chevauchement = False
        while i < len(features1) and features1[i]['stop'] < feat_h['start']:
            i += 1 # Potentiel man
        k = i
        while k < len(features1) and features1[k]['start'] <= feat_h['stop']:
            # En tant que bon kiffeur de SQL, j'ai réuni dans les features de TOUTES les séquences en une liste...
            # Alors il faut que je m'assure que les comparaisons se fassent sur le MEME chromosome !!
            if are_overlapping(feat_h['start'], feat_h['stop'], features1[k]['start'], features1[k]['stop'], threshold) and feat_h['name'] == features1[k]['name'] :
                chevauchement = True
                break
            k += 1
        if not chevauchement:
            liste_added.append(feat_h)

    return {
        'Missed': liste_missed,
        'Added': liste_added,
        'Common': liste_common
    }


def get_score(start1: int, stop1: int, start2: int, stop2: int) -> float:
    """
        Prend en paramètre les coordonnées des exons, et renvoie un score d'identité
    """
    ratio1 = stop1 - start1 
    ratio2 = stop2 - start2 
    return min(ratio1,ratio2) / max(ratio1,ratio2)


def compare_both_feats(feat1: dict, feat2: dict) -> dict :
    """
        Prend en paramètre deux features.
        Renvoie un rapport de type dict pour résumer la comparaison des deux features.
    """
    if feat1['complement'] if 'complement' in feat1.keys() else feat1['sens'] == feat2['complement'] if 'complement' in feat2.keys() else feat2['sens']: 
        # On s'assure qu'on est sur le même sens
        # Je pense pas qu'Helixer aurait pu se tromper sur le sens mais on sait jamais
        # Par ailleurs Location est sous cette forme [liste_start, liste_stop], on a toujours le même nombre de start que de stop
        ref = feat1['location']
        helixer = feat2['location']

        # Exons
        calculs = []
        compteur_c = 0
        # Je vais chercher chaque exon se chevauchant histoire de pouvoir calculer le % d'identité
            # ref[0] : liste de tous les starts de la réf
            # ref[1] : liste de tous les stops de la ref
        for i in range(len(ref[0])):
            for j in range(len(helixer[0])):
                if are_overlapping(ref[0][i], ref[1][i], helixer[0][j], helixer[1][j], 0):
                    compteur_c += 1
                    calculs.append(get_score(ref[0][i], ref[1][i],helixer[0][j], helixer[1][j])) 
                    # Normalement ça donne un résultat entre 0 et 1 pour chaque exon
        
        for c in range(max(max(len(ref[0]), len(helixer[0])), compteur_c)-min(max(len(ref[0]), len(helixer[0])), compteur_c)):
            calculs.append(0) # Les exons manquants chez une annotation donnent un score nul
        if len(calculs) != 0:
            score = round((sum(calculs)/len(calculs))*100, 2) # On fait une moyenne pour le % d'identité
        else:
            score = 0

        # Introns 
        introns_ref, introns_hel  = [], []
        for i in range(len(ref[0]) -1): # Les introns de référence
            introns_ref.append((ref[1][i], ref[0][i+1]))
        for i in range(len(helixer[0]) -1): # Les introns de Helixer
            introns_hel.append((helixer[1][i], helixer[0][i+1]))
        
        compteur_i = 0 # Compteur d'introns annotés des deux côtés
        for intron in introns_ref:
            for intron_h in introns_hel :
                if intron == intron_h:
                    compteur_i +=1
                    break
        
        # Les types de différences entre les annotations
        if score == 100:
            if ref[0][0] != helixer[0][0] or ref[1][-1] != helixer[1][-1]: # Start ou Stop différents
                type_c = "~"
            else:
                type_c = "="
        else : # Score d'identité en dessous de 100%
            # On va comparer en fonction des introns recouverts ou matchés
            if len(introns_ref) == len(introns_hel) and compteur_i == len(introns_ref):
                type_c = "~" # Les exons ne sont pas pareils mais les chaines d'introns le sont
            elif compteur_i == len(introns_hel) and compteur_c == len(helixer[0]): 
                type_c = "c" # Les introns de Helixer sont contenus dans l'annotation de référence
            elif compteur_i == len(introns_ref) and compteur_c == len(ref[0]): 
                type_c = "k" # Les introns de Helixer contiennent l'annotation de référence
            elif compteur_c == max(len(ref[0]), len(helixer[0])) :
                type_c = "m" # Tous les introns de référence sont recouverts / matchés
            elif compteur_c >= 1:
                type_c = "j" # Au moins un match entre les exons
            else :
                type_c = "n" # Tous les introns de référence ne sont pas matchés
            
        # On rajoute le tout dans le rapport
        idf = 'id_feat'
        return {'Chromosome': feat1['name'], 'Gène Ref': feat1['id_feat'] if not "." in feat1['id_feat'][:-2] else feat1['id_feat'][:-2],
                'Gène Ref ID': feat1['id'], 'Helixer': feat2[idf],'Helixer ID': feat2['id'],'Isoforme Ref': feat1['id_feat'], 'Isoforme Ref ID': feat1['id']
                , 'Isoformes': [] ,'Identité': score, 'Type': type_c}


def compare_all_feats(dir, seuil : int, type_f: str, overlaps: list[tuple]) -> list[dict]:
    """
        Prend en compte un % seuil et un type de superposition, et la liste des features en commun.  
        Renvoie une liste de dictionnaire dans lequel on retrouve :  
        - le feature de référence  
        - le feature de comparaison
        - le score d'identité
        - le type d'égalité  

        La liste est filtrée selon le seuil et le type mis en paramètre.  
        Plus le score est proche de 100, plus les features sont similaires.
    """
    seuil = max(min(seuil, 100), 0) # Pas en dessous de 0 et pas au dessus de 100
    type_f = type_f if type_f in "=~ckmjn" else "="  # Un type par défaut si le type est mal précisé

    rapport = [] # Rapport non filtré
    rapport1 = [] # Rapport filtré

    # Etablissement d'un rapport sans filtre
    for comparaison in overlaps :
        resultat = compare_both_feats(comparaison[0], comparaison[1])
        if resultat: # Peut renvoyer None si les deux feats n'ont pas le même sens !
            rapport1.append(resultat)

    # Filtrage
    liste_opps = set() # Si c'est dedans, on en veut pas pour plus tard
    file = "StatsSAMECDS.txt" if type_f == "=" else "StatsALMOSTSAMECDS.txt"
    with open(f"{dir}/{file}", "w") as stats:
        for element in rapport1:
            if element['Helixer'] in liste_opps:
                continue  # BYE BYE LE DOUBLON
            isoformes = []
            # On cherche tous les isoformes qui chevauchent un Helixer
            for element2 in rapport1:
                if element['Helixer'] == element2['Helixer']:
                    liste_opps.add(element2['Helixer'])
                    isoformes.append(element2) # Obtention de l'isoforme
            result = max(isoformes, key=lambda d: d['Identité']) # Le meilleur transcrit alternatif
            isoformes.remove(result)
            result['Isoformes'] = [iso['Isoforme Ref'] for iso in isoformes] # Une liste des transcrits alternatifs

            # Liste des features selon les filtres
            if type_f != "=" :
                if result['Identité'] >= seuil and result['Type'] == type_f:
                    rapport.append(result)
            else :
                if result['Type'] == type_f: # Quand on est dans le cas "=", tous les scores sont à 100
                    rapport.append(result)

            stats.write(f"{result['Type']}\n")
    return rapport


def comparing_UTR(comparaison_cds: dict, species_id) -> dict:
    """
        Prend en paramètre deux features et renvoie les différences dans les régions 5' UTR et 3'UTR
    """
    feat1 = comparaison_cds['Isoforme Ref']
    expert_CDS= getFeatwithIdfAndType(feat1, "CDS", species_id)
    expert_mRNA = getFeatwithIdfAndType(feat1, "mRNA", species_id)

    feat2 = comparaison_cds['Helixer']
    helixer_CDS = getFeatwithIdfAndType(feat2, "helixer_CDS", species_id) 
    helixer_mRNA = getFeatwithIdfAndType(feat2, "helixer_mRNA", species_id)
    
    # On compare les deux ARNm et on obtient leurs régions UTR
    comparaison = compare_both_feats(expert_mRNA, helixer_mRNA)
    
    expert5 = sequenceDAO.getSequenceDNAFromTo(expert_mRNA['id_seq'], expert_mRNA['start'], expert_CDS['start'])
    expert3 = sequenceDAO.getSequenceDNAFromTo(expert_mRNA['id_seq'], expert_CDS['stop'], expert_mRNA['stop'])
    helixer5 = sequenceDAO.getSequenceDNAFromTo(helixer_mRNA['id_seq'], helixer_mRNA['start'], helixer_CDS['start'])
    helixer3 = sequenceDAO.getSequenceDNAFromTo(helixer_mRNA['id_seq'], helixer_CDS['stop'], helixer_mRNA['stop'])
    if expert_CDS['complement'] == -1:
        expert5 = rev_compl(expert5)
        expert3 = rev_compl(expert3)
        helixer5 = rev_compl(expert5)
        helixer3 =  rev_compl(expert3)

    expert5 = "" if len(expert5) == 1 else expert5
    expert3 = "" if len(expert3) == 1 else expert3

    from difflib import SequenceMatcher
    # Pas de substitution jusqu'à preuve du contraire, que des indels
    def diff_utrs(a, b):
        matcher = SequenceMatcher(None, a, b)
        diff_a = []
        diff_b = []
        for tag, i1, i2, j1, j2 in matcher.get_opcodes():
            if tag == "equal":
                continue

            diff_a.append(a[i1:i2])
            diff_b.append(b[j1:j2])

        return "".join(diff_a).replace("".join(diff_b), ""), "".join(diff_b).replace("".join(diff_a), "")

    # On cherche à avoir les différences entre 5'UTR 
    utr5h, utr5e = diff_utrs(helixer5, expert5)
    # Ici en 3'UTR
    utr3h, utr3e = diff_utrs(helixer3, expert3)

    return {'Chromosome': expert_mRNA['name'], 'Gène': comparaison_cds['Gène Ref'] ,'Référence': feat1, 'Helixer': feat2, 
        'Isoformes': comparaison_cds['Isoformes'],'Identité': comparaison['Identité'], 'Type': comparaison['Type'], 
        '5 UTR Helixer': utr5h, '3 UTR Helixer': utr3h, '5 UTR Expert': utr5e, '3 UTR Expert': utr3e, '5Expert Len' : len(expert5), 
        '3Expert Len' : len(expert3), '5Helixer Len' : len(helixer5), '3Helixer Len' : len(helixer3)}


def compare_every_UTR(liste_equals: list[dict], species_id) -> list[dict]:
    """
        Prend en paramètre une liste de features dont les annotations sont exactement
        ou quasiment identiques. 
        Renvoie un rapport sous forme de liste avec les différences dans les régions UTR.
    """
    liste_all = []
    for element in liste_equals :
        liste_all.append(comparing_UTR(element, species_id))
    return liste_all


def protein_analysis(type : str, liste_added : list[dict], fasta_file: str, dir: str, dl_new_ver : bool):
    """
        Prend en compte la liste des gènes ajoutés par Helixer et un fichier .fasta contenant toutes ses protéines prédites.  
        Retourne un fichier tabulé avec  :
        - L'id_feat du CDS
        - La taille du CDS
        - Le nombre d'exons du CDS
        - Le nom de la protéine générée (Son alias trouvé dans Uniprot, sinon celui dans le fichier .fasta)
        - L'annotation fonctionnelle de cette protéine si existante.  
        !! Il faut avoir la dernière version de la base Uniprot:  
        Le téléchargement et l'installation de la base se font directement si dl_new_ver est True.  
        !! Et oui encore un disclaimer... : Il faut avoir Blast installé localement !
    """
    if dl_new_ver: # Si on veut update la base ou l'installer de zéro
        # J'installe et je dezip tout ce dont j'ai besoin, je crée la base blast et je supprime les .gz alors inutiles
        os.makedirs("Projet/Uniprot", exist_ok=True)
        os.system("rm -f Projet/Uniprot/*")
        os.system("wget -O 'Projet/Uniprot/swissprot_db_ref.xml.gz' https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz")
        os.system("gunzip -f 'Projet/Uniprot/swissprot_db_ref.xml.gz'")
        os.system("rm 'Projet/Uniprot/swissprot_db_ref.xml.gz'")
        os.system("wget -O 'Projet/Uniprot/swissprot_db.fasta.gz' 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'")
        os.system("gunzip -f 'Projet/Uniprot/swissprot_db.fasta.gz'")
        os.system("rm 'Projet/Uniprot/swissprot_db.fasta.gz'")
        os.system(f"makeblastdb -in 'Projet/Uniprot/swissprot_db.fasta' -dbtype prot -parse_seqids -out 'Projet/Uniprot/swissprot'")

    # OKAY GARMIN !! BLAST MOI TOUT CA ET EXPLOSE MON PC
    os.system(f"blastp -query {fasta_file} -db 'Projet/Uniprot/swissprot' -out {dir}/results_{type}.tsv \
              -outfmt '6 qseqid sseqid pident length qcovs evalue bitscore' -matrix BLOSUM45 \
              -num_threads 12 -evalue 1e-5 -max_target_seqs 5")

    # Il est temps de créer le fichier
    with open(f"{dir}/Protein_Analysis_{type}.tsv", "w") as proteins, open(f"{dir}/results_{type}.tsv", "r") as results:
            proteins.write("Gene\tCDSLength\tExonsNb\tProtein\tGO Annot\tFunction\tALength\tQCover\tScore\tEvalue\tPIdent\tSpecies\n")
            import xml.etree.ElementTree as ET # J'importe qu'ici car on utilisera pas souvent protein_analysis()...
            NS = {"up": "https://uniprot.org/uniprot"}

            uniprot = {}
            for event, elem in ET.iterparse("Projet/Uniprot/swissprot_db_ref.xml", events=("end",)):
                if elem.tag != "{https://uniprot.org/uniprot}entry":
                    continue
                # ACCESSION
                acc = elem.findtext("up:accession", namespaces=NS)
                # NAME
                name = elem.findtext("up:protein/up:recommendedName/up:fullName", namespaces=NS)
                # FULL NAME IF NO NAME...
                if name is None:
                    name = elem.findtext("up:protein/up:submittedName/up:fullName", namespaces=NS)
                # SPECIES
                species = elem.findtext("up:organism/up:name[@type='scientific']", namespaces=NS)
                # FUNCTION
                function = ""
                for c in elem.findall("up:comment", namespaces=NS):
                    if c.attrib.get("type") == "function":
                        function = " ".join(c.itertext()).strip()
                # GO (Coucou Carène)
                go_terms = [
                    db.attrib["id"]
                    for db in elem.findall("up:dbReference", namespaces=NS)
                    if db.attrib.get("type") == "GO"
                ]

                uniprot[acc] = {
                    "name": name,
                    "species" : species,
                    "function": function,
                    "go": go_terms
                }

                elem.clear()
            
            # Manipuler les dict me permet ici d'éviter des boucles imbriquées qui tueraient le PC..
            # Je fais en sorte que best[0] fasse : ['Protein'] = {subject, bitscore...} tout en retenant le meilleur hit de cette protéine 

            # Passons aux résultats blast
            best = {} # Plusieurs hits par protéine, je vais juste garder le meilleur (je me base sur la meilleure evalue et le bitscore)
            for line in results:
                q, s, pid, alength, qcovs, evalue, bitscore = line.strip().split("\t")
                # Merci Franck (les noms des protéines query ont pour format 'id_feat'_'id')
                s = s.split("|")[1]
                evalue = float(evalue)
                bitscore = float(bitscore)
                hit = (bitscore, -evalue, s)
                if q not in best or hit > best[q]["key"]:
                    best[q] = {
                        "sseqid": s,
                        "bitscore": bitscore,
                        "pident" : pid,
                        "evalue": evalue,
                        "Alength": alength,
                        "Qcover" : qcovs,
                        "key": hit
                    }
            # On peut enfin écrire 
            for added in liste_added: # added = {'id_feat', 'length', 'exons'...}
                label, func, go, score, identity, evalue, Alength, QCover, species = "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"   # Tu chantes ou quoi
                if f"{added['id_feat']}_{added['id']}" in best :
                    prot = best[f"{added['id_feat']}_{added['id']}"]
                    subject, score, identity, evalue, Alength, QCover = prot['sseqid'], prot['bitscore'], prot['pident'], prot['evalue'], prot['Alength'], prot['Qcover']
                    if subject in uniprot: 
                        label, func, go, species = uniprot[subject]['name'], uniprot[subject]['function'], ";".join(uniprot[subject]['go']), uniprot[subject]['species']
                        # Gene\     Length\    ExonsNb\    Protein\   GO Annot\     Function\   Alength\    QCover\   Score\   Evalue\   PIdent        
                proteins.write(f"{added['id_feat']}\t{added['length']}\t{added['exons']}\t{label}\t{go}\t{func}\t{Alength}\t{QCover}\t{score}\t{evalue}\t{identity}\t{species}\n")


def TE_blast(liste_added : list[dict], dir : str, name : str, dl: bool = False):
    """
        Prends en compte la liste des gènes rajoutés par Helixer,  
        Execute un BlastN pour détecter les matchs avec les ET du génome donné.  
        Renvoie un fichier tabulé avec :  
        - Le nom du gène ajouté
        - Ses coordonnées dans le génome
        - Les ET avec qui il match
        - Les informations relatives à ces ET
    """
    # Les fichiers .gff3 pour plus tard... avec les coordonnées de TEs
    gff_original = f"Projet/Assemblies/TES/{name}.gff3"
    gff_corrige = f"{dir}/TES/{name}.fixed.gff3"

    # Création d'un fichier .fasta avec toutes les séquences des gènes ajoutés par Helixer...
    with open(f"{dir}/seqs/added_TE_DNA.fasta", "w") as seqFile:
            for added in liste_added:
                fullFeature = featureDAO.getFeatureForId(added['id'])
                sequence = sequenceDAO.getSequenceDNAFromTo(fullFeature['sid'], fullFeature['start'], fullFeature['stop'])
                seqFile.write(">" + added['id_feat'] + f"|{added['id']}" + "\n")
                seqFile.write(sequence + "\n")

    if dl:
        os.makedirs(f"{dir}/TES", exist_ok=True)
        os.system(f"rm -f {dir}/TES/*")

        # On lance bedtools pour faire le lien .gff3 ~ .fasta
        # En gros, je veux créer ma propre base de ET en extrayant leurs séquences correspondantes via leurs coordonnées depuis le fichier .fasta 

        # Remplace tout ce qui n'est pas 'chr' par 'chr' au tout début de chaque ligne
        os.system(
            f"""awk -F'\\t' 'BEGIN{{OFS="\\t"}}
            NR==1 {{print; next}}
            {{
                if ($1 ~ /^[Ss]caffold/) {{
                    print
                    next
                }}

                if ($1 ~ /[Mm ]$/) {{
                    $1 = "chrM"
                }}
                else if ($1 ~ /[Cc ]$/) {{
                    $1 = "chrC"
                }}
                else {{
                    sub(/^[^0-9]*0*/, "chr", $1)
                }}

                print
            }}
            ' {gff_original} > {gff_corrige}"""
        )       
        
        # Création d'un fichier .bed avec les infos du .gff3 créé avant
        os.system(f"""
        awk 'BEGIN{{OFS="\\t"}}
        {{
            split($9,a,";");
            sub("ID=","",a[1]);
            print $1,$4-1,$5,a[1],".",$7;
        }}' {gff_corrige} > {dir}/TES/{name}.bed
        """)

        # On normalise les noms des chromosomes même dans le .fa
        os.system(
            f"""sed -E '
            /^>/ {{
                s/^>([^ ]*)[Mm](.*)$/>chrM\\2/
                s/^>([^ ]*)[Cc](.*)$/>chrC\\2/
                s/^>[^0-9 ]*0*([0-9]+)(.*)$/>chr\\1\\2/
            }}
            ' Projet/Assemblies/{name}.fa > Projet/Assemblies/{name}1.fa"""
        )
        os.system(
            f"""sed -i 's/^>chrChrM/>chrM/' Projet/Assemblies/{name}1.fa""")

        # On prend les séquences correspondantes aux features présents dans le .bed     
        os.system(f"samtools faidx Projet/Assemblies/{name}1.fa")
        os.system(f"bedtools getfasta -fi Projet/Assemblies/{name}1.fa -bed {dir}/TES/{name}.bed -fo {dir}/TES/{name}_TEs.fasta -nameOnly")

        # Roooh mais il saoule là à faire des blastdb
        cmd = (
            f"sed 's/::/_/g' {dir}/TES/{name}_TEs.fasta | "
            "awk '/^>/{"
            "h=$0; sub(/^>/,\"\",h); "
            "gsub(/_Chr[0-9A-Za-z_]*:[0-9]+-[0-9]+/,\"\",h); "
            "gsub(/__+/,\"_\",h); "
            "h=substr(h,1,45); "
            "print \">\"h\"_\"++i; next}1' "
            f"> {dir}/TES/blast_TEs1.fasta"
        )
        os.system(cmd)

    os.system(f"makeblastdb -in '{dir}/TES/blast_TEs1.fasta' -dbtype nucl -parse_seqids -out '{dir}/TES/blast_TEs1'")

    # OKAY GARMIN !! BLAST MOI TOUT CA ET EXPLOSE MON PC
    os.system(f"blastn -query {dir}/seqs/added_TE_DNA.fasta -db '{dir}/TES/blast_TEs1' -out {dir}/TEs_results.tsv \
              -outfmt '6 qseqid sseqid pident length qcovs evalue bitscore'")
    
    # Dict pour pouvoir faire le lien .gff3 ~ résultats Blast
    te_index = {}
    with open(gff_corrige, "r") as gff_file:
        for line in gff_file:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            
            # Extraction des infos utiles de la ligne
            chrom, _, feature_type, start_pos, stop_pos, _, strand, _, attributes = parts
            te_id = attributes.split(";")[0].replace("ID=", "")
            te_index[te_id] = (chrom, feature_type, start_pos, stop_pos, strand)


    # Passons au fichier tabulé
    with open(f"{dir}/found_TEs.tsv", "w") as TEs, open(f"{dir}/TEs_results.tsv", "r") as results:
            TEs.write("Gene\tChr\tStrand\tStart\tStop\tTE_ID\tTEChr\tTEStrand\tTEType\tTEStart\tTEStop\tALength\tQCover\tScore\tEvalue\tPIdent\n")
            for line in results: # added = {'id_feat', 'length', 'exons'...}
                q, s, pid, alength, qcovs, evalue, bitscore = line.strip().split("\t")
                id = q.split("|")[1]
                true_feat = featureDAO.getFeatureForId(id)
                sens = "+" if true_feat['sens'] == 1 else "-" if true_feat['sens'] == -1 else "."
                gene, chromosome, start, stop, strand = true_feat['idf'], true_feat['chromosome'], true_feat['start'], true_feat['stop'], sens
                
                te_found = s.split("_")[0]
                if te_found in te_index:
                    chr_name, te_type, testart, testop, testrand = te_index[te_found]
                else:
                    chr_name, te_type, testart, testop, testrand = "Unknown", "Unknown", "NA", "NA", "."
                label, TEchr, TEtype, TEStart, TEStop, TEstrand, Alength, QCover, score, eval, identity = s, chr_name, te_type, testart, testop, testrand, alength, qcovs, bitscore, evalue, pid
                TEs.write(f"{gene}\t{chromosome}\t{strand}\t{start}\t{stop}\t{label}\t{TEchr}\t{TEstrand}\t{TEtype}\t{TEStart}\t{TEStop}\t{Alength}\t{QCover}\t{score}\t{eval}\t{identity}\n")


def bed_to_gff_compare(dir, input, compare_to):
    """
        Prends un .BED de CDS et le convertit en .gff3
        Le compare à un .gff3 existant via bedtools.
        (Ici on l'utilise pour les TEs)
    """
    output = input.replace(".bed", ".gff3")
    
    # On convertit le .BED query en .gff3
    os.system(f"""awk 'BEGIN {{
        OFS="\\t";
        print "##gff-version 3"
    }}
    {{
        chr=$1

        if (chr ~ /^[Ss]caffold/) {{
        }}
        else if (chr ~ /[Mm ]$/) {{
            chr = "chrM"
        }}
        else if (chr ~ /[Cc ]$/) {{
            chr = "chrC"
        }}
        else {{
            sub(/^[^0-9]*0*/, "chr", chr)
        }}

        strand = ($6 == "-1") ? "-" : "+"

        print chr, "Helixer", "CDS", $2+1, $3, ".", strand, ".", "ID=" $4
    }}' {input} > {output}""")

    # Le .gff3 sujet qu'on utilise comme élément de comparaison est trié : on conserve que les chromosomes annotés 
    os.system(f"""
    awk 'BEGIN{{OFS="\\t"}}
    $1 ~ /^chr/ || $1 ~ /^[0-9]+$/ {{print > "{dir}/TES/chromosomes.gff3"}}
    $1 !~ /^chr/ && $1 !~ /^[0-9]+$/ {{print > "{dir}/TES/scaffolds.gff3"}}' {compare_to}
    """)
    
    # On fait 2 fichiers d'overlaps entre le query et le sujet
    os.system(f"bedtools intersect -a {output} -b {dir}/TES/chromosomes.gff3 > {input.replace(".bed", "_output.gff3")}") # 1bp suffit
    os.system(f"bedtools intersect -a {output} -b {dir}/TES/chromosomes.gff3 -f 0.5 -r > {input.replace(".bed", "_output_cover.gff3")}")
    # Le query et le sujet se recouvrent au moins de moitié

    os.system(f"""
        awk 'NR==FNR {{
            if ($0 ~ /ID=/) {{
                split($0, f, "ID=")
                ids[f[2]] = 1
            }}
            next
        }}
        FNR==1 {{print; next}}
        ($1 in ids)
        ' {input.replace(".bed", "_output_cover.gff3")} \
        {dir}/found_TEs.tsv \
        > {dir}/filtered_TEs.tsv
    """)
    # J'y pense mais en vrai Bedtools aurait pu servir au tout début...
    # Mais le format pour les analyses d'après nécessitaient ma pipeline ! ...

    os.system(f"Rscript Projet/Rapport/CodesFig/TEsFINALBOSS.r {dir}")


def TE_density(dir, name2, name):
    """
        Renvoie la densité des TEs à proximité de gènes.
    """
    # On convertit le .gff3 en .bed
    gff_file = f"{dir}/TES/{name2}.fixed.gff3"
    output_bed = f"{dir}/Bed/TE.bed"

    os.system(f"cut -f1,2 Projet/Assemblies/{name2}1.fa.fai > {dir}/TES/genome.chrom.sizes")
    if name2 == "ARAPORT11":
        os.system(f"""
        awk -F'\\t' '
        BEGIN{{OFS="\\t"}}
        $3=="transposable_element" || $3=="transposable_element_gene" {{
            print $1,$4-1,$5
        }}' {gff_file} > {output_bed}
        """)
    else :
        os.system(f"""
        awk -F'\\t' '
        BEGIN{{OFS="\\t"}}{{
            print $1,$4-1,$5
        }}' {gff_file} > {output_bed}
        """)

    # On prend le .bed et on met en place les fenêtres
    input_bed0 = f"{dir}/Bed/{name}_commonEXP.bed"
    input_bed1 = f"{dir}/Bed/{name}_commonHEL.bed"
    input_bed2 = f"{dir}/Bed/{name}_added.bed"
    input_bed3 = f"{dir}/Bed/{name}_missed.bed"
    inputs = [input_bed0, input_bed1, input_bed2, input_bed3]

    output_bed0 = f"{dir}/Bed/commonEXP_window.bed"
    output_bed1 = f"{dir}/Bed/commonHEL_window.bed"
    output_bed2 = f"{dir}/Bed/added_window.bed"
    output_bed3 = f"{dir}/Bed/missed_window.bed"
    outputs = [output_bed0, output_bed1, output_bed2, output_bed3]


    genome_sizes = f"{dir}/TES/genome.chrom.sizes"
    window_size = 2000

    for i in range(4): # Pour tous nos cas 
        os.system(f"""
        bedtools slop \
            -i {inputs[i]} \
            -g {genome_sizes} \
            -b {window_size} \
            > {outputs[i]}
        """)
        # Le but est de chercher les overlaps avec la fenêtre donnée
        os.system(f"""
        bedtools intersect \
            -a {outputs[i]} \
            -b {output_bed} \
            -wo \
            > {outputs[i].replace("window.bed", "overlap.txt")}
        """)       
    os.system(f"Rscript Projet/Rapport/CodesFig/densityTE.r {dir}/Bed/")


def PLM_analysis(dir, get_all_feats : list[dict], JASPAR : bool = False, john_mail: str = "67@gmail.com"):
    """
        Ouvre PLMViewer et cherche les PLM chez Helixer.
    """
    os.makedirs(f"{dir}/PLM", exist_ok=True)

    with open(f"{dir}/PLM/HelPLM.fasta", "w") as seqFile:
        for hel in get_all_feats:
            fullFeature = featureDAO.getFeatureForId(hel['id'])
            # TSS : -1000; TSS ;+500
            start = fullFeature['start']
            stop = fullFeature['stop']
            seq = sequenceDAO.getSequencesById(fullFeature['sid'])['length']

            start1 = start - 1000 if start - 1000 >= 0 else 1
            stop1 = stop + 500 if stop + 500 <= seq else seq
            sequence = sequenceDAO.getSequenceDNAFromTo(fullFeature['sid'], start1, stop1)
            if fullFeature['sens'] == -1:
                sequence = rev_compl(sequenceDAO.getSequenceDNAFromTo(fullFeature['sid'], start1, stop1))
            seqFile.write(">" + hel['id_feat'] + f"|{hel['id']}" + "\n")
            seqFile.write(sequence + "\n") 

    # Il faut installer Selenium pour ce faire !
    from selenium import webdriver
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.common.exceptions import TimeoutException
    from selenium.webdriver.common.keys import Keys
    from selenium.webdriver.chrome.options import Options
    import time

    options = Options()
    options.add_experimental_option("detach", True)
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    driver = webdriver.Chrome(options=options)
    try:
        driver.get("https://plmview.ips2.universite-paris-saclay.fr/?go=new")
    except TimeoutException:
        driver.refresh()

    plm_coustou = "Projet/PLM_Vero.txt"
    # Upload Helixer
    driver.find_element(By.CSS_SELECTOR,"label[for='_easyui_radiobutton_6']").click()
    # Je cherche le bouton pour upload
    seq_input1 = driver.find_element(By.CSS_SELECTOR, "input[id='filebox_file_id_2']")
    seq_input1.send_keys(os.path.abspath(f"{dir}/PLM/HelPLM.fasta"))

    # PLM truc bidule qu'on doit upload
    if JASPAR :
        arrows = driver.find_elements(By.CSS_SELECTOR, "a.combo-arrow")
        arrows[2].click()  
        time.sleep(0.5)
        elem = driver.find_element(By.ID, "motifDefaultFile_easyui_combobox_i3_2")
        elem.click()
    else:
        driver.find_element(By.CSS_SELECTOR, "label[for='_easyui_radiobutton_8']").click()
        seq_input2 = driver.find_element(By.CSS_SELECTOR, "input[id='filebox_file_id_3']")
        seq_input2.send_keys(os.path.abspath(plm_coustou))
        
    # Nom et mail et run (j'espère)
    driver.find_element(By.CSS_SELECTOR, "input[id='_easyui_textbox_input9']").send_keys("crakshay") # Oui bah je sais pas qui mettre
    driver.find_element(By.CSS_SELECTOR, "input[id='_easyui_textbox_input10']").send_keys(john_mail) # Bon oui voilà c'est mon mail perso
    time.sleep(15)
    driver.find_element(By.XPATH, "//a[span[span[text()='Run PLM']]]").click() # RUN
    time.sleep(2)
    


def exporting_results(species_id: int, dir: str, do_all: bool, UTR_analysis: bool, overall_RNA: bool,
                      categories_analysis: bool, TE_analysis: bool, 
                      plm_view: bool, john_mail: str, download_db: bool,
                      protein_comparison: bool, busco: bool, generate_seqs: bool):
    """
        Execute l'analyse entière pour une espèce donnée.
    """
    def progress_bar(progress, total, length=40):
        """
            Une barre de progression faite à l'arrache...
        """
        percent = progress / total
        filled = int(length * percent)
        bar = "\033[32m█\033[0m" * filled + "-" * (length - filled)
        print(f"\r[{bar}] ({percent*100:.1f}%)", end=" |")

    if do_all: # Tout faire
        generate_seqs, categories_analysis, TE_analysis, protein_comparison, plm_view, UTR_analysis, overall_RNA, busco = True, True, True, True, True, True, True, True
    nb_functions = sum((UTR_analysis, overall_RNA,
                      categories_analysis, TE_analysis, 
                      plm_view, protein_comparison, busco, generate_seqs)) +2
    completed = 0

    # On entame l'analyse 
    name = speciesDAO.getSpeciesById(species_id)['name'].replace(" ", "_")
    name2 = speciesDAO.getSpeciesById(species_id)['name']
    os.makedirs("Projet/Rapport", exist_ok=True)
    os.makedirs(dir, exist_ok = True)
    os.makedirs(f"{dir}/seqs", exist_ok = True)
    os.makedirs(f"{dir}/figs", exist_ok = True)
    os.makedirs(f"{dir}/BlastAnalysisADDED", exist_ok = True)
    os.makedirs(f"{dir}/Stats", exist_ok = True)
    os.makedirs(f"{dir}/Bed", exist_ok = True)
    os.makedirs(f"{dir}/UTRAnalysis", exist_ok = True)
    os.makedirs(f"{dir}/Analysis", exist_ok = True)

    print(f"Analyse de {name2}.")
    progress_bar(1, 100)

    #not_in_common = find_missingAnnotationV1(species_id, "CDS")
    experts = featureDAO.getAllFeatsForSpecies(species_id, "CDS")
    helixers = featureDAO.getAllFeatsForSpecies(species_id, "helixer_CDS")
    expertsRNA = featureDAO.getAllFeatsForSpecies(species_id, "mRNA")
    helixerRNA = featureDAO.getAllFeatsForSpecies(species_id, "helixer_mRNA")
    not_in_common = comparing_annotations(experts, helixers)
    commons = [duos[1] for duos in not_in_common['Common']]
    helixer_commons = [duos[0] for duos in not_in_common['Common']]
    pfam = featureDAO.getAllFeatsForSpecies(species_id, "PFAM")

    completed +=1
    progress_bar(completed/nb_functions*100, 100) 

    ## On les intègre à des fichiers différents
    # Ici pour les stats
    with open(f"{dir}/Stats/Commons_{name}.txt", "w") as dual:
        dual.write("Reference\tHelixer\tLength\tNbExons\n")
        for duo in not_in_common['Common']:
            other = duo[1]
            dual.write(f"{duo[0]['id_feat']}\t{other['id_feat']}\t{other['length']}\t{other['exons']}\n")
    with open(f"{dir}/Stats/Added_{name}.txt", "w") as more:
        more.write("Helixer\tLength\tNbExons\n")
        for duo in not_in_common['Added']:
            more.write(f"{duo['id_feat']}\t{duo['length']}\t{duo['exons']}\n")
    with open(f"{dir}/Stats/Missed_{name}.txt", "w") as less:
        less.write("Reference\tLength\tNbExons\n")
        for duo in not_in_common['Missed']:
            less.write(f"{duo['id_feat']}\t{duo['length']}\t{duo['exons']}\n")
    
    # On passe aux fichiers BED
    with open(f"{dir}/Bed/{name}_missed.bed", "w") as missed:
        for element in not_in_common['Missed']:
            m = f"{normalize_chr(element['name'])}\t{element['start']}\t{element['stop']}\t{element['id_feat']}\t0\t{element['complement']}\n"
            missed.write(m)
    with open(f"{dir}/Bed/{name}_added.bed", "w") as added:
        for element in not_in_common['Added']:
            a = f"{normalize_chr(element['name'])}\t{element['start']}\t{element['stop']}\t{element['id_feat']}\t0\t{element['complement']}\n"
            added.write(a)
    with open(f"{dir}/Bed/{name}_commonHEL.bed", "w") as comm:
        for element in commons:
            c = f"{normalize_chr(element['name'])}\t{element['start']}\t{element['stop']}\t{element['id_feat']}\t0\t{element['complement']}\n"
            comm.write(c)
    with open(f"{dir}/Bed/{name}_commonEXP.bed", "w") as comm1:
        for element in helixer_commons:
            c = f"{normalize_chr(element['name'])}\t{element['start']}\t{element['stop']}\t{element['id_feat']}\t0\t{element['complement']}\n"
            comm1.write(c)

    # Big résumé
    with open(f"Projet/Rapport/resume.txt", "a") as g:
        if species_id == 20 :
            name2 = "ARAPORT11"
        if species_id == 34 :
            name2 = "TAIR12"
        if species_id == 23:
            name2 = "B73"
        g.write(f"{name2}\t{len(not_in_common['Missed'])}\t{len(not_in_common['Added'])}\t{len(not_in_common['Common'])}\t{len(experts)}\t{len(helixers)}\t{len(expertsRNA)}\t{len(helixerRNA)}\n")

    completed +=1
    progress_bar(completed/nb_functions*100, 100) 

    # Protéines pour added, missed et commons
    def fasta_files(type, set, type_feat = "", versus = False):
        if type_feat != "mRNA" :
            os.makedirs(f"{dir}/seqs/{type}", exist_ok = True)
            for gene in set:
                prot = calcProtein(gene['id'])

                protName = gene['id_feat'].split('.')[0]
                logger.debug(protName)
                logger.debug(prot)

                seqName = gene['id_feat']       
                seqName += "_"+str(gene['id']) if not versus else ""

                # Création du fichier de sequence
                seqFile = open(os.path.join(f"{dir}/seqs/{type}", gene['id_feat'].replace("(", '_').replace(")","_")+"_"+str(gene['id'])+".fasta"), "w")
                seqFile.write(">"+seqName+"\n")
                for seq in prot:
                    logger.debug(seq)
                    seqFile.write(seq+"\n")
                seqFile.close()
            os.system(f"find {dir}/seqs/{type}/ -maxdepth 1 -name '*.fasta' -exec cat {{}} + > {dir}/all_proteins_{name}_{type}.fasta")

    if generate_seqs: # Génération des fichiers .fasta protéiques pour chaque catégorie
        fasta_files("Added", not_in_common['Added'])
        fasta_files("Missed", not_in_common['Missed'])
        fasta_files("CommonREF", helixer_commons)
        fasta_files("CommonHELIXER", commons)
    
        completed +=1
        progress_bar(completed/nb_functions*100, 100) 

    if species_id == 34: # Analyse spéciale TAIR12
        # ARAPORT11 VS TAIR12
        # Pour ARAPORT11
        experts11 = featureDAO.getAllFeatsForSpecies(20, "CDS")
        experts11RNA = featureDAO.getAllFeatsForSpecies(20, "mRNA")

        def generate_report(annot1, annot2, type_f, versus):
            versus.write("id_feat\tTAIRStart\tTAIRStop\tTAIRLength\tTAIRExons\tARAStart\tARAStop\tARALength\tARAExons\n")
            ids_B = {d['id_feat'] : d for d in annot2}
            liste_versusTAIR, liste_versusARA = [], []
            for cds in annot1:
                if cds['id_feat'] in ids_B:
                    id_feat = cds['id_feat']
                    liste_versusARA.append(ids_B[id_feat]) # ARAPORT
                    liste_versusTAIR.append(cds) # TAIR12
                    start, stop, length, exons = ids_B[id_feat]['start'],ids_B[id_feat]['stop'], ids_B[id_feat]['length'], ids_B[id_feat]['exons']
                    versus.write(f"{cds['id_feat']}\t{cds['start']}\t{cds['stop']}\t{cds['length']}\t{cds['exons']}\t{start}\t{stop}\t{length}\t{exons}\n")
            fasta_files(f"Versus_TAIR", liste_versusTAIR, type_f, True)
            fasta_files(f"Versus_ARA", liste_versusARA, type_f, True)
        
        with open(f"{dir}/Analysis/TAIR12vsARAPORT11_CDS.txt", "w") as CDSversus:
            generate_report(experts, experts11, "CDS", CDSversus)

        with open(f"{dir}/Analysis/TAIR12vsARAPORT11_mRNA.txt", "w") as mRNAversus:
            generate_report(expertsRNA, experts11RNA, "mRNA", mRNAversus)

        # TAIR12 TEs VS Added
        if TE_analysis:
            TEs = featureDAO.getAllFeatsForSpecies(species_id, "mobile_element")
            TEVerdict = comparing_annotations(not_in_common['Added'], TEs)

            with open(f"{dir}/Analysis/TEsTAIR12.txt", "w") as false_pos_tes :
                for element in TEVerdict['Common']:
                    other1 = element[1]['id_feat']
                    false_pos_tes.write(f"{element[0]['id_feat']}\t{other1}\n")

            # TAIR12 repeat regions VS Added
            repeat_regions = featureDAO.getAllFeatsForSpecies(species_id, "repeat_region")
            repVerdict = comparing_annotations(not_in_common['Added'], repeat_regions)

            with open(f"{dir}/Analysis/repTAIR12.txt", "w") as false_pos_rp :
                for element in repVerdict['Common']:
                    other1 = element[1]['id_feat']
                    false_pos_rp.write(f"{element[0]['id_feat']}\t{other1}\n")    

            completed +=1
            progress_bar(completed/nb_functions*100, 100) 
    
    if categories_analysis: # Analyse de chaque catégorie
        # Analyse protéique
        protein_analysis("Added", not_in_common['Added'], f"{dir}/all_proteins_{name}_Added.fasta", dir, download_db)
        protein_analysis("Missed", not_in_common['Missed'], f"{dir}/all_proteins_{name}_Missed.fasta", dir, download_db)
        protein_analysis("Common", commons, f"{dir}/all_proteins_{name}_CommonHELIXER.fasta", dir, download_db)

        # Analyse PFAM
        pfam_for_expert = comparing_annotations(pfam, not_in_common['Missed'])
        pfam_for_helixer = comparing_annotations(pfam, not_in_common['Added'])
        expert2 = [x[1] for x in not_in_common['Common']]
        pfam_for_both = comparing_annotations(pfam, expert2)

        seuil = 50
        qcover = 60
        save_dir = os.path.abspath(f"{dir}/")
        plus = f"{dir}/Stats/Added_{name}.txt"
        moins = f"{dir}/Stats/Missed_{name}.txt"
        commun = f"{dir}/Stats/Commons_{name}.txt"
        filea = f"{dir}/Protein_Analysis_Common.tsv" # Common
        fileb = f"{dir}/Protein_Analysis_Added.tsv" # Added
        filec = f"{dir}/Protein_Analysis_Missed.tsv" # Missed

        isoform = False
        if not isoform:
            helixer = set([d[1]["id_feat"].rsplit(".", 1)[0] for d in pfam_for_helixer['Common']])
            expert = set([d[1]["id_feat"].rsplit(".", 1)[0] for d in pfam_for_expert['Common']]) 
            both = set([d[1]["id_feat"].rsplit(".", 1)[0] for d in pfam_for_both['Common']]) 
        else:
            helixer = set([d[1]["id_feat"] for d in pfam_for_helixer['Common']])
            expert = set([d[1]["id_feat"] for d in pfam_for_expert['Common']]) 
            both = set([d[1]["id_feat"] for d in pfam_for_both['Common']])    


        input1, input2, input3 = len(helixer), len(expert), len(both)
        os.system(f"Rscript Projet/Rapport/CodesFig/AddedMissedCommon.r {qcover} {seuil} {name2} {save_dir} {plus} {moins} {commun} {input1} {input2} {input3} {filea} {fileb} {filec} {isoform}")
        os.system(f"Rscript Projet/Rapport/CodesFig/ProteinAnalysisSpecies.r {name} {save_dir} {qcover} {fileb} {pfam}")
        
        # Pas mis EggNOG ici 

        completed +=1
        progress_bar(completed/nb_functions*100, 100) 

    if TE_analysis : # Analyse des éléments transposables
        TE_blast(not_in_common['Added'], dir, name2, download_db) 
        bed_to_gff_compare(dir, f"{dir}/Bed/{name}_added.bed", f"{dir}/TES/{name2}.fixed.gff3")
        TE_density(dir, name2, name)
        completed +=1
        progress_bar(completed/nb_functions*100, 100) 

    if protein_comparison: # Comparaison Helixer VS Référence au niveau protéique
        from Bio import SeqIO
        from Bio.Align import PairwiseAligner, substitution_matrices
        import matplotlib.pyplot as plt
        from collections import Counter

        file1 = f"{dir}/all_proteins_{name}_CommonHELIXER.fasta"
        file2 = f"{dir}/all_proteins_{name}_CommonREF.fasta"
        tsv_file = f"{dir}/Stats/Commons_{name}.txt" 

        # Load sequences
        def load_fasta_by_prefix(path):
            lookup = {}
            for r in SeqIO.parse(path, "fasta"):
                prefix = "_".join(r.id.split("_")[:-1])
                lookup[prefix] = str(r.seq).rstrip("*")
            return lookup

        helix_seqs = load_fasta_by_prefix(file1)  
        tair_seqs  = load_fasta_by_prefix(file2)  

        # id -> HELIXER_prefix
        pairs = []
        with open(tsv_file) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                tair_id, helix_prefix = parts[0], parts[1]
                if tair_id in tair_seqs and helix_prefix in helix_seqs:
                    pairs.append((tair_id, tair_seqs[tair_id], helix_prefix, helix_seqs[helix_prefix]))


        # Aligner setup
        blosum62 = substitution_matrices.load("BLOSUM62")
        aligner = PairwiseAligner()
        aligner.substitution_matrix = blosum62
        aligner.mode = "global"
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5

        identities = []

        for tair_id, s_tair, helix_id, s_helix in pairs:
            if s_tair == s_helix:
                identities.append(100.0)
                continue

            aln = aligner.align(s_tair, s_helix)[0]

            matches = sum(
                1 for (s1, e1), (s2, e2) in zip(aln.aligned[0], aln.aligned[1])
                for a, b in zip(s_tair[s1:e1], s_helix[s2:e2]) if a == b
            )
            length = sum(e1 - s1 for s1, e1 in aln.aligned[0])

            identities.append(100 * matches / length if length > 0 else 0)

        # Plot
        def classify(x):
            if x == 100:   return "100%"
            elif x >= 90:  return "90–99%"
            elif x >= 70:  return "70–89%"
            elif x >= 50:  return "50–69%"
            else:          return "<50%"

        categories = [classify(x) for x in identities]
        counts = Counter(categories)
        order = ["<50%", "50–69%", "70–89%", "90–99%", "100%"]
        values = [counts.get(k, 0) for k in order]
        total = len(identities)

        fig, ax = plt.subplots(figsize=(8, 5))
        bars = ax.bar(order, values)

        for bar, val in zip(bars, values):
            pct = 100 * val / total
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(values) * 0.01,
                f"{val}\n({pct:.1f}%)",
                ha="center", va="bottom", fontsize=9
            )
        ax.set_ylim(0, max(values) * 1.12)
        ax.set_xlabel("Protein identity class")
        ax.set_title(
            f"Helixer vs {name2} Protein Identity Distribution\n"
            f"n = {len(pairs):,} (with isoforms)",
            fontsize=12
        )
        ax.set_ylabel("Number of protein pairs")
        plt.tight_layout()
        plt.savefig(f'{dir}/figs/ProtHELIXERVSREF.png')

        completed +=1
        progress_bar(completed/nb_functions*100, 100) 

    if plm_view: # Analyse des PLM chez Helixer Added
        experts_only = True
        noUTR_only = True

        # Ici on cherche à prendre les Helixer chevauchant avec des Experts sans 5'UTR défini
        if not experts_only:
            if noUTR_only :
                noUTR = []
                for duos in not_in_common['Common'] :
                    mRNA = getFeatwithIdfAndType(duos[0]['id_feat'].split(";")[0], "mRNA", species_id)
                    if mRNA is None :
                        continue
                    if duos[0]['complement'] == 1:
                        if duos[0]['start'] - mRNA['start'] == 0 :
                            mRNAHelixer = getFeatwithIdfAndType(duos[1]['id_feat'].split(";")[0], "helixer_mRNA", species_id)
                            if duos[1]['start'] - mRNAHelixer['start'] > 0 :
                                noUTR.append(mRNAHelixer)
                    if duos[0]['complement'] == -1:
                        if mRNA['stop'] - duos[0]['stop'] == 0 :
                            mRNAHelixer = getFeatwithIdfAndType(duos[1]['id_feat'].split(";")[0], "helixer_mRNA", species_id)
                            if mRNAHelixer['stop'] - duos[1]['stop'] > 0 :
                                noUTR.append(mRNAHelixer)
                PLM_analysis(dir, noUTR + not_in_common['Added'], False, john_mail) 

            # Ici on cherche à prendre les Helixer chevauchant avec des Experts avec 5'UTR défini    
            if not noUTR_only :
                withUTR = []
                for duos in not_in_common['Common'] :
                    mRNA = getFeatwithIdfAndType(duos[0]['id_feat'].split(";")[0], "mRNA", species_id)
                    if mRNA is None :
                        continue
                    if duos[0]['complement'] == 1:
                        if duos[0]['start'] - mRNA['start'] > 0 :
                            mRNAHelixer = getFeatwithIdfAndType(duos[1]['id_feat'].split(";")[0], "helixer_mRNA", species_id)
                            withUTR.append(mRNAHelixer)
                    if duos[0]['complement'] == -1:
                        if mRNA['stop'] - duos[0]['stop'] > 0 :
                            mRNAHelixer = getFeatwithIdfAndType(duos[1]['id_feat'].split(";")[0], "helixer_mRNA", species_id)
                            withUTR.append(mRNAHelixer)
                PLM_analysis(dir, withUTR, False, john_mail) 
        
        if experts_only:
            if noUTR_only :
                noUTR = []
                for cds in experts:
                    mRNA = getFeatwithIdfAndType(cds['id_feat'].split(";")[0], "mRNA", species_id)
                    if mRNA is None:
                        continue
                    if cds['complement'] == 1:
                        if cds['start'] - mRNA['start'] == 0 :
                            noUTR.append(mRNA)
                    if cds['complement'] == -1:
                        if mRNA['stop'] - cds['stop'] == 0 :
                            noUTR.append(mRNA)
                PLM_analysis(dir, noUTR, False, john_mail) # Fait pour les experts

            if not noUTR_only :
                withUTR = []
                for cds in experts:
                    mRNA = getFeatwithIdfAndType(cds['id_feat'].split(";")[0], "mRNA", species_id)
                    if mRNA is None :
                        continue
                    if cds['complement'] == 1:
                        if cds['start'] - mRNA['start'] > 0 :
                            withUTR.append(mRNA)
                    if cds['complement'] == -1:
                        if mRNA['stop'] - cds['stop'] > 0 :
                            withUTR.append(mRNA)
                PLM_analysis(dir, withUTR, False, john_mail) 
        #PLM_analysis(dir, not_in_common['Added'], noUTR, True, john_mail, True) # Jaspar

        completed +=1
        progress_bar(completed/nb_functions*100, 100) 

    if UTR_analysis: # Comparaison Helixer VS Référence au niveau UTR
        # Comparaison de toutes les annotations communes
        for_equal = compare_all_feats(dir, 100, "=", not_in_common['Common'])
        for_almost_equal = compare_all_feats(dir, 85, "~", not_in_common['Common'])
        progress_bar(50,100)
        #print(" Comparaison entre annotations communes.\n")

        # Comparaison des régions UTR
        compare_UTR_e = compare_every_UTR(for_equal, species_id) if len(for_equal) != 0 else []
        compare_UTR_ae = compare_every_UTR(for_almost_equal, species_id) if len(for_almost_equal) != 0 else []
        progress_bar(60,100)
        #print(" Comparaison des régions UTR.\n")

        # On passe aux fichiers TSV
        with open(f"{dir}/UTRAnalysis/{name}_UTR_ForSameCDS.tsv", "w") as UTR:
            UTR.write("Chromosome\tGène\tIsoformes\tRéférence\tHelixer\tScore Identité\tType Egalité\t3 UTR Expert Len\t5 UTR Expert Len\t3 UTR Helixer Len\t5 UTR Helixer Len\t5 UTR Helixer\t3 UTR Helixer\t5 UTR Expert\t3 UTR Expert\n")
            if len(compare_UTR_e) != 0 :
                for element in compare_UTR_e:
                    u = f"{element['Chromosome']}\t{element['Gène']}\t{element['Isoformes']}\t{element['Référence']}\t \
                    {element['Helixer']}\t{element['Identité']}\t{element['Type']}\t{element['3Expert Len']}\t{element['5Expert Len']}\t{element['5Helixer Len']}\t{element['3Helixer Len']}\t{element['5 UTR Helixer']}\t{element['3 UTR Helixer']}\t{element['5 UTR Expert']}\t{element['3 UTR Expert']}\n"
                    UTR.write(u)
        progress_bar(90,100)
        #print(" Intégration aux fichiers TSV.\n")

        with open(f"{dir}/UTRAnalysis/{name}_UTR_ForAlmostSameCDS.tsv", "w") as UTR1:
            UTR1.write("Chromosome\tGène\tIsoformes\tRéférence\tHelixer\tScore Identité\tType Egalité\t3 UTR Expert Len\t5 UTR Expert Len\t3 UTR Helixer Len\t5 UTR Helixer Len\t5 UTR Helixer\t3 UTR Helixer\t5 UTR Expert\t3 UTR Expert\n")
            if len(compare_UTR_ae) != 0 :
                for element in compare_UTR_ae:
                    u1 = f"{element['Chromosome']}\t{element['Gène']}\t{element['Isoformes']}\t{element['Référence']}\t \
                    {element['Helixer']}\t{element['Identité']}\t{element['Type']}\t{element['3Expert Len']}\t{element['5Expert Len']}\t{element['5Helixer Len']}\t{element['3Helixer Len']}\t{element['5 UTR Helixer']}\t{element['3 UTR Helixer']}\t{element['5 UTR Expert']}\t{element['3 UTR Expert']}\n"
                    UTR1.write(u1)
        progress_bar(100,100)
        #print(" It's so over ahh")

        completed +=1
        progress_bar(completed/nb_functions*100, 100) 

    if overall_RNA: # Comparaison taille et nombre exons dans mRNA
        with open(f"{dir}/Stats/ExpertsLENEXONS.txt", "w") as ex, open(f"{dir}/Stats/HelixerLENEXONS.txt", "w") as he:
            ex.write("id_feat\tLength\tNbExons\tUTR5\tUTR3\n")
            he.write("id_feat\tLength\tNbExons\tUTR5\tUTR3\n")
            for exp in experts:
                if len(exp['id_feat'].split(".")) > 1 and exp['id_feat'].split(".")[1] != "1":
                    continue
                expert_mRNA = getFeatwithIdfAndType(exp['id_feat'], "mRNA", species_id)
                if expert_mRNA is None:
                    continue
                if exp['complement'] == 1:
                    expert5 = sequenceDAO.getSequenceDNAFromTo(expert_mRNA['id_seq'], expert_mRNA['start'], exp['start'])
                    expert3 = sequenceDAO.getSequenceDNAFromTo(expert_mRNA['id_seq'], exp['stop'], expert_mRNA['stop'])
                else:
                    expert3 = sequenceDAO.getSequenceDNAFromTo(expert_mRNA['id_seq'], expert_mRNA['start'], exp['start'])
                    expert5 = sequenceDAO.getSequenceDNAFromTo(expert_mRNA['id_seq'], exp['stop'], expert_mRNA['stop'])
                ex.write(f"{exp['id_feat']}\t{exp['length']}\t{exp['exons']}\t{len(expert5)}\t{len(expert3)}\n")


            for hel in helixers:
                helixer_mRNA = getFeatwithIdfAndType(hel['id_feat'], "helixer_mRNA", species_id)
                if helixer_mRNA is None :
                    continue
                if hel['complement'] == 1:
                    helixer5 = sequenceDAO.getSequenceDNAFromTo(helixer_mRNA['id_seq'], helixer_mRNA['start'], hel['start'])
                    helixer3 = sequenceDAO.getSequenceDNAFromTo(helixer_mRNA['id_seq'], hel['stop'], helixer_mRNA['stop'])
                else:
                    helixer3 = sequenceDAO.getSequenceDNAFromTo(helixer_mRNA['id_seq'], helixer_mRNA['start'], hel['start'])
                    helixer5 = sequenceDAO.getSequenceDNAFromTo(helixer_mRNA['id_seq'], hel['stop'], helixer_mRNA['stop'])
                he.write(f"{hel['id_feat']}\t{hel['length']}\t{hel['exons']}\t{len(helixer5)}\t{len(helixer3)}\n")

        completed +=1
        progress_bar(completed/nb_functions*100, 100)

    if busco: # Analyse BUSCO pour chaque protéome
        file0 = f"{dir}/all_proteins_{name}_CommonREF.fasta"
        file1 = f"{dir}/all_proteins_{name}_CommonHELIXER.fasta"
        file2 = f"{dir}/all_proteins_{name}_Added.fasta"
        file3 = f"{dir}/all_proteins_{name}_Missed.fasta"


        os.system(f"busco -i {file0} \
            -o {dir}/Analysis/busco_ref_common \
            -l embryophyta_odb10 \
            -m protein -c 8 > /dev/null 2>&1")


        os.system(f"busco -i {file1} \
            -o {dir}/Analysis/busco_helixer_common \
            -l embryophyta_odb10 \
            -m protein -c 8 > /dev/null 2>&1")

        os.system(f"busco -i {file2} \
            -o {dir}/Analysis/busco_helixer_added \
            -l embryophyta_odb10 \
            -m protein -c 8 > /dev/null 2>&1")

        os.system(f"busco -i {file3} \
            -o {dir}/Analysis/busco_ref_specific \
            -l embryophyta_odb10 \
            -m protein -c 8 > /dev/null 2>&1")
        
        completed +=1
        progress_bar(completed/nb_functions*100, 100)


if __name__ == '__main__':
    """
        id name
        34 Arabidopsis thaliana TAIR12
        23 Zea mays B73
        30 Arabidopsis lyrata
        20 Arabidopsis thaliana
        22 Brachypodium distachyon
        26 Malus domestica
        28 Rosa chinensis Old Blush
        24 Solanum lycopersicum
    """
    
    
    # Set parameters
    parser = argparse.ArgumentParser(description='Search differences between and Helixer and Reference annotations')
    parser.add_argument('--dir', type=str, dest='dir',  required=True, default='Data', help='Destination directory')
    parser.add_argument('--gender', type=str,  required=True, dest='gender', help='Schéma à attaquer')
    parser.add_argument('--speId', type=int,  required=True, dest='speId', default=-1, help='Source species Id')
    parser.add_argument('--doAll', type=bool,  required=False, dest='doAll', default=False, help='Whole analysis')
    parser.add_argument('--UTRAnalysis', type=bool,  required=False, dest='UTRAnalysis', default=False, help='UTR analysis within Helixer')
    parser.add_argument('--RNAAnalysis', type=bool,  required=False, dest='RNAAnalysis', default=False, help='Length and Exons Helixer VS Reference')
    parser.add_argument('--CatAnalysis', type=bool,  required=False, dest='CatAnalysis', default=False, help='Added/Missed/Commons functionnal analysis')
    parser.add_argument('--TEAnalysis', type=bool,  required=False, dest='TEAnalysis', default=False, help='Analysis of TEs predicted by Helixer')
    parser.add_argument('--PLMView', type=bool,  required=False, dest='PLMView', default=False, help='Detection of motifs in 5UTR of Helixer')
    parser.add_argument('--Email', type=str,  required=False, dest='Email', default="67@gmail.com", help='Email address required for PLMView')
    parser.add_argument('--makeDB', type=bool,  required=False, dest='makeDB', default=False, help='Creates a BlastP DB')
    parser.add_argument('--proteinCompare', type=bool,  required=False, dest='proteinCompare', default=False, help='Protein Helixer VS Protein Reference')
    parser.add_argument('--BUSCO', type=bool,  required=False, dest='BUSCO', default=False, help='Run BUSCO for Added/Missed/Common')
    parser.add_argument('--generate', type=bool,  required=False, dest='generate', default=False, help='Generate protein .fasta files for each category')

    # Read parameters
    args = parser.parse_args()

    # Set Schema
    FLAGdb.setGender(args.gender)

    id = args.speId
    dir = args.dir

    all = args.doAll
    utr = args.UTRAnalysis
    rna = args.RNAAnalysis
    cat = args.CatAnalysis
    te = args.TEAnalysis
    plmview = args.PLMView
    plmmail = args.Email
    dl = args.makeDB
    protein = args.proteinCompare
    busco = args.BUSCO
    generate = args.generate

    exporting_results(
        id,
        dir,
        all,
        utr,
        rna,
        cat,
        te,
        plmview,
        plmmail,
        dl,
        protein,
        busco,
        generate
    )
    