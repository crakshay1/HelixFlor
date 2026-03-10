# Inclure des librairies GBOT
import os
import sys

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
from FLAGdbServer.apis.sequence import sequenceDAO                  # Module pour interroger des séquences


#################################################################################
def are_overlapping(start1: int, stop1: int, start2: int, stop2: int, seuil: int) -> bool:
    """
        Prend en paramètre les starts et stops de deux éléments.
        Renvoie True si les deux éléments se chevauchent, sinon False.
    """
    return min(stop1, stop2) - max(start1, start2) > seuil # Si c'est négatif -> Pas de chevauchement


def find_missingAnnotationV1(id : int, feature : str) -> dict :
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
        print(sequence['accession'])

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


def get_score(start1: int, stop1: int, start2: int, stop2: int) -> float:
    """
        Prend en paramètre les coordonnées des exons, et renvoie un score d'identité
    """
    ratio1 = stop1 - start1 
    ratio2 = stop2 - start2 
    return min(ratio1,ratio2) / max(ratio1,ratio2)


def compare_both_feats(feat1: dict, feat2: dict) -> dict:
    """
        Prend en paramètre deux features.
        Renvoie un rapport de type dict pour résumer la comparaison des deux features.
    """
    if feat1['complement'] == feat2['complement'] : # On s'assure qu'on est sur le même sens
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
                    calculs.append(get_score(ref[0][i], ref[1][i],helixer[0][j], helixer[1][j])) # Normalement ça donne un résultat entre 0 et 1 pour chaque exon
        
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
        
        # Les types de différencences entre les annotations
        if score == 100:
            if ref[0][0] != helixer[0][0] or ref[1][-1] != helixer[1][-1]:
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
        return {'Gène Ref': feat1['id_feat'], 'Helixer': feat2['id_feat'], 'Isoforme Ref': "" , 'Isoformes': [] ,'Identité': score, 'Type': type_c}
    

def compare_all_feats(id : int, feat_type : str, seuil : int, type: str) -> list[dict]:
    """
        Prend en paramètre l'id de l'espèce et le type de features qu'on veut comparer.
        Prend aussi un compte un % seuil et un type de superposition.
        Renvoie une liste de dictionnaire dans lequel on retrouve :
        le feature de référence, le feature de comparaison, et le score d'identité, le type d'égalité.
        La liste est filtrée selon le seuil et le type mis en paramètre.
        Plus le score est proche de 100, plus les features sont similaires.
    """
    seuil = max(min(seuil, 100), 0) # Pas en dessous de 0 et pas au dessus de 100
    type = type if type in "=~ckmjn" else "="  # Un type par défaut si le type est mal précisé

    overlaps = obtain_both_feats(id, feat_type) # Je prends tous les features qui se chevauchent
    rapport = [] # Rapport non filtré
    rapport1 = [] # Rapport filtré

    # Etablissement d'un rapport sans filtre
    for comparaison in overlaps :
        resultat = compare_both_feats(comparaison[0], comparaison[1])
        if resultat: # Peut renvoyer None si les deux feats n'ont pas le même sens !
            rapport1.append(resultat)

    # Filtrage
    liste_opps = set() # Si c'est dedans, on en veut pas pour plus tard
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
        result['Isoforme Ref'] = result['Gène Ref'] # On veut juste le nom du gène
        result['Gène Ref'] = result['Gène Ref'] if not "." in result['Gène Ref'] else result['Gène Ref'][:-2]
        result['Isoformes'] = [iso['Gène Ref'] for iso in isoformes] # Une liste des transcrits alternatifs
        #print(result)

        # Liste des features selon les filtres
        if type != "=" :
            if result['Identité'] >= seuil and result['Type'] == type:
                rapport.append(result)
        else :
            if result['Type'] == type: # Quand on est dans le cas "=", tous les scores sont à 100
                rapport.append(result)
    return rapport
     

def comparing_UTR(feat1: str, feat2: str) -> dict:
    """
        Prend en paramètre deux features et renvoie les différences dans les régions 5' UTR et 3'UTR
    """
    expert_mRNA = featureDAO.getFeatByIdfAndType(feat1, "mRNA")
    helixer_mRNA = featureDAO.getFeatByIdfAndType(feat2, "helixer_mRNA")
    expert_CDS= featureDAO.getFeatByIdfAndType(feat1, "CDS")
    helixer_CDS = featureDAO.getFeatByIdfAndType(feat2, "helixer_CDS")
    
    comparaison = compare_both_feats(expert_mRNA, helixer_mRNA)
    expert5 = sequenceDAO.getSequenceDNAFromTo(expert_mRNA['id_seq'], expert_mRNA['start'], expert_CDS['start'])
    expert3 = sequenceDAO.getSequenceDNAFromTo(expert_mRNA['id_seq'], expert_CDS['stop'], expert_mRNA['stop'] )
    helixer5 = sequenceDAO.getSequenceDNAFromTo(helixer_mRNA['id_seq'], helixer_mRNA['start'], helixer_CDS['start'])
    helixer3 = sequenceDAO.getSequenceDNAFromTo(helixer_mRNA['id_seq'], helixer_CDS['stop'], helixer_mRNA['stop'] )
    utr5 =""
    utr3 = ""

    return {'Référence': expert_mRNA['id_feat'], 'Helixer': helixer_mRNA['id_feat'],'Identité': comparaison['Identité'], 'Type': comparaison['Type'], '5 UTR': utr5, '3 UTR': utr3}

########################## TEST #########################################
id = 20
espece = speciesDAO.getSpeciesById(id)
liste_seq = sequenceDAO.getSequencesBySpeciesId(id)
liste_CDS = featureDAO.getAllFeatsForSequence(liste_seq[0]['id'], "CDS") # Ouuuf j'accède à tous les CDS
liste_CDS_Helixer = featureDAO.getAllFeatsForSequence(liste_seq[1]['id'], "helixer_CDS") # Pareil pour tous les Helixer CDS
liste_missed, liste_added = [], []

"""
test = find_missingAnnotationV1(id, "CDS")
added = test['Added']
missed = test['Missed']
cds = test['Expert']
helixer = test['Helixer']
espece = speciesDAO.getSpeciesById(id)
#print(f"Pour le génome de {espece['name']}, \n On a {cds} CDS prédits par les experts dont {len(missed)} ratés par Helixer.\n On a {helixer} CDS prédits par Helixer dont {len(added)} CDS ajoutés par rapport à l'annotation experte.")


in_both = obtain_both_feats(id, "CDS")
#print(f"Pour le génome de {espece['name']}, \n On a {len(in_both)} CDS prédits par les experts et par Helixer.")


stat = compare_all_feats(id, "CDS", 1, "=")
print(stat)
"""
print(liste_CDS_Helixer)
print("\n")
print(featureDAO.getFeatureForId(628))
print("\n")
print(featureDAO.getFeatByIdfAndType("Arabidopsis_thaliana_ChrC_000020", "helixer_mRNA"))
print(sequenceDAO.getSequenceDNAFromTo(7, 1, 300))
