import argparse
import os
import sys

# Ajout du chemin de recherche des librairies (Elle se trouve dans un repertoire au dessus)
p = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, p)

from FLAGdbServer import setting                                    # module d'initialisation des variables
from FLAGdbServer.tools import logger                               # Module d'initialisation des logs
from FLAGdbServer.database import FLAGdb                            # Module de connection a la base
from FLAGdbServer.apis.species import speciesDAO, SpeciesNotFound   # Module pour interroger des especes
from FLAGdbServer.apis.feature import featureDAO                    # Module pour interroger des features
from FLAGdbServer.apis.sequence import sequenceDAO                  # Module pour interroger des séquences
from ToolsServer.extractSequence import calcProtein


if __name__ == '__main__':
    # Set parameters
    parser = argparse.ArgumentParser(description='Search helixer difference between reference annotation')
    parser.add_argument('--dir', type=str, dest='dir',  required=True, default='Data', help='destination directory')
    parser.add_argument('--gender', type=str,  required=True, dest='gender', help='schema a attaquer')
    parser.add_argument('--speId', type=int,  required=True, dest='speId', default=-1, help='source species Id')

    args = parser.parse_args()

    dir = args.dir
    gender = args.gender
    speId = args.speId

    FLAGdb.setGender(args.gender)

    cptHelix = 0
    cptPFAM = 0
    for file in os.listdir(dir) : 
        nom, ext = os.path.splitext(file)
        if ext != '.fasta' :
            continue

        (gene, id) = nom.rsplit('_',1)
        cptHelix += 1
        helixer_CDS = featureDAO.getFeatureForId(id)
        pfams = featureDAO.getFeatureForSequence( helixer_CDS['sid'], "PFAM",  helixer_CDS['start']-1,  helixer_CDS['stop']+1)


        print( cptHelix, "Gene :", gene, " id : ",id ," on sequence : ", helixer_CDS['sid'], " position : ", helixer_CDS['start'], "..", helixer_CDS['stop'], " pfams :", len(pfams))
        if len(pfams) != 0 :
            cptPFAM +=1

    print("Helixer total sans annotaiton : ",cptHelix, " avec PFAM :", cptPFAM)