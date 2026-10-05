from flask_restx import Namespace, Resource, fields
from flask import render_template, request, redirect
from flask import make_response
from flask_restx import reqparse

from FLAGdbServer.database import FLAGdb
from FLAGdbServer.database import QueryMaster
from FLAGdbServer.tools import logger
from FLAGdbServer.apis.sequence import sequenceDAO

#from setting import *
from ..setting import homeURL
from ..setting import URL
from ..setting import URLPATH
from ..setting import SUBURIPATH
import datetime

# TODO modifier cela pour que flask ne retourne que ce qui est necessaire en fonction du type de l'objet
# exp paralogues necessaire que pour CDS!!! 

feature_ns = Namespace('Feature', description='Feature operations')
genecard_ns = Namespace('GeneCard',description='Feature operations')

feature = feature_ns.model('Feature', {
    'id': fields.Integer(readOnly=True, description='Id for the Feature'),
    'idf': fields.String(required=True,
                             description='id_feat for the Feature'),
    'start': fields.Integer(required=True,
                            description='Start for the Feature'),
    'type': fields.String( description='Type for the Feature'),
    'stop': fields.Integer(required=True, description='Stop for the Feature'),
    'sens': fields.Integer(required=True, description='Sens for the Feature'),
 #   'color': fields.Integer(required=True, description='Feature Color '),
    'te_type': fields.List(fields.String(required=False, description='')),

    'p': fields.List(
                  fields.List(fields.Integer,
                              required=True,
                              description='Positions for the Feature')),
    'alternatif': fields.List(fields.Integer,
                                required=True,
                                description='number of alternatif gene'),
    # masterqual (mq) is type dependant for _features
    # CDS => product.
    # FST => class.
    'comments': fields.String(required=False, description='Comments for the Feature'),
    'sid': fields.Integer(required=True, description='Sequence Id'),
    'mq': fields.List(fields.String(required=False, description='')),
    'repeat': fields.List(fields.String(required=False, description='is repeat')),
#    'has_paralogues': fields.Boolean(required=False, description='is there some paralogs'),
#    'has_orthologues': fields.Boolean(required=False, description='is there some orthologs'),
    'color': fields.Integer(
                                description='Color'),
    })




parser = reqparse.RequestParser()
parser.add_argument('sequenceId', type=int)
parser.add_argument('type', type=str)
parser.add_argument('userid', type=int)
parser.add_argument('project', type=str)
parser.add_argument('start', type=int)
parser.add_argument('stop', type=int)

parser2 = reqparse.RequestParser()
parser2.add_argument('idList', type=str)

def chunkSequence(seq):
        n = 80
        # remove space bug pfam;
        seq = seq.replace(' ', '')
        chunks = [seq[i:i+n] for i in range(0, len(seq), n)]
        return chunks

class FLAGDBFeatureDAO:
    def __init__(self):
        self.list = None

    def addFeatureToProject(self, userid, project, sequenceid,\
            complement, comments, typefeat, starts, stops, idfeat):
        """ Function to add feature to edit projects 
        """

        location = ''
        if complement == -1 : 
            location='complement('
        location+="join("
        i=0
        for start in starts:
            location+=str(start)+".."+str(stops[i])+","
            i+=1
        location = location[:-1]
        location += ')'
        if complement == -1:
            location+=")"
        # Features
        query = 'insertFeatures'
        featureQuery = QueryMaster[query].replace('$seqid$', str(sequenceid) )\
                            .replace('$type$', typefeat )\
                            .replace('$loc$', location )\
                            .replace('$sens$', str(complement) )\
                            .replace('$start$', str(starts[0]) )\
                            .replace('$stop$', str(stops[-1]) )\
                            .replace('$date$', str(datetime.date.today()) )\
                            .replace('$idf$', str(idfeat) )\
                            .replace('$dnaid$', str('') )\
                            .replace('$comm$', comments)\
                            .replace('$souid$', str(-1))

        if userid==-10:
            schema = '\"'+project+'\"'
            response = FLAGdb.query('SELECT ID FROM "'+project+'".sources where source=\'USER\'')
            sourceId = response.fetchone()[0]
            featureQuery = featureQuery.replace('-1)', str(sourceId)+')')


        else : 
            schema = 'user_'+project
            featureQuery = featureQuery.replace('source_id)', 'source_id, user_id)')
            featureQuery = featureQuery.replace('-1)', '-1,'+str(userid)+')')

        featureQuery = featureQuery.replace('$schema$', schema+'.')
        
        logger.debug(featureQuery)
        response = FLAGdb.query(featureQuery)
        featureId = response.fetchone()[0]
        # locations
        i=0
        for start in starts:
            query = 'insertLocations'
            locationQuery = QueryMaster[query].replace('$fid$', str(featureId) )\
                    .replace('$sens$', str(complement) )\
                    .replace('$start$', str(start) )\
                    .replace('$stop$', str(stops[i]) )
            i+=1
            locationQuery = locationQuery.replace('$schema$', schema+'.')

            response = FLAGdb.query(locationQuery)


    def removeFeatureToProject(self, userid, project, featureid) :
        # eliminer les locations
        FLAGdb.query("delete from user_"+project+".locations where feature_id="+str(featureid))
        # eliminer les qualifiers
        FLAGdb.query("delete from user_"+project+".qualifiers where feature_id="+str(featureid))
        # eliminer le feature
        FLAGdb.query("delete from user_"+project+".features where id="+str(featureid))


    def getFeatureListByIds(self, f_ids):
        """ Fonction qui recupere la liste des Features a partir d'une liste d'IDs.

        """
        logger.debug("getFeatureListByIds")
        
         
        teforfly = ''
        if FLAGdb._dbSchema=='droso.':
            teforfly="""
             CASE 
                WHEN f.type = 'TE' then ( select ARRAY[type, subtype] from $schema$et_description where name=id_feat)
             END as te_type,
             """
        # On retire les [] du str
        query = QueryMaster['getFeatureListByIds']\
            .replace('$teforfly$', teforfly)\
            .replace('$ids$', str(f_ids)[1:-1])

        cursor = FLAGdb.query(query)

        featList = []
        row = cursor.fetchone()
        while row is not None:
            featList.append(row[0])
            row = cursor.fetchone()

        if len(featList) == 0:
            #    raise NoResultFound("prout");
            feature_ns.abort(404, "Not Found")

        logger.debug(featList)
        FLAGdb.queryclose()
        return featList

    def getFeatureListByIdFeats(self, f_ids):
        logger.debug("getFeatureListByIdFeats")
        

        query = QueryMaster['getFeatureListByIdFeats']\
            .replace('$ids$', str(f_ids)[1:-1])  # On retire les [] du str

        cursor = FLAGdb.query(query)

        featList = []
        row = cursor.fetchone()
        while row is not None:
            featList.append(row[0])
            row = cursor.fetchone()

        if len(featList) == 0:
            #    raise NoResultFound("prout");
            feature_ns.abort(404, "Not Found")

        logger.debug(featList)
        FLAGdb.queryclose()
        return featList
    
    def getCDSListByIdFeats(self, f_ids):
        logger.debug("getCDSListByIdFeats")
        

        query = QueryMaster['getCDSListByIdFeats']\
            .replace('$ids$', str(f_ids)[1:-1])  # On retire les [] du str
        cursor = FLAGdb.query(query)

        featList = []
        row = cursor.fetchone()
        while row is not None:
            featList.append(row[0])
            row = cursor.fetchone()

        if len(featList) == 0:
            #    raise NoResultFound("prout");
            feature_ns.abort(404, "Not Found")

        logger.debug(featList)
        FLAGdb.queryclose()
        return featList


    def getCDSListById(self, f_ids):
        logger.debug("getCDSListById")
        

        query = QueryMaster['getCDSListById']\
            .replace('$ids$', str(f_ids)[1:-1])  # On retire les [] du str
        cursor = FLAGdb.query(query)

        featList = []
        row = cursor.fetchone()
        while row is not None: 
            featList.append(row[0])
            row = cursor.fetchone()

        if len(featList) == 0:
            #    raise NoResultFound("prout");
            feature_ns.abort(404, "Not Found")

        logger.debug(featList)
        FLAGdb.queryclose()
        return featList

    def getInteractById(self, f_id, depth=0):
        logger.debug("getInteractForId")
        

        query = QueryMaster['getInteractForId']\
            .replace('$id$', str(f_id)) \
            .replace('$depth$', str(depth)) 
        cursor = FLAGdb.query(query)

        
        row = cursor.fetchone()
       
        return row
    
    def getisInteractId(self, f_id):
        logger.debug("isInteract")
        

        query = QueryMaster['isInteract']\
            .replace('$id$', str(f_id))
        cursor = FLAGdb.query(query)

        
        row = cursor.fetchone()
       
        return row


    def getAllFeatsForSequence(self, id, type ):
        logger.debug("getAllFeatsForSequence")

        query = QueryMaster['getAllFeatsForSequence'].replace('$id$', str(id)).replace('$type$', str(type))
        cursor = FLAGdb.query(query)
        self.list = cursor.fetchone()[0]
        if self.list is None:
            feature_ns.abort(404, "id {} doesn't exist".format(id))
        logger.debug(self.list)
        FLAGdb.queryclose()
        return self.list
    
    def getAllFeatsForSpecies(self, id, type ):
        logger.debug("getAllFeatsForSpecies")

        query = QueryMaster['getAllFeatsForSpecies'].replace('$id$', str(id)).replace('$type$', str(type))
        cursor = FLAGdb.query(query)
        self.list = cursor.fetchone()[0]
        if self.list is None:
            feature_ns.abort(404, "id {} doesn't exist".format(id))
        logger.debug(self.list)
        FLAGdb.queryclose()
        return self.list
    
    def getFeatByIdfAndType(self, idf, type, species_id) :
        logger.debug("getFeatByIdfAndType")

        query = QueryMaster['getFeatByIdfAndType'].replace('$idf$', str(idf)).replace('$type$', str(type)).replace('$sp$', str(species_id))
        cursor = FLAGdb.query(query)
        row = cursor.fetchone()

        if row is None:
            #    raise NoResultFound("prout"); # prout
            feature_ns.abort(404, "Not Found")

        list = row[0]
        logger.debug(list)
        FLAGdb.queryclose()
        return list

    
    # def getList(self):
    #     logger.debug("getSequences")
    #     query = QueryMaster['getSequences']
    #     logger.debug(query)
    #     cursor = FLAGdb.query(query)
    #     self.list = cursor.fetchone()[0]
    #     logger.debug(self.list);
    #     FLAGdb.queryclose()
    #     return self.list

    def getProductForFeature(self, id,):
        logger.debug("getProductForFeature")
        

        query = QueryMaster['getProductForFeature'].replace('$id$', str(id))
        cursor = FLAGdb.query(query)
        row = cursor.fetchone()

        if row is None:
            #    raise NoResultFound("prout");
            logger.debug(row)
            return None

        product = row[0]
        FLAGdb.queryclose()
        return chunkSequence(product)

    def getLocationForFeature(self, id,):
        logger.debug("getLocationForFeature")
        

        query = QueryMaster['getLocationForFeature'].replace('$id$', str(id))
        cursor = FLAGdb.query(query)
        row = cursor.fetchone()

        if row is None:
            #    raise NoResultFound("prout");
            logger.debug(row)
            return None

        location = row[0]
        FLAGdb.queryclose()
        return location

    
    def getUserFeatureForSequence(self, id, userid, project, start, stop):
        logger.debug("getUserFeaturesForSequence")

        query = QueryMaster['getUserFeaturesForSequence']\
            .replace('$id$', str(id))\
            .replace('$uid$', str(userid))\
            .replace('$start$', str(start))\
            .replace('$stop$',  str(stop))\
            .replace('$schema$', "user_"+project+".")

        logger.debug(query)
        list = []
        #
        cursor = FLAGdb.query(query)
        row = cursor.fetchone()
        list = row[0]
        logger.debug(list)
        FLAGdb.queryclose()
        if list is None:
            list = []
        return list
       
      
    def getFeatureForSequence(self, id, type, start, stop):
        logger.debug("getFeaturesForSequence")
        
        #cursor = FLAGdb.query(
        #    "select species_id from $schema$sequences where id="+str(id))

        #row = cursor.fetchone()

        teforfly = ''

        if FLAGdb._dbSchema=='droso.':
            teforfly="""
             CASE 
                WHEN '$type$' = 'TE' then ( select ARRAY[type, subtype] from $schema$et_description where name=id_feat)
             END as te_type,
             """

        query = QueryMaster['getFeaturesForSequence']\
            .replace('$teforfly$', teforfly)\
            .replace('$id$', str(id))\
            .replace('$type$', type)\
            .replace('$start$', str(start))\
            .replace('$stop$',  str(stop))
                    #    .replace('$speid$', str(row[0]))\

        cursor = FLAGdb.query(query)
        row = cursor.fetchone()

        if row is None:
            #    raise NoResultFound("prout");
            feature_ns.abort(404, "Not Found")

        self.list = row[0]
        logger.debug(self.list)
        FLAGdb.queryclose()
        if self.list is None:
            self.list = []
        return self.list

    def getFeatureForId(self, idx):
        logger.debug("getFeatureById")
        

        teforfly=""
        if FLAGdb._dbSchema=='droso.':
            teforfly="""
             CASE 
                WHEN f.type = 'TE' then ( select ARRAY[  name, type, subtype, file, refid, classification  ] from $schema$et_description where name=id_feat)
             END as te_type,
             """

        query = QueryMaster['getFeatureById'].replace('$teforfly$', teforfly)\
            .replace('$id$', str(idx))

        cursor = FLAGdb.query(query)
        row = cursor.fetchone()

        if row is None:
            #    raise NoResultFound("prout");
            feature_ns.abort(404, "Not Found")

        list = row[0]
        logger.debug(list)
        FLAGdb.queryclose()
        return list

    def getOrthologuous(self, idfeat):  
        logger.debug("getOrthologuous")
       

        query = QueryMaster['getOrthologuesForId']\
            .replace('$id$', str(idfeat))
           
        cursor = FLAGdb.query(query)

        row = cursor.fetchone()
        

        logger.debug(row)
        FLAGdb.queryclose()
        return ( row)

    def getOrthologuousAndLinks(self, idfeat):
        logger.debug("getOrthologuous")
       

        query = QueryMaster['getOrthologuesId']\
            .replace('$id$', str(idfeat))
           
        cursor = FLAGdb.query(query)

        row = cursor.fetchone()
        

        logger.debug(row)
        FLAGdb.queryclose()
        return ( row)


    def hasOrthologuous(self, idfeat):
        logger.debug("hasOrthologuous")
       

        query = QueryMaster['hasOrthologuesForId']\
            .replace('$id$', str(idfeat))
           
        cursor = FLAGdb.query(query)

        row = cursor.fetchone()

        logger.debug(row)
        FLAGdb.queryclose()
        if row[0] > 0 :
            return True
        else:
            return False

    def hasHomologuous(self, idfeat):
        logger.debug("hasHomologuous")
       

        query = QueryMaster['hasHomologuesForId']\
            .replace('$id$', str(idfeat))
           
        cursor = FLAGdb.query(query)

        row = cursor.fetchone()

        logger.debug(row)
        FLAGdb.queryclose()
        return row[0]


    def hasParaloguous(self, idfeat):
        logger.debug("hasParaloguous")
       

        query = QueryMaster['hasParaloguesForId']\
            .replace('$id$', str(idfeat))
           
        cursor = FLAGdb.query(query)

        row = cursor.fetchone()

        logger.debug(row)
        FLAGdb.queryclose()
        if row[0] > 0 :
            return True
        else:
            return False




    def getParaloguous(self, idfeat):
        logger.debug("getParaloguous")
       

        query = QueryMaster['getParaloguesForId']\
            .replace('$id$', str(idfeat))
           
        cursor = FLAGdb.query(query)

        row = cursor.fetchone()
    
        logger.debug(row)
        FLAGdb.queryclose()
        return ( row)


    def getFeaturesIdWithTypeIdFeat(self, type, idfeat):
        logger.debug("getFeaturesIdWithTypeIdFeat")
       

        query = QueryMaster['getFeaturesIdWithTypeIdFeat']\
            .replace('$id$', str(idfeat))\
            .replace('$t$', type)
        cursor = FLAGdb.query(query)

        rows = []
        row = cursor.fetchone()
        while row is not None:
            rows.append(row[0])
            row = cursor.fetchone()

        if len(rows) == 0:
            #    raise NoResultFound("prout");
            feature_ns.abort(404, make_response(render_template('404.html')))

        logger.debug(rows)
        FLAGdb.queryclose()
        return (cursor.rowcount, rows)

    def getFeatureWithTypeID(self, type, id):
        logger.debug("getFeatureWithTypeID")
       

        query = QueryMaster['getFeatureWithTypeID']\
            .replace('$id$', str(id))\
            .replace('$type$', str(type))

        cursor = FLAGdb.query(query)

        row = cursor.fetchone()

        if row is None:
            #    raise NoResultFound("prout");
            feature_ns.abort(404, "Not Found")

        logger.debug(row[0])
        FLAGdb.queryclose()
        return row[0]



    def getFeaturesIdWithIdFeat(self, idfeat):
        """
            recuperer la liste des features pour un idfeat donné.
            Ex pour une PFAM on aura plusieurs features car on a une
            multilocalisation pour un CDS on aura 2 features le mRNA et le CDS
        """
        logger.debug("getFeaturesIdWithIdFeat")
       

        query = QueryMaster['getFeaturesIdWithIdFeat']\
            .replace('$id$', str(idfeat))

        cursor = FLAGdb.query(query)

        rows = []
        row = cursor.fetchone()
        while row is not None:
            rows.append(row[0])
            row = cursor.fetchone()

        if len(rows) == 0:
            #    raise NoResultFound("prout");
            feature_ns.abort(404,  idfeat + ' not found in GBOT')

        logger.debug(rows)
        FLAGdb.queryclose()
        return (cursor.rowcount, rows)

        


featureDAO = FLAGDBFeatureDAO()

@feature_ns.response(404, 'Feature not found')
@feature_ns.route('/OrthologuousAndLinks/<int:featureid>')
class getOrtologuousAndLinks(Resource):
    ''' Get the list of Features '''
    @feature_ns.doc('list of Features')
    def get(self, featureid):
        ''' List all Features '''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        jsonlist = featureDAO.getOrthologuousAndLinks(featureid)
       
        return jsonlist

@feature_ns.response(404, 'Feature not found')
@feature_ns.route('/hasOrthologuous/<int:featureid>')
class hasOrthologuousForId(Resource):
    ''' Get the list of Features '''
    @feature_ns.doc('list of Features')
    def get(self, featureid):
        ''' List all Features '''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        boolean  = featureDAO.hasOrthologuous(featureid)
       
        return boolean


@feature_ns.response(404, 'Feature not found')
@feature_ns.route('/hasParaloguous/<int:featureid>')
class hasParaloguousForId(Resource):
    ''' Get the list of Features '''
    @feature_ns.doc('list of Features')
    def get(self, featureid):
        ''' List all Features '''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        boolean  = featureDAO.hasParaloguous(featureid)
       
        return boolean

@feature_ns.response(404, 'Feature not found')
@feature_ns.route('/hasHomologuesForId/<int:featureid>')
class hasHomologuesForId(Resource):
    ''' Get the list of Features '''
    @feature_ns.doc('list of Features')
    def get(self, featureid):
        ''' List all Features '''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        return featureDAO.hasHomologuous(featureid)
       


@feature_ns.response(404, 'Feature not found')
@feature_ns.route('/Orthologuous/<id_feat>')
class getOrthologuous(Resource):
    ''' Get the list of Features '''
    @feature_ns.doc('list of Features')
    def get(self, id_feat):
        ''' List all Features '''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        print(id_feat)
        list = featureDAO.getCDSListByIdFeats("['"+id_feat+"%']") 

        print(list)

        jsonlist = featureDAO.getOrthologuous(list[0]['id'])
       
        return jsonlist



@feature_ns.response(404, 'Feature not found')
@feature_ns.route('/Orthologuous/<int:featureid>')
class getOrthologuous(Resource):
    ''' Get the list of Features '''
    @feature_ns.doc('list of Features')
    def get(self, featureid):
        ''' List all Features '''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        jsonlist = featureDAO.getOrthologuous(featureid)
       
        return jsonlist


@feature_ns.response(404, 'Feature not found')
@feature_ns.route('/Paraloguous/<int:featureid>')
class getParaloguous(Resource):
    ''' Get the list of Features '''
    @feature_ns.doc('list of Features')
    def get(self, featureid):
        ''' List all Features '''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        jsonlist = featureDAO.getParaloguous(featureid)
       
        return jsonlist


@feature_ns.response(404, 'Feature not found')
@feature_ns.route('/List/')
class FeatureWithList(Resource):
    ''' Get the list of Features '''
    @feature_ns.doc('list of Features')
    def post(self):
        ''' List all Features '''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        args = feature_ns.payload['idList']
        featList = featureDAO.getCDSListByIdFeats(args)
        # Extract just CDS
        feats = []
        featsId = {}
        for f in featList:
            featsId[f['idf']]=1
            if f['type'] == 'CDS':
                feats.append(f)
        error = []
        for arg in args:
            if arg not in featsId:
                error.append(arg)
        
        print(error)

        return {'features': feats, 'error': error }


@feature_ns.response(404, 'Feature not found')
@feature_ns.route('/Interact/<int:featureid>')
class Interact(Resource):
    ''' Get the list of Features '''
    @feature_ns.doc('list of Features')
    def get(self, featureid):
        ''' List all Features '''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        featList = featureDAO.getInteractById(featureid)
      

        return featList


@feature_ns.response(404, 'Feature not found')
@feature_ns.route('/isInteract/<int:featureid>')
class isInteract(Resource):
    ''' Get the list of Features '''
    @feature_ns.doc('list of Features')
    def get(self, featureid):
        ''' List all Features '''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        featList = featureDAO.getisInteractId(featureid)
      

        return featList



@feature_ns.response(404, 'Feature not found')
@feature_ns.route('/ListId/')
class FeatureWithList(Resource):
    ''' Get the list of Features '''
    @feature_ns.doc('list of Features')
    def post(self):
        ''' List all Features '''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        args = feature_ns.payload['idList']
        featList = featureDAO.getCDSListById(args)
        # Extract just CDS
        feats = []
        featsId = {}
        for f in featList:
            featsId[f['idf']]=1
            if f['type'] == 'CDS':
                feats.append(f)
        error = []
        for arg in args:
            if arg not in featsId:
                error.append(arg)
        
        print(error)

        return {'features': feats, 'error': error }


@feature_ns.route('/Add/')
class FeatureAdd(Resource):
    ''' Get add Features '''
    @feature_ns.doc('Add Features')
    def post(self):
        ''' List all Features '''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        userid = feature_ns.payload['userid']
        project = feature_ns.payload['project']
        sequenceid = feature_ns.payload['sequenceid']
        complement = feature_ns.payload['complement']
        comments = feature_ns.payload['comments']
        typefeat = feature_ns.payload['type']
        starts = feature_ns.payload['starts']
        stops = feature_ns.payload['stops']
        idfeat = feature_ns.payload['idfeat']

        logger.debug("USERID : "+str(userid))    
        logger.debug("project : "+str(project))    
        logger.debug("sequenceid : "+str(sequenceid))    
        logger.debug("complement : "+str(complement))    
        logger.debug("comments : "+str(comments))    
        logger.debug("typefeat : "+str(typefeat))    
        logger.debug("starts : "+str(starts))    
        logger.debug("stops : "+str(stops))    
        logger.debug("idfeat : "+str(idfeat))    
        featureDAO.addFeatureToProject(userid, project, sequenceid,\
            complement, comments, typefeat, starts, stops, idfeat) 
        return {}

@feature_ns.route('/Remove/')
class FeatureRemove(Resource):
    ''' Get Remove Features '''
    @feature_ns.doc('Remove Features')
    def post(self):
        ''' List all Features '''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        userid = feature_ns.payload['userid']
        project = feature_ns.payload['project']
        featureid = feature_ns.payload['featureid']
        
        logger.debug("USERID : "+str(userid))    
        logger.debug("project : "+str(project))    
        logger.debug("featureid : "+str(featureid))    
        
        featureDAO.removeFeatureToProject(userid, project, featureid) 
        return {}




@feature_ns.response(404, 'Feature not found')
@feature_ns.route('/<int:featureid>')
class Feature(Resource):
    @feature_ns.doc('get a Feature')
    @feature_ns.produces(["application/json", "text/html"])
    def get(self, featureid):
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)
        import json

        feature = featureDAO.getFeatureForId(featureid)
        if 'LocScores' in feature['qualifiers']:
            feature['qualifiers']['LocScores'][0] = json.loads(feature['qualifiers']['LocScores'][0])

        if "html" in request.headers["Accept"]:
            return self.html(feature, featureid)
        return feature

    def html(self, feature, featureid):
        '''html function for sequences'''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        resp = make_response(render_template('html/moreInfo/index.html',
                                             feature=feature,
                                             id=featureid,
                                             product=featureDAO
                                             .getProductForFeature(featureid)))
        # On informe dans le header que le retour sera du html
        resp.headers['Content-Type'] = 'text/html'
        return resp


@feature_ns.response(404, 'Feature not found')
@feature_ns.route('/product/<int:featureid>')
class FeatureProduct(Resource):
    @feature_ns.doc('get a Product for Feature')
    @feature_ns.produces(["application/json"])
    def get(self, featureid):
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        chunks = featureDAO.getProductForFeature(featureid)
        return chunks



@feature_ns.response(404, 'Feature not found')
@feature_ns.route('/location/<int:featureid>')
class FeatureLocation(Resource):
    @feature_ns.doc('get a Location for Feature')
    @feature_ns.produces(["application/json"])
    def get(self, featureid):
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        loc = featureDAO.getLocationForFeature(featureid)
        return loc




@feature_ns.response(404, 'Feature not found')
@feature_ns.route('/feat/<type>/<idfeat>')
class Featurefeat(Resource):
    @feature_ns.doc('get a Feature with type and id_feat')
    @feature_ns.produces(["application/json"])
    def get(self, type, idfeat):
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        feat = featureDAO.getFeatureWithTypeID(type, idfeat)
        return feat


@feature_ns.route('/<type>/<id_feat>')
@feature_ns.response(404, 'Feature not found')
class FeatureWithTypeAndIDFeat(Resource):
    @feature_ns.doc('get a Feature With Type and Id_feat')
    @feature_ns.produces(["application/json", "text/html"])
    def get(self, type, id_feat):
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        # ATTENTIONS PEUT Y AVOIR PLUSIEURS ID
        (nb, f_ids) = featureDAO.getFeaturesIdWithTypeIdFeat(type, id_feat)
        # Si on a qu'un Feature pour ce id_feat
        if nb == 1:
            feature = featureDAO.getFeatureForId(f_ids[0])
            if "html" in request.headers["Accept"]:
                return self.html(feature, f_ids[0])
            return feature
        else:

            features = featureDAO.getFeatureListByIds(f_ids)
            if "html" in request.headers["Accept"]:
                return self.htmlList(features, type, id_feat)
            return features

    def html(self, feature, featureid):
        '''html function for sequences'''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        # reparation du retour html avec le template
        resp = make_response(render_template('html/moreInfo/index.html',
                                             feature=feature,
                                             id=featureid))
        # On informe dans le header que le retour sera du html
        resp.headers['Content-Type'] = 'text/html'
        return resp

    def htmlList(self, features, type, id_feat):
        '''html function for sequences'''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        # reparation du retour html avec le template
        resp = make_response(render_template('html/moreInfo/featureList.html',
                                             featList=features,
                                             type=type,
                                             id_feat=id_feat))
        # On informe dans le header que le retour sera du html
        resp.headers['Content-Type'] = 'text/html'
        return resp
#############################

@genecard_ns.route('/<id_feat>')
@feature_ns.route('/<id_feat>')
@feature_ns.response(404, 'Feature not found')
class FeatureIDFeat(Resource):
    @feature_ns.doc('get a Feature With Id_feat')
    @feature_ns.produces(["application/json", "text/html"])
    def get(self, id_feat):
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        mapIt = False
        if request.args.get('map') == '1':
            mapIt = True
        # ATTENTIONS PEUT Y AVOIR PLUSIEURS ID
        (nb, f_ids) = featureDAO.getFeaturesIdWithIdFeat(id_feat)
        # Si on a qu'un Feature pour ce id_feat
        if nb == 1:
            feature = featureDAO.getFeatureForId(f_ids[0])
            if "html" in request.headers["Accept"]:
              
                if 'GeneCard' in  request.base_url: 
                            url = 'html/moreInfo/GeneCard/indexCDSmRNA.html'
                            resp = make_response(
                                render_template(url,
                                                cds=feature,
                                                mrna=feature,
                                                id=feature['id'],
                                                id_feat=feature['idf'],
                                                product=featureDAO
                                                .getProductForFeature(
                                                    feature['id']),
                                                sequence = chunkSequence(sequenceDAO.getSequenceDNAFromTo(
                                                    feature['sid'],
                                                    feature['start'],
                                                    feature['stop'])),                                                    
                                                homeURL=homeURL,
                                                webURL=URL+URLPATH+"/"+SUBURIPATH+"/index.html"
                                                ))
                            resp.headers['Content-Type'] = 'text/html'
                            return resp
                return self.html(feature, f_ids[0], mapIt)
            return featureDAO.getFeatureListByIds(f_ids)
        else:
            features = featureDAO.getFeatureListByIds(f_ids)
            
            if "html" in request.headers["Accept"]:
                # reparation du retour html avec le template
                #  ICICICI!!!!
                if ((features[0]['type'] == 'CDS') and
                   (features[1]['type'] == 'mRNA')):
                    CDS = featureDAO.getFeatureForId(features[0]['id'])
                    mRNA = featureDAO.getFeatureForId(features[1]['id'])
                    if mapIt:
                        return redirect('/'+SUBURIPATH+'/index.html?id_feat=' +
                                        str(features[0]['idf']+'&fId=' +
                                            str(features[0]['id'])))
                    else:
                        url = 'html/moreInfo/indexCDSmRNA.html'
                        if 'GeneCard' in  request.base_url: 
                            url = 'html/moreInfo/GeneCard/indexCDSmRNA.html'
                        logger.debug(features[0])
                        resp = make_response(
                                render_template(url,
                                                cds=CDS,
                                                mrna=mRNA,
                                                id=features[0]['id'],
                                                id_feat=features[0]['idf'],
                                                product=featureDAO
                                                .getProductForFeature(
                                                    features[0]['id']),
                                                homeURL=homeURL,
                                                webURL=URL+URLPATH+"/"+SUBURIPATH+"/index.html"
                                                ))
    # On informe dans le header que le retour sera du html
                        resp.headers['Content-Type'] = 'text/html'
                        return resp
                return self.htmlList(features, id_feat, mapIt)
            return features

    def html(self, feature, featureid, mapIt=False):
        '''html FeatureID'''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        # reparation du retour html avec le template
        if mapIt:
            return redirect('/'+SUBURIPATH+'/index.html?id_feat=' +
                            str(feature['id_feat']+'&fId=' +
                            str(featureid)))

#            return redirect('/flagdb/index.html?speciesidx=' +
#                            str(speciesIdx))
#                                speciesIdx=request.form.get('speciesIdx',
#                                                            str(speciesIdx))))
        else:
            resp = make_response(render_template('html/moreInfo/index_genecard.html',
                                                 feature=feature,
                                                 homeURL=homeURL,
                                                 id=featureid), 
                                                 )
        # On informe dans le header que le retour sera du html
        resp.headers['Content-Type'] = 'text/html'
        return resp

    def htmlList(self, features, id_feat, mapIt=False):
        '''html function for Feature List'''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        # reparation du retour html avec le template
        resp = make_response(render_template('html/moreInfo/featureList.html',
                                             featList=features,
                                             id_feat=id_feat,
                                              homeURL=homeURL,
                                             mapIt=mapIt),
                                            )
        # On informe dans le header que le retour sera du html
        resp.headers['Content-Type'] = 'text/html'
        return resp


@feature_ns.route('/')
@feature_ns.response(404, 'Feature not found')
class FeatureList(Resource):
    ''' Get the list of Features '''
    @feature_ns.doc('list of Features')
    @feature_ns.expect(parser)
    @feature_ns.marshal_with(feature, code=201)
    def post(self):
        ''' List all Features for a Sequence ID from Start to Stop'''
        name = request.cookies.get('gender')
        FLAGdb.setGender(name)

        args = parser.parse_args()
        # is it for user features?
        logger.debug(args)
        if args['userid'] != None:
            featureList = featureDAO.getUserFeatureForSequence(
                                                        args['sequenceId'],
                                                       args['userid'],
                                                       args['project'],
                                                       args['start'],
                                                       args['stop'])
        else:
            featureList = featureDAO.getFeatureForSequence(args['sequenceId'],
                                                       args['type'],
                                                       args['start'],
                                                       args['stop'])
        return featureList
