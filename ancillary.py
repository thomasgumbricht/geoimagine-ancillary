'''
Created on 23 feb. 2018

@author: thomasgumbricht
'''

from os import path, makedirs, remove, system
from sys import exit
from geoimagine.kartturmain import Composition, RegionLayer
import geoimagine.support.karttur_dt as mj_dt 
from geoimagine.ancillary import ancillary_import_v73 as ancillary_import
from geoimagine.gis.gis import GetVectorProjection
from shutil import copyfile, copyfileobj

class AncilComposition(Composition): 
    def __init__(self,compD):
        """The constructor requires a dict {datadir, datafile, compyright,title,accessdate,theme,subtheme,label,version,
        dataset,product,datapath,metapath,dataurl,metaurl}.""" 
        for key, value in compD.items():
            if not hasattr(self, key):
                setattr(self, key, value) 

    def SetDataFiletype(self,filetype):
        self.datafiletype = filetype
   
    def SetPath(self,mainpath):
        self.FP = path.join(mainpath,self.datadir)
        
class AncillaryLayer(RegionLayer):
    'Common base class for Ancillary spatial data of different formats'
    def __init__(self, comp, rawComp, region, acqdate, process):
        self.comp = comp
        self.layertype = 'ancillary'
        self.comp = comp
        self.rawComp = rawComp   
        self.acqdate = acqdate
        self.ancilid = '%(p)s.%(d)s' %{'p':process.dsid,'d':rawComp.datalayer}  

        self.datadir = rawComp.datadir
        self.metadir = rawComp.metadir

        if len(rawComp.accessdate) >= 4:
            self.accessdate = mj_dt.yyyymmddDate(rawComp.accessdate) 
        else:
            self.accessdate = mj_dt.Today()
        self.createdate = mj_dt.Today()
        self.FPN = False
        if self.comp.hdr == 'ers':
            FN = '%(f)s.%(d)s.%(e)s' %{'f':self.rawComp.datafile,'d':self.comp.dat,'e':self.comp.hdr}
        else:
            if len(self.comp.hdr) > 0:
                FN = '%(f)s.%(e)s' %{'f':self.rawComp.datafile,'e':self.comp.hdr}
            elif len(self.dat) > 0:
                FN = '%(f)s.%(e)s' %{'f':self.rawComp.datafile,'e':self.comp.dat}
            elif len(self.comp.hdr) == 0 and len(self.comp.dat) == 0:
                #The data is in a folder, e.g. ArcVew raster format
                FN = self.datafile
            else:
                print (self.datafile,self.dat,self.hdr)
                exit('FN error in ancillary')
        #For some reason os.path.join does not work here
        FP = '%s/%s' %(self.comp.mainpath,self.rawComp.datadir)
        #FP = os.path.join(ancilC.mainpath,self.datadir)  
        FPN = path.join(FP,FN)
        self.FPN = FPN
        #add the path except for volume and mai folder + to add to db
        dbFP = path.join(self.rawComp.datadir,FN)
        dbFP = '../%(d)s' %{'d':dbFP}
        self.dbFP = dbFP
        if self.comp.hdr in ['zip']:
            self.zip = 'zip'
        elif self.comp.hdr in ['tar.gz','tar']:
            self.zip = 'tar.gz'
        else:
            self.zip = False 
            
    def _UnZip(self):
        import zipfile
        zipFP,zipFN = path.split(self.FPN)
        tempFP = path.join(zipFP,'ziptmp')
        self.tempFP = tempFP
        if not path.isdir(tempFP):
            makedirs(tempFP)
        zipF = zipfile.ZipFile(self.FPN, "r")
        self.FPN = False
        #compstr is the filetype to look for in the zipfile

        compstr = '.%(e)s' %{'e':self.comp.dat}
        for fname in zipF.namelist():
            fname = path.basename(fname) 
            # skip directories
            if not fname:
                continue
            #get the fname components
            stem,ext = path.splitext(fname)
            fnameout = '%(s)s%(e)s' %{'s':stem, 'e':ext.lower()}
            #check if thsi file is of the data type expected
            if ext.lower() == compstr:
                if self.FPN: 
                    exitstr = 'EXITING - ancullary zip archice %(s)s contains multiple data files, you must unzip and give data filenames' %{'s':zipFN}
                    exit(exitstr)
                else:
                    #Change hdr type
                    self.comp.hdr = self.comp.dat
                    self.FPN = path.join(tempFP,fnameout)
            # copy file (taken from zipfile's extract)
            source = zipF.open(fname)
            target = file(path.join(tempFP, fnameout), "wb")
            with source, target:
                copyfileobj(source, target)
                if self.verbose:
                    print ('shutil.copyfileobj',source, target)
        if not self.FPN:
            exit('Exiting, no supported file type found in the zip file')
        if not path.isfile(self.FPN):
            exit('Something wrong with the unzipping of ancillary data')
      
    def _UnTar(self):
        import tarfile
        if self.comp.hdr == 'tar.gz':
            tarFPN = path.splitext(self.FPN)[0]
            tarFP,tarFN = path.split(tarFPN)
            if not path.isfile(tarFPN):
                cmd = 'cd %(d)s; gunzip -dk %(s)s' %{'d':tarFP, 's':path.split(self.FPN)[1]}
                os.system(cmd)
        else:
            tarFPN = self.FPN   
        tempFP = path.join(tarFP,'tartmp')
        self.tempFP = tempFP
        if not path.isdir(tempFP):
            makedirs(tempFP) 
        #compstr is the filetype to look for in the zipfile
  
        compstr = '.%(e)s' %{'e':self.comp.dat}
        tarF = tarfile.TarFile(tarFPN, "r")
        self.FPN = False
        for fname in tarF.getnames():
            fname = path.basename(fname) 
            # skip directories
            if not fname:
                continue
            #get the fname components
            ext = path.splitext(fname)[1]
            #check if this file is of the data type expected
            if ext.lower() == compstr:
                if self.FPN: 
                    exitstr = 'EXITING - ancillary zip archice %(s)s contains multiple data files, you must unzip and give data filenames' %{'s':tarFN}
                    exit(exitstr)
                else:
                    #Change hdr type
                    self.comp.hdr = self.comp.dat
                    self.FPN = path.join(tempFP,fname)
            tarF.extract(fname, tempFP)
        if self.zip == 'tar.gz':
            #remove the tar file
            remove(tarFPN)
        if not path.isfile(self.FPN):
            exit('Someting wrong with the unzipping of ancillary data')
            
class ProcessAncillary:
    'class for addint ancillary data'   
    def __init__(self, process,session,verbose):
        self.verbose = verbose
        self.process = process            
        #direct to subprocess
        print (self.process.proc.processid)

        if self.process.proc.processid == 'organizeancillary':
            self._OrganizeAncillary(session)
        elif self.process.proc.processid == 'OrganizeGrace':
            #This is the same as organizeancillary but from Grace
            self._OrganizeAncillary(session)
        elif self.process.proc.processid == 'anciltxttodb':
            self._AncilTxtToDb(session)
        else:
            print (self.process.proc.processid)
            BKAJHF
            exit('Ancillary process not recognized in ProcessAncillary')
                       
    def _OrganizeAncillary(self,session):

        self.dsid  = '%(i)s.%(c)s.%(v)s.%(r)s' %{'i':self.process.proc.paramsD['instid'],'c':self.process.proc.paramsD['dsname'],'v':self.process.proc.paramsD['dsversion'],'r':self.process.proc.paramsD['regionid']}
        self.dsplotid = '%(i)s.%(c)s.%(r)s' %{'i':self.process.proc.paramsD['instid'],'c':self.process.proc.paramsD['dsname'],'r':self.process.proc.paramsD['regionid']}
        #check that the region exists as default, otherwise it must first be created
        rec = session._SelectDefaultRegion(self.process.proc.paramsD['regionid'])

        if rec == None:
            exitstr ='Organizing ancillary or specimen data requires an existing region: %s does not exist' %(self.process.params.regionid)
            exit(exitstr)
        print (self.process.dstLayerD)
        if len(self.process.dstLayerD) == 0:
            exitstr = 'EXITING, no locus defined in Ancillary.OrganizeAncillary'
            exit(exitstr)
        for locus in self.process.dstLayerD:
            print ('locus',locus)
            if len(self.process.dstLayerD[locus]) == 0:
                exitstr = 'EXITING, no dates defined in Ancillary.OrganizeAncillary'
                exit(exitstr)
            
            for datum in self.process.dstLayerD[locus]:
                print ('datum',datum)
                print ('self.process.dstLayerD',self.process.dstLayerD)
                print ('self.process.dstLayerD[locus][datum]',self.process.dstLayerD[locus][datum])
                if len(self.process.dstLayerD[locus][datum]) == 0:
                    exitstr = 'EXITING, no compositions defined in Ancillary.OrganizeAncillary'
                    print (exitstr)
                    BALLE
                    exit(exitstr)

                for comp in self.process.dstLayerD[locus][datum]:
                    #print ('comp',comp)
                    self.dstLayer = self.process.dstLayerD[locus][datum][comp]
                    #print ('self.dstLayer',self.dstLayer)
                    #print ('self.dstLayer.FPN',self.dstLayer.FPN)
                    #print ('self.dstLayer.locus.locus',self.dstLayer.locus.locus)

                    if not self.dstLayer._Exists() or self.process.proc.overwrite:         
                        #Create the source layer
                        self._SrcLayer(comp)
                        if not path.isfile(self.srcFPN):
                            warnstr = 'The ancillary source file %(fpn)s can not be found, skipping' %{'fpn':self.srcFPN}
                            print (warnstr)
                            BALLE
                            continue

                        if self.dstLayer.comp.celltype.lower() in ['none','csv','txt']:
                            self._ImportText()
                        elif self.dstLayer.comp.celltype.lower() == 'vector':
                            self._ImportVector()
                        else:
                            self._ImportRaster() 
                            
                    self._UpdateDB(self.dstLayer,comp,session)  
  
    def _SrcLayer(self,comp):
        '''
        '''
        self.srcrawD = self.process.proc.srcraw.paramsD[comp]

        print ('self.process.srcpath.hdrfiletype',self.process.srcpath.hdrfiletype)
  
        if self.process.proc.srcpathD['hdrfiletype'][0] == '.':
            ext = self.process.proc.srcpathD['hdrfiletype']
        else:
            ext = '.%s' %(self.process.proc.srcpathD['hdrfiletype'])
        self.srcFN = '%(fn)s%(e)s' %{'fn':self.srcrawD['datafile'],'e':ext}            
        self.srcFP = path.join('/Volumes',self.process.proc.srcpathD['volume'], self.srcrawD['datadir'])
        self.srcFPN = path.join(self.srcFP,self.srcFN)

    def _ImportVector(self):
        if self.verbose:
            printstr = '    importing vector: %(srcfpn)s\n        to: %(dstfpn)s' %{'srcfpn':self.srcFPN, 'dstfpn':self.dstLayer.FPN}
            print (printstr)
            
        spatialRef = GetVectorProjection(self.srcFPN)
        gdalcmd = '/Library/Frameworks/GDAL.framework/Programs/ogr2ogr -skipfailures'
        if not spatialRef.epsg:
            gdalcmd = ' %(s1)s -a_srs EPSG:%(epsg)d' %{'s1':gdalcmd, 'epsg': int(self.process.parameters.epsg)}
        else:
            if spatialRef.epsg == 4326:
                pass
            else:
                gdalcmd = ' %(s1)s -t_srs EPSG:4326 ' %{'s1':gdalcmd}
        gdalcmd = ' %(s1)s %(dst)s %(src)s' %{'s1':gdalcmd, 'dst': self.dstLayer.FPN, 'src':self.srcFPN}
        system(gdalcmd)
        #subprocess.check_call(gdalcmd)

    def _ImportRaster(self):
        if self.verbose:
            print ('    Importing raster')

        if self.process.srcpath.hdrfiletype.lower() == 'lis':
            '''This is a very special format, only applies to gghydro'''
            ancillary_import.GGHtranslate(self.Lin.FPN,self.Lout.FPN,self.Lout.comp.celltype,self.Lout.comp.cellnull,False)
            
        elif self.process.srcpath.hdrfiletype.lower() == '1x1':
            '''This is a very special format, only applies to stillwell'''
            ancillary_import.StillwellTranslate(self.Lin.FPN,self.Lout.FPN,self.Lout.comp.celltype,self.Lout.comp.cellnull,False)
            
        elif self.process.srcpath.datfiletype.lower() == 'trmm':
            '''This is a very special format, only applies to TRMM data with north to the right'''
            ancillary_import.TRMMTranslate(self.Lout.comp, self.Lin.FPN,self.Lout.FPN,False)
            
        elif self.process.params.importdef.lower() == 'grace':
            '''This is a very special format, only applies to TRMM data with north to the right'''
            #ancillary_import.TRMMTranslate(self.Lout.comp, self.Lin.FPN,self.Lout.FPN,False)
            if not self.dstLayer.comp.celltype == 'Float32':
                BALLE
            if not self.dstLayer.comp.cellnull == 32767:
                print (self.dstLayer.comp.cellnull)
                BALLE
  
            ancillary_import.GRACETranslate(self.srcFPN, self.dstLayer.FPN, self.dstLayer.comp, False)
        else:
            print (self.Lin.FPN,self.Lout.FPN)
            print (self.Lin.comp.dat)
            print (self.Lin.comp.band)
            BALLE
        
            
    def _ImportText(self):
        if self.verbose:
            printstr = '    importing text file: %(srcfpn)s\n        to: %(dstfpn)s' %{'srcfpn':self.srcFPN, 'dstfpn':self.dstLayer.FPN}
            print (printstr)
        copyfile(self.srcFPN,self.dstLayer.FPN)

    def _StringReplace(self,strObj,compo,searchStr,replaceStr):
        if strObj == 'folder' and hasattr(compo, 'folder'):
            compo.folder = compo.folder.replace(searchStr,replaceStr,1)
        if strObj == 'band' and hasattr(compo, 'band'):
            compo.band = compo.band.replace(searchStr,replaceStr,1)
        if strObj == 'prefix' and hasattr(compo, 'prefix'):
            compo.prefix = compo.prefix.replace(searchStr,replaceStr,1)
        if strObj == 'suffix' and hasattr(compo, 'suffix'):
            compo.suffix = compo.suffix.replace(searchStr,replaceStr,1)
        if strObj == 'yyyydoy' and hasattr(compo, 'yyyydoy'):
            compo.yyyydoy = compo.yyyydoy.replace(searchStr,replaceStr,1)
        if strObj == 'datadir' and hasattr(compo, 'datadir'):
            compo.datadir = compo.datadir.replace(searchStr,replaceStr,1)
        if strObj == 'datafile' and hasattr(compo, 'datafile'):
            compo.datafile = compo.datafile.replace(searchStr,replaceStr,1)
        if strObj == 'metapath' and hasattr(compo, 'metapath'):
            compo.metapath = compo.metapath.replace(searchStr,replaceStr,1)
        if strObj == 'dataurl' and hasattr(compo, 'dataurl'):
            compo.dataurl = compo.dataurl.replace(searchStr,replaceStr,1)
        if strObj == 'metaurl' and hasattr(compo, 'metaurl'):
            compo.metaurl = compo.metaurl.replace(searchStr,replaceStr,1)
        if strObj == 'dataset' and hasattr(compo, 'dataset'):
            compo.dataset = compo.dataset.replace(searchStr,replaceStr,1)
        if strObj == 'title' and hasattr(compo, 'title'):
            compo.title = compo.title.replace(searchStr,replaceStr,1)
        if strObj == 'label' and hasattr(compo, 'label'):
            compo.label = compo.label.replace(searchStr,replaceStr,1)
           
    def _AncilTxtToDb(self,session):
        import csv
        for locus in self.process.srcLayerD:
            print ('locus',locus)
            for datum in self.process.srcLayerD[locus]:
                print ('datum',datum)
                print ('self.process.dstLayerD',self.process.srcLayerD)
                print ('self.process.dstLayerD[locus][datum]',self.process.srcLayerD[locus][datum])
                for comp in self.process.srcLayerD[locus][datum]:

                    self.srcLayer = self.process.srcLayerD[locus][datum][comp]
                    print (self.srcLayer)
                    print (self.srcLayer.FPN)
                    print (self.process.srcLayerD[locus][datum][comp].comp.band)

                    print (self.process.params.template)

                    if not path.isfile(self.srcLayer.FPN):
                        warnstr = 'The ancillary source file %(fpn)s can not be found, skipping' %{'fpn':self.srcLayer.FPN}
                        print (warnstr)
                        BALLE
                        continue
                    queryL = []
                    if self.process.params.template == 'climateindex':
                        with open(self.srcLayer.FPN) as f:
                            reader = csv.reader(f, delimiter=' ', skipinitialspace = True)
                            startyaer, endyear = next(reader)
                            for row in reader:
                                y = row[0]
                                for m in range(1,13):
                                    acqdate = mj_dt.yyyy_mm_dd_Date(y,m,1)
                                    acqdatestr = mj_dt.DateToStrDate(acqdate)[0:6]
                                    value = row[m]
                                    if y == endyear and value[0:3] == '-99':
                                        continue
                                    queryL.append({'index':self.process.srcLayerD[locus][datum][comp].comp.band, 'acqdate':acqdate,'acqdatestr':acqdatestr,'value':value})
                                if y == endyear:
                                    break

                    session._InsertClimateIndex(queryL)
      
    def _UpdateDB(self, layer, comp, session):
        
        srcrawD = self.process.proc.srcraw.paramsD[comp]
        if self.process.proc.userProj.system == 'specimen':
            system = 'specimen'
        else:
            system = 'ancillary'

        session._ManageAncilDS(system, self.process.proc.paramsD, self.dsid, self.process.proc.overwrite, self.process.proc.delete)
        session._InsertCompDef(layer.comp,srcrawD['title'],srcrawD['title'])
        session._LinkDsCompid(self.dsid, layer.comp.compid, self.process.proc.overwrite, self.process.proc.delete)
        session._InsertCompProd(layer.comp)

        '''TGTODO ADD function
        if self.Lout.comp.hdr not in ['shp','.shp','csv','.csv']:
            ConnAncillary.ManageAncillGeo(self.Lout,self.process.delete,self.process.overwrite)
        '''
        session._InsertLayer(layer, self.process.proc.overwrite, self.process.proc.delete)
        #ConnAncillary.ManageAncilLayer(self.Lout,delete,overwrite) 
        #ConnAncillary.ManageAncilMeta(self.Lin,delete,overwrite)