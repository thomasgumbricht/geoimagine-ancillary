'''
Created on 29 Aug 2018

@author: thomasgumbricht
'''

from osgeo import ogr,osr
from sys import exit
import os, glob, zipfile
import geoimagine.gis.mj_gis_v80 as mj_gis
from operator import itemgetter

def MGRSorginalNGAtilesOS(FP,dstFN):
    '''Expects a folder with all NGA original tiles assembled per Grid Zone Designator in shapefiles
    '''
    x = 0
    dstFPN = os.path.join(FP,dstFN)
    for FPN in glob.glob(os.path.join(FP, '*.zip')):
        if 'arctic' in FPN.lower():
            continue
        FN = os.path.split(FPN)[1]
        #print ('FN',FN)
        zipFPN = os.path.join(FP,FN)
        #print ('zipFPN',zipFPN)

        zipFP = os.path.splitext(zipFPN)[0]
        #print (zipFP)
        '''
        zip_ref = zipfile.ZipFile(zipFPN, 'r')
        zip_ref.extractall(zipFP)
        zip_ref.close()
        '''
        #shpFN = '%(fn)s.shp' %{'fn':os.path.splitext(FN)[0]}
        #print (shpFN)
        #srcFPN = os.path.join(zipFP,shpFN)
        #print (srcFPN)
        for shpFPN in glob.glob(os.path.join(zipFP, '*.shp')):

            #Open and process the shp file
            if x == 0:
                cmd = '/Library/Frameworks/GDAL.framework/Versions/2.2/Programs/ogr2ogr -skipfailures %(dst)s %(src)s' %{'dst':dstFPN, 'src':shpFPN}
                #print (cmd)
                os.system(cmd)
                #BALLE
            else:
                cmd = '/Library/Frameworks/GDAL.framework/Versions/2.2/Programs/ogr2ogr -skipfailures -append %(dst)s %(src)s' %{'dst':dstFPN, 'src':shpFPN}
                os.system(cmd)

        '''
        else:
            exitstr = 'File not found %(f)s' %{'f':srcFPN}
            exit(exitstr)
        '''
        x+=1
        
def MGRSorginalNGAtiles(FP,dstFN):
    '''Expects a folder with all NGA original tiles assembled per Grid Zone Designator in shapefiles
    '''
    from math import sqrt
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dstFPN = os.path.join(FP,dstFN)
    if os.path.exists(dstFPN):
        driver.DeleteDataSource(dstFPN)
    dstDS = driver.CreateDataSource(dstFPN)
    dst_srs = osr.SpatialReference()
    dst_srs.ImportFromEPSG(4326)
    dstLayer = dstDS.CreateLayer('tiles', dst_srs, geom_type=ogr.wkbPolygon)
    dst_srs = mj_gis.MjProj()
    dst_srs.SetFromEPSG(4326)
    mgrs_field = ogr.FieldDefn("MGRS", ogr.OFTString)
    mgrs_field.SetWidth(8)
    dstLayer.CreateField(mgrs_field) 
    
    gzd_field = ogr.FieldDefn("GZD", ogr.OFTString)
    gzd_field.SetWidth(8)
    dstLayer.CreateField(gzd_field) 

    for FPN in glob.glob(os.path.join(FP, '*.zip')):
        if 'arctic' in FPN.lower():
            continue
        FN = os.path.split(FPN)[1]
        #print ('FN',FN)
        zipFPN = os.path.join(FP,FN)
        #print ('zipFPN',zipFPN)

        zipFP = os.path.splitext(zipFPN)[0]
        print (zipFP)
        '''
        zip_ref = zipfile.ZipFile(zipFPN, 'r')
        zip_ref.extractall(zipFP)
        zip_ref.close()
        '''

        for root, dirs, files in os.walk(zipFP):
            shpFPN = ''
            for name in files:
                if name.endswith("shp"):
                    shpFPN = os.path.join(root, name)

            if os.path.isfile(shpFPN):
                #print ('    root',root)
                #print ('    ',shpFPN)
                srcDS = driver.Open(shpFPN, 0)
                srcLayer = srcDS.GetLayer()
                spatialRef = srcLayer.GetSpatialRef()
                src_srs = mj_gis.MjProj()
                src_srs.SetProj(spatialRef)
                mgrsL = []
                layerDefinition = srcLayer.GetLayerDefn()
                fieldNameL =[]
                for i in range(layerDefinition.GetFieldCount()):
                    fieldNameL.append( layerDefinition.GetFieldDefn(i).GetName() )
                for feature in srcLayer:
                    if 'MGRS' in fieldNameL:
                        MGRS = feature.GetField("MGRS")
                    elif ('GZD' in fieldNameL and '100kmSQ_ID' in fieldNameL):
                        MGRS = '%(gzd)s%(EN)s' %{'gzd':feature.GetField("GZD"), 'EN':feature.GetField("100kmSQ_ID")}
                    else:
                        print (fieldNameL)
                        UNKNOWN
                    srcGeom = mj_gis.Geometry()
                    srcGeom.GeomFromFeature(feature)
                    
                    if srcGeom.shapelyGeom.geom_type == 'MultiPolygon':
                        print ( '    ', MGRS, srcGeom.shapelyGeom.geom_type )
                        #print (srcGeom.shapelyGeom)
                        polygeom = mj_gis.Geometry()
                        polygeom.MultPolyToSinglePoly(srcGeom.shapelyGeom)
                        #print (polygeom.shapelyGeom)
                        #BALLE
                        #for poly in srcGeom.shapelyGeom:
                        print (polygeom.shapelyGeom)
                        edgeL = list(polygeom.shapelyGeom.exterior.coords)
                    else:
                        edgeL = list(srcGeom.shapelyGeom.exterior.coords)
                        
                    boundsL = srcGeom.BoundsToPoly()

                    cL = ['ul','ur','lr','ll']
                    cornerD = {}
                    cornerL = list(boundsL.exterior.coords)

                    for i, c in enumerate (cL):
                        dist = 100000
                        for j, d in enumerate (edgeL):
                            testdist = sqrt((cornerL[i][0]-d[0])**2+(cornerL[i][1]-d[1])**2)
                            if testdist < dist:
                                dp = j
                                dist = testdist
                        if dp >= len(edgeL) or dist == 100000:
                            print ('MGRS',MGRS)
                            print ('boundsL',boundsL)
                            print ('dist',dist)
                            print ('dp',dp)
                            print ('edgeL',edgeL)
                            print ('edgeL',len(edgeL))
                            BALLE
                        cornerD[c] = edgeL[dp]
                    #print (cornerD)
                    '''
                    print ('corner',cornerL)
                    FISK
                    print ('MGRS',MGRS)
                    print ('len',len(cornerL))
                    
                    
                    if len(cornerL) > 5:
                        tol = 1
                        factor = 1.1
                        while len(cornerL) > 5:
                            tol *= factor
                            corners = srcGeom.shapelyGeom.simplify(tol, preserve_topology=False)
                            if not corners.is_empty:
                                cornerL = list(corners.exterior.coords)
                            if tol > 1000000:
                                factor = 0.9

                   
                    coords = list(corners.exterior.coords)
                    if len(coords) == 5:
                        coords.pop(4)
                    if len(coords) != 4:
                        print ('coords',coords)
                        ERROR
    
                    ysortedCoords = sorted(coords, key=itemgetter(1))
                    minyCoords = ysortedCoords[0:2]
                    maxyCoords = ysortedCoords[2:4]
                    cornerD = {}
                    if minyCoords[0][0] < minyCoords[1][0]:
                        cornerD['ll'] = minyCoords[0]
                        cornerD['lr'] = minyCoords[1]
                    else:
                        cornerD['ll'] = minyCoords[1]
                        cornerD['lr'] = minyCoords[0]
                                    
                    if maxyCoords[0][0] < maxyCoords[1][0]:
                        cornerD['ul'] = maxyCoords[0]
                        cornerD['ur'] = maxyCoords[1]
                    else:
                        cornerD['ul'] = maxyCoords[1]
                        cornerD['ur'] = maxyCoords[0]
                    
                    #print (cornerD)
                    '''

                    cornerPtL = (cornerD['ul'],cornerD['ur'],cornerD['lr'],cornerD['ll'])
                    
                    cornerGeom = mj_gis.Geometry()
                    cornerGeom.PointsToPolygonGeom(cornerPtL)
                    #print (cornerGeom.shapelyGeom)
                    #cornerGeom.SetShapelyGeom(corners)
                    
                    #print ('corners',corners)
                    #cornerGeom = mj_gis.Geometry()
                    #cornerGeom.SetShapelyGeom(corners)

                    #BALLE
                    #bounds = srcGeom.BoundsToPoly() 
                    #print (bounds)
                    '''
                    featureDefn = dstLayer.GetLayerDefn()
                    transfeat = ogr.Feature(featureDefn)       
                    transfeat.SetGeometry(srcGeom.ogrGeom)
                    transfeat.SetField("MGRS", 'orig')
                            
                    transfeat.SetField("GZD", 'orig')
                    dstLayer.CreateFeature(transfeat)
                    transfeat = None
                    
                    
                    featureDefn = dstLayer.GetLayerDefn()
                    transfeat = ogr.Feature(featureDefn) 
                    #      
                    transfeat.SetGeometry(cornerGeom.ogrGeom)
                    transfeat.SetField("MGRS", 'bounds')
                            
                    transfeat.SetField("GZD", 'bounds')
                    dstLayer.CreateFeature(transfeat)
                    transfeat = None
                            
                    srcDS = None 
                    dstDS = None      
                            
                    BALLE
                    '''
                    dstGeom = src_srs.ReprojectGeom(cornerGeom, dst_srs)
                    
                    #bounds = dstGeom.BoundsToPoly() 
                    #print (bounds)

                    if 'MGRS' in fieldNameL:
                        MGRS = feature.GetField("MGRS")
                    elif ('GZD' in fieldNameL and '100kmSQ_ID' in fieldNameL):
                        MGRS = '%(gzd)s%(EN)s' %{'gzd':feature.GetField("GZD"), 'EN':feature.GetField("100kmSQ_ID")}
                    else:
                        print (fieldNameL)
                        UNKNOWN
                    if MGRS in mgrsL:
                        #continue
                        pass
                    else:
                        boundsL = list(dstGeom.shapelyGeom.exterior.coords) 
                        print ( len(boundsL),boundsL, abs(int(round(boundsL[0][0]))))
                        if MGRS[0:2] == '01' and abs( int( round( boundsL[0][0] ) ) ) == 180:
                            print ('    fixing')
                            boundsL[0] = (-180,boundsL[0][1])
                            boundsL[3] = (-180,boundsL[3][1])
                            if len(boundsL) > 4: 
                                boundsL[4] = (-180,boundsL[4][1])
                            
                        if MGRS[0:2] == '60' and abs(int(round(boundsL[0][1]))) == 180:
                            boundsL[1] = (180,boundsL[1][1])
                            boundsL[2] = (180,boundsL[2][1])
                        lonL = [p[0] for p in boundsL]
                        latL = [p[1] for p in boundsL]
                        if boundsL[0][0] > boundsL[1][0]:
                            printstr = 'Skipping %(mgrs)s' %{'mgrs':MGRS}  
                            #if MGRS[0:2] == '01':
                            #    boundsL[0][0] = boundsL[3][0] = boundsL[4][0] = -180
                            print (printstr)
                            print (boundsL)
                            BALLE

                        elif boundsL[1][0] - boundsL[0][0] > 6.0001:
                            printstr = 'Skipping %(mgrs)s' %{'mgrs':MGRS}
                            
                            print (printstr)
                            print (boundsL[0][0], boundsL[1][0])
                            print (boundsL)
 
                            testGeom = mj_gis.Geometry()
                            testGeom.GeomFromFeature(feature)
                            print (testGeom.shapelyGeom)

                            BALLE
                        elif round(min(lonL),2) < -180 or  round(max(lonL),2) > 180:
                            printstr = 'Skipping %(mgrs)s' %{'mgrs':MGRS}
                            
                            print (printstr)
                            print (boundsL)
                            BALLE
                        elif round(min(latL),2) < -90 or  round(max(latL),2) > 90:
                            printstr = 'Skipping %(mgrs)s' %{'mgrs':MGRS}
                            
                            print (printstr)
                            print (boundsL)
                            BALLE
                            
                        else:
    
                            featureDefn = dstLayer.GetLayerDefn()
                            transfeat = ogr.Feature(featureDefn)
                            
                            transfeat.SetGeometry(dstGeom.ogrGeom)
                            
                            transfeat.SetField("MGRS", MGRS)
                            
                            transfeat.SetField("GZD", feature.GetField("GZD"))
                            dstLayer.CreateFeature(transfeat)
                            transfeat = None
                            mgrsL.append(MGRS)
                # Save and close DataSource
                srcDS = None
    dstDS = None
                
def MGRSorginalNGAtilesV2(FP,dstFN):
    '''Expects a folder with all NGA original tiles assembled per Grid Zone Designator in shapefiles
    '''
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dstFPN = os.path.join(FP,dstFN)
    if os.path.exists(dstFPN):
        driver.DeleteDataSource(dstFPN)
    dstDS = driver.CreateDataSource(dstFPN)
    dst_srs = osr.SpatialReference()
    dst_srs.ImportFromEPSG(4326)
    dstLayer = dstDS.CreateLayer('tiles', dst_srs, geom_type=ogr.wkbPolygon)
    dst_srs = mj_gis.MjProj()
    dst_srs.SetFromEPSG(4326)
    mgrs_field = ogr.FieldDefn("MGRS", ogr.OFTString)
    mgrs_field.SetWidth(8)
    dstLayer.CreateField(mgrs_field) 
    
    gzd_field = ogr.FieldDefn("GZD", ogr.OFTString)
    gzd_field.SetWidth(8)
    dstLayer.CreateField(gzd_field) 
    
    area_field = ogr.FieldDefn("AREAKM2", ogr.OFTReal)
    area_field.SetWidth(8)
    area_field.SetPrecision(2)
    dstLayer.CreateField(area_field) 
    
    easting_field = ogr.FieldDefn("EASTING", ogr.OFTString)
    easting_field.SetWidth(12)
    dstLayer.CreateField(easting_field) 
    
    northing_field = ogr.FieldDefn("NORTHING", ogr.OFTString)
    northing_field.SetWidth(12)
    dstLayer.CreateField(northing_field) 
    
    maxwidth_field = ogr.FieldDefn("MAXWIDTH", ogr.OFTInteger)
    maxwidth_field.SetWidth(8)
    dstLayer.CreateField(maxwidth_field)
    
    maxheight_field = ogr.FieldDefn("MAXHEIGHT", ogr.OFTInteger)
    maxheight_field.SetWidth(8)
    dstLayer.CreateField(maxheight_field)
    
    for FPN in glob.glob(os.path.join(FP, '*.zip')):
        if 'arctic' in FPN.lower():
            continue
        FN = os.path.split(FPN)[1]
        zipFPN = os.path.join(FP,FN)

        zipFP = os.path.splitext(zipFPN)[0]
        print (zipFP)
        '''
        zip_ref = zipfile.ZipFile(zipFPN, 'r')
        zip_ref.extractall(zipFP)
        zip_ref.close()
        '''
        for root, dirs, files in os.walk(zipFP):
            shpFPN = ''
            for name in files:
                if name.endswith("shp"):
                    shpFPN = os.path.join(root, name)

            if os.path.isfile(shpFPN):
                srcDS = driver.Open(shpFPN, 0)
                srcLayer = srcDS.GetLayer()
                spatialRef = srcLayer.GetSpatialRef()
                src_srs = mj_gis.MjProj()
                src_srs.SetProj(spatialRef)
                mgrsL = []
                layerDefinition = srcLayer.GetLayerDefn()
                fieldNameL =[]
                for i in range(layerDefinition.GetFieldCount()):
                    fieldNameL.append( layerDefinition.GetFieldDefn(i).GetName() )
                for feature in srcLayer:
                    if 'MGRS' in fieldNameL:
                        MGRS = feature.GetField("MGRS")
                    elif ('GZD' in fieldNameL and '100kmSQ_ID' in fieldNameL):
                        MGRS = '%(gzd)s%(EN)s' %{'gzd':feature.GetField("GZD"), 'EN':feature.GetField("100kmSQ_ID")}
                    else:
                        print (fieldNameL)
                        exit()
                    srcGeom = mj_gis.Geometry()
                    srcGeom.GeomFromFeature(feature)
                    
                    if srcGeom.shapelyGeom.geom_type == 'MultiPolygon': 
                        #Special routine converting MultiPolygon to Polygon
                        polygeom = mj_gis.Geometry()
                        polygeom.MultPolyToSinglePoly(srcGeom.shapelyGeom)
                        edgeL = list(polygeom.shapelyGeom.exterior.coords)
                    else:
                        edgeL = list(srcGeom.shapelyGeom.exterior.coords)
                    edgeGeom = mj_gis.Geometry()
                    edgeGeom.PointsToPolygonGeom(edgeL)
                    areakm2 = int(round(edgeGeom.shapelyGeom.area/1000000))
                    
                    #Retrieve the maximum width and the maximum height
                    xL = [p[0] for p in edgeL]
                    yL = [p[1] for p in edgeL]
                    maxwidth = int(round( max(xL)-min(xL) ) )
                    maxheight = int(round( max(yL)-min(yL) ) )

                    dstGeom = src_srs.ReprojectGeom(edgeGeom, dst_srs)

                    if 'MGRS' in fieldNameL:
                        MGRS = feature.GetField("MGRS")
                    elif ('GZD' in fieldNameL and '100kmSQ_ID' in fieldNameL):
                        MGRS = '%(gzd)s%(EN)s' %{'gzd':feature.GetField("GZD"), 'EN':feature.GetField("100kmSQ_ID")}
                    else:
                        print (fieldNameL)
                        exit()
                        
                    if 'EASTING' in fieldNameL:
                        easting = feature.GetField("EASTING")
                    elif 'Easting' in fieldNameL:
                        easting = feature.GetField("Easting")
                    elif 'EASTINT' in fieldNameL:
                        easting = feature.GetField("EASTINT")
                    elif 'E3' in fieldNameL:
                        easting = feature.GetField("E3")
                    else:
                        print (fieldNameL)
                        exit()
   
                    if MGRS in mgrsL:
                        #continue
                        pass
                    else:
                        boundsL = list(dstGeom.shapelyGeom.exterior.coords) 
                        
                        fixBoundsL = []
                        #
                        #Force the edges longitude = -180 and 180
                        if MGRS[0:2] == '01':
                            for bounds in boundsL:
                                if round( bounds[0],4 ) == -180.0000:
                                    #print(bounds )
                                    fixBoundsL.append( (-180.0,bounds[1]))
                                elif int(round(bounds[0])) == 180:
                                    fixBoundsL.append( (-180.0,bounds[1]))
                                else:
                                    fixBoundsL.append(bounds)
                        elif MGRS[0:2] == '60':
                            for bounds in boundsL:
                                if round( bounds[0],4 ) == 180.0000:
                                    #print(bounds )
                                    fixBoundsL.append( (180.0,bounds[1]))
                                elif int(round(bounds[0])) == -180:
                                    fixBoundsL.append( (180.0,bounds[1]))
                                else:
                                    fixBoundsL.append(bounds)
                        else:
                            fixBoundsL = boundsL
                        #Set the geometry                              
                        dstGeom = mj_gis.Geometry()
                        dstGeom.PointsToPolygonGeom(fixBoundsL)    
                        #Create the feature
                        featureDefn = dstLayer.GetLayerDefn()
                        transfeat = ogr.Feature(featureDefn)
                        
                        transfeat.SetGeometry(dstGeom.ogrGeom)
                        
                        transfeat.SetField("MGRS", MGRS)
                        
                        transfeat.SetField("GZD", feature.GetField("GZD"))
                        
                        transfeat.SetField("AREAKM2", areakm2)
                        
                        transfeat.SetField("EASTING", easting)
                        
                        transfeat.SetField("NORTHING", feature.GetField("NORTHING"))
                        
                        transfeat.SetField("MAXWIDTH", maxwidth)
                        
                        transfeat.SetField("MAXHEIGHT", maxheight)
                        
                        dstLayer.CreateFeature(transfeat)
                        transfeat = None
                        mgrsL.append(MGRS)
                # Save and close DataSource
                srcDS = None
    # Save and close DataSource
    dstDS = None
    
def MGRSShapeWithin(srcFPN, dstFPN, clipFPN, centroid):
    #Open the clip filekal
    driver = ogr.GetDriverByName("ESRI Shapefile")
    clipDS = driver.Open(clipFPN, 0)
    clipLayer = clipDS.GetLayer()
    clipFeature = clipLayer.GetNextFeature()
    clipGeom = mj_gis.Geometry()
    clipGeom.GeomFromFeature(clipFeature)
    
    if os.path.exists(dstFPN):
        driver.DeleteDataSource(dstFPN)
    dstDS = driver.CreateDataSource(dstFPN)
    
    srcDS = driver.Open(srcFPN, 0)
    srcLayer = srcDS.GetLayer()
    spatialRef = srcLayer.GetSpatialRef()
    if centroid:
        dstLayer = dstDS.CreateLayer('MGRS', spatialRef, geom_type=ogr.wkbPoint)
    else:
        dstLayer = dstDS.CreateLayer('MGRS', spatialRef, geom_type=ogr.wkbPolygon)
    field_name = ogr.FieldDefn("MGRS", ogr.OFTString)
    field_name.SetWidth(8)
    dstLayer.CreateField(field_name)    
        
    for feature in srcLayer:
        srcGeom = mj_gis.Geometry()
        srcGeom.GeomFromFeature(feature)
        #Check if overlap or inside
        if centroid:
            if srcGeom.ShapelyWithin(clipGeom):
                print ('copy feature',feature.GetField("MGRS"))
                # Create the feature in the layer (geojson)
                dstLayer.CreateFeature(feature)
                # Destroy the feature to free resources
                feature.Destroy()   
            else:
                print ('skip feature',feature.GetField("MGRS"))
        else:
            centroidGeom = mj_gis.Geometry()
            centroidGeom.SetShapelyGeom(srcGeom.shapelyGeom.centroid)
            if centroidGeom.ShapelyWithin(clipGeom):
                print ('copy feature',feature.GetField("MGRS"))
                # Create the feature in the layer (geojson)
                dstLayer.CreateFeature(feature)
                # Destroy the feature to free resources
                feature.Destroy()   
            else:
                print ('skip feature',feature.GetField("MGRS"))

    dstDS.Destroy()
    clipDS.Destroy()
    srcDS.Destroy()
                
if __name__ == "__main__":
    '''
    FP ='/Volumes/karttur2tb/ANCILIMPORT/NGA/MGRS_100kmSQ_ID'
    #FP ='/Volumes/karttur2tb/ANCILIMPORT/NGA/test'
    dstFN = 'mgrs_100km_polys.shp'
    MGRSorginalNGAtilesV2(FP,dstFN)
    '''
    
    
    srcFPN = '/Volumes/karttur2tb/ANCILIMPORT/NGA/MGRS_100kmSQ_ID/mgrs_100km_polys.shp'
    dstFPN = '/Volumes/karttur2tb/ANCILIMPORT/NGA/MGRS_100kmSQ_ID/NGA_mgrs_100km_polys.shp'
    clipFPN = '/Volumes/karttur2tb/ANCILIMPORT/ESA/sentinel_tiles/sentinel_2a.shp'
    
    '''
    srcFPN = '/Volumes/karttur2tb/ANCILIMPORT/ESA/sentinel_tiles/sentinel_index_polygons.shp'
    dstFPN = '/Volumes/karttur2tb/ANCILIMPORT/ESA/sentinel_tiles/sentinel_mgrs_polys.shp'
    clipFPN = '/Volumes/karttur2tb/ANCILIMPORT/ESA/sentinel_tiles/sentinel_2a.shp'
    '''
    centroid = False
    #MGRSShapeWithin identifies all mgrs tiles for which Sentinel products are created
    MGRSShapeWithin(srcFPN,dstFPN,clipFPN,centroid)
    
    
    
        
        
        

        