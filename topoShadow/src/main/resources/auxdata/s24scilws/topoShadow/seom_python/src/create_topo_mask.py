# -*- coding: utf-8 -*-
'''
SEOM Project create topography mask
@author:     Jonathan Guinet (CS-SI)
@copyright:  2019 CS-SI. All rights reserved.
@contact:    jonathan.guinet@c-s.fr 
'''

import os
import time
import sys

import gdal
from gdalconst import GA_ReadOnly
import ogr

from argparse import ArgumentParser
import logging

from osgeo import osr
import numpy as np
import subprocess

import xml.etree.cElementTree as _ET

logging.basicConfig()
logger = logging.getLogger('logger')


def process_mtd(metadata):
    product_info = {}
    e = _ET.parse(metadata)
    for elt in e.iter():
        if elt.tag == "SENSING_TIME":
            product_info[elt.tag] = elt.text
        if elt.tag == "HORIZONTAL_CS_CODE":
            product_info[elt.tag] = elt.text
        if elt.tag == "ULX":
            product_info[elt.tag] = int(elt.text)            
        if elt.tag == "ULY":
            product_info[elt.tag] = int(elt.text)
        if elt.tag == "TILE_ID":
            product_info[elt.tag] = elt.text
        product_info['NROWS'] = 10980
        product_info['NCOLS'] = 10980
        product_info['XDIM'] = 10.0
        product_info['YDIM'] = -10.0                  
    return product_info

def get_granule_mtd(product_dir):
    granule_dir = os.path.join(product_dir,"GRANULE")
    dirs=os.listdir(granule_dir)
    for granule_subdir in dirs :
        granule_mtd_dir =  os.path.join(granule_dir,granule_subdir)
        for filename in os.listdir(granule_mtd_dir):
            if filename.endswith(".xml"):
                return os.path.join(granule_mtd_dir,filename)

def getDEMdir(conf):
    tag = 'dem.root.directory'
    dem_root_dir = None
    with open (conf, "r") as fileHandler:
    # Read each line in loop
        for line in fileHandler:
        # As each line (except last one) will contain new line character, so strip that
            if(line.startswith(tag)):
                dem_root_dir = line.split('=')[1][:-1]
    if(dem_root_dir is None):
        raise Exception("Error : could not found DEM root dir in rugged conf")
    #print dem_root_dir
    for root, ___, files in os.walk(dem_root_dir):
        for filename in files:
            if filename.endswith(".tif"):
                return os.path.join(root, filename)
            
    

def convert_into_DEM_CRS(conf,product_UL,product_LR):
    dem_tile_image = getDEMdir(conf)
    logging.info(dem_tile_image)
        #read image 
    DEM_dataset = gdal.Open(dem_tile_image, GA_ReadOnly)
    if DEM_dataset is None:
        msg = 'error in DEM reading'
        raise Exception(msg)

    geotransform = DEM_dataset.GetGeoTransform()
    #print 'geo ', geotransform
    DEM_dataset = None
    #round to geotansform
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    deltaX = np.floor((product_UL[0] - originX)/pixelWidth)
    UL_DEM_X = originX + deltaX * pixelWidth
    deltaY = np.floor((product_UL[1] - originY)/pixelHeight)
    UL_DEM_Y = originY + deltaY * pixelHeight   

    size_X = int(np.ceil ((product_LR[0] - UL_DEM_X)/pixelWidth))
    size_Y = int(np.ceil ((product_LR[1] - UL_DEM_Y)/pixelHeight))
    return ([UL_DEM_X,UL_DEM_Y],[pixelWidth,pixelHeight],[size_X,size_Y])

def compute_grid(rugged_config, grid_image, date,EPSG, UL, resolution,size):
    """
    """

    try:
        # Create location compute_grid
        out = subprocess.call(['java', '-Djava.library.path={}'.format(os.environ['GDAL_LIBRARY_PATH']), '-jar', os.environ['SEOM_JAVA'], 
                rugged_config,"{}".format(date), grid_image, "{}".format(UL[0]), "{}".format(UL[1]), "{}".format(resolution[0]), "{}".format(resolution[1]), "{}".format(size[0]),"{}".format(size[1]),  "{}".format(EPSG)])
    except Exception:
        raise Exception("Error : could not create location grid from Rugged")
    if not os.path.exists(grid_image):
        raise Exception("Error : could not create location grid from Rugged")
    


    
def create_shapefile(output_mask,output):
    src_ds = gdal.Open( output_mask )
    if src_ds is None:
        logging.error('Unable to open %s', output_mask)
        sys.exit(1)
    proj_data = src_ds.GetProjectionRef()
    srcband = src_ds.GetRasterBand(1)
   
    if (os.path.exists(output)):
        os.remove(output)
    logging.info("creating shapefile...")
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource( output)
    srs = osr.SpatialReference()
    srs.ImportFromWkt( proj_data )
    dst_layer = dst_ds.CreateLayer(output, srs = srs )

    dst_fieldname = 'DN'

    fd = ogr.FieldDefn( dst_fieldname, ogr.OFTInteger )
    dst_layer.CreateField( fd )
    dst_field = 0

    gdal.Polygonize( srcband, srcband, dst_layer, dst_field, [], callback = None)


def create_mask_fromgrid(input_distance,output_mask,threshold):
    
    grid_dataset = gdal.Open(input_distance, GA_ReadOnly)
    if grid_dataset is None:
        msg = 'error in grid reading'
        raise Exception(msg)

    geotransform = grid_dataset.GetGeoTransform()
    grid_band = grid_dataset.GetRasterBand(1)
    proj_data = grid_dataset.GetProjectionRef()
    sizeX = grid_dataset.RasterXSize
    sizeY = grid_dataset.RasterYSize   
    grid_data = grid_band.ReadAsArray().astype(np.float)
    grid_data_mask = grid_data > threshold

    
    if(os.path.exists(output_mask)):
        os.remove(output_mask)


    driver = gdal.GetDriverByName('GTiff')

    dataset = driver.Create(output_mask, sizeX,sizeY,1,gdal.GDT_Byte)

    if(dataset is None):
        msg = "file %s can't be created " % output_mask
        logging.error(msg)

    dataset.SetGeoTransform(geotransform)
    dataset.SetProjection(proj_data)

    dataset.GetRasterBand(1).WriteArray(grid_data_mask)
    dataset.GetRasterBand(1).SetNoDataValue(0)
        #dataset.GetRasterBand(2).WriteArray(data_y)
        #dataset.GetRasterBand(2).SetNoDataValue(nodata)        
        
    if (os.path.exists(output_mask) == False):
        msg = "file %s has not be created " % output_mask
        logging.error(msg)    
        

def create_topo_mask():
    """ mask_creation_at_DEM_resolution.
        Les options de traitement sont:
    """
    prg = "create_topo_mask"


    logger.info('')
    logger.info('------------------------------')
    logger.info(prg)
    logger.info('------------------------------')

    logger.info('')

    parser = ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                        default=False, help="Verbose mode")

    parser.add_argument("ruggedConf",help="rugged configuration file")
    parser.add_argument("product",help="product (.SAFE root dir) or L1C/L2A granule MTD (.xml)")
    parser.add_argument("--ULLR",type=int, nargs=4, help="user defined footprint [ULX ULY LRX LRY]")
    parser.add_argument("--EPSG",type=int, help="EPSG code of --ULLR")
    parser.add_argument("output", help="output mask filename (vector data in ESRI shaprefile format)")
    parser.add_argument("tmp",help="working dir")
    parser.add_argument("--threshold",type=float,default=0.1,help="detection threshold")   
    args = parser.parse_args()

    # Processing time
    t0 = time.time()

    
    # Read product (Vesuve area)
    if (args.product.endswith(".xml")):
        metadata = args.product
    else:
        metadata = get_granule_mtd(args.product) 
    
    product_infos = process_mtd(metadata)


    product_width = product_infos['NCOLS']
    product_height = product_infos['NROWS']
    
    ULX_carto = product_infos['ULX']
    ULY_carto = product_infos['ULY']
    pixel_size_x = product_infos['XDIM']
    pixel_size_y = product_infos['YDIM']  
    LRX_carto = ULX_carto + pixel_size_x *product_width
    LRY_carto = ULY_carto + pixel_size_y *product_height

    epsg_code = int(product_infos['HORIZONTAL_CS_CODE'].split(':')[1])
    if args.ULLR is not None:
        epsg_code = args.EPSG
        ULLR = args.ULLR
        ULX_carto = ULLR[0]
        ULY_carto = ULLR[1]
        LRX_carto = ULLR[2]
        LRY_carto = ULLR[3]
    

    inspatialRef = osr.SpatialReference()
    inspatialRef.ImportFromEPSG(epsg_code)
    outspatialRef = osr.SpatialReference()
    outspatialRef.ImportFromEPSG(4326)
    coordTrans = osr.CoordinateTransformation(inspatialRef, outspatialRef) 

    product_UL =  coordTrans.TransformPoint(ULX_carto, ULY_carto)
    #print ' product UL', product_UL
    product_LR =  coordTrans.TransformPoint(LRX_carto, LRY_carto)
    #print 'product LR', product_LR          
    gdal.AllRegister()
    UL_DEM,res_DEM,size_DEM = convert_into_DEM_CRS(args.ruggedConf,product_UL,product_LR)

    #print 'UL_DEM ',UL_DEM
    #print 'size_DEM', size_DEM 
       
    product_name = product_infos["TILE_ID"]
    logging.info('product_name: %s', product_name)
    reference_date = product_infos['SENSING_TIME']
    logging.info('reference_date: %s', reference_date)   
    logging.info('compute grid')
    grid_image = os.path.join(args.tmp, "gridDistance_{}.tif".format(product_name))
    
    #GRID CRS is WGS84 EPSG code is 4326
    EPSG = 4326
    compute_grid(args.ruggedConf, grid_image, reference_date,EPSG, UL_DEM, res_DEM,size_DEM)
   
    #create mask from grid
    mask = os.path.join(args.tmp, "mask_{}.tif".format(product_name))
    create_mask_fromgrid(grid_image,mask,args.threshold)
    create_shapefile(mask,args.output)

    
    logging.info('Done. Execution time: %f seconds',time.time() - t0)

if __name__ == "__main__":
        
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    create_topo_mask()
