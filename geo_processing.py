# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 10:12:46 2022

@author: ArthurDiniz
"""

from osgeo import gdal, osr, ogr
import geopandas as gpd
import pandas as pd
import glob
import numpy as np
import os

def raster2vector(in_path,out_path):
    src_ds = gdal.Open( in_path )

    srcband = src_ds.GetRasterBand(1)
    dst_layername = 'oilpalm_HarvestedAreaHectares'
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource( out_path )

    sp_ref = osr.SpatialReference()
    sp_ref.SetFromUserInput('EPSG:4326')

    dst_layer = dst_ds.CreateLayer(dst_layername, srs = sp_ref )

    fld = ogr.FieldDefn("DN", ogr.OFTInteger)
    dst_layer.CreateField(fld)
    dst_field = dst_layer.GetLayerDefn().GetFieldIndex("DN")

    gdal.Polygonize(srcband, None, dst_layer, dst_field, [], callback=None )

    src_ds = None
    dst_layername = None
    drv = None
    dst_ds = None
    sp_ref = None
    dst_layer = None
    fld = None
    dst_field = None

    vet_base = gpd.read_file(out_path)
    vet_base = vet_base.loc[vet_base.DN==1]
    vet_base.to_file(out_path)

def create_buffer(line_add,out_add,distance_km=1,src ='EPSG:4674'):
    gdf = gpd.read_file(line_add)
    gdf.crs = src
    distance_degree = distance_km/111.12
    gdf = gdf.to_crs(src)
    gdf.geometry = gdf.geometry.buffer(distance_degree, 6)
    gdf.to_file(out_add)  

def join_shape(shape_add,out_add):
    #shape_add = r'C:\Arthur\Niras\dados\Entrega 2\Restrições Socioambientais\teste'
    add_list = glob.glob(shape_add+'/*.shp')
    out_gdf =gpd.GeoDataFrame()
    for n in add_list:
        gdf = gpd.read_file(n)
        out_gdf=out_gdf.append(gdf)
    out_gdf.to_crs("EPSG:4326")
    out_gdf.to_file(out_add+'/merged.shp')
        

def rasterize_reference(in_shp_file_name,in_raster_file_name,out_raster_file_name,src=4326,burn_val=1):
    mat = get_tif_array(in_raster_file_name)
    gt = get_tif_gt(in_raster_file_name)
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(in_shp_file_name, 1) #1 is read/write
    
    shp = ogr.Open(in_shp_file_name)

    lyr = shp.GetLayer()
    
    mat2 = mat/mat
    mat[mat!=-9999]=0
    mat[mat==-9999]=-10000
    mat2 = mat2+mat
    mat2[mat2==1]=0
    
    create_tif(out_raster_file_name,gt,mat2)
    data = gdal.Open(out_raster_file_name,gdal.GA_Update)
    
    gdal.RasterizeLayer(data, [1], lyr, None,burn_values=[burn_val])
    data = None
    lyr = None
    shp = None
    dataSource = None
    driver = None


def get_tif_array(tif_add):
    data = gdal.Open(tif_add)
    data_array = data.ReadAsArray()
    return data_array

def get_tif_gt(tif_add):
    data = gdal.Open(tif_add)
    gt = data.GetGeoTransform()
    return gt


def create_tif(tif_add,gt,tif_matrix,no_data=-9999,srs=4326):
    Y,X = tif_matrix.shape
    
    driver = gdal.GetDriverByName("GTiff")
    outRaster = driver.Create(tif_add, X, Y, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform(gt)
    outBand = outRaster.GetRasterBand(1)
    outBand.SetNoDataValue(no_data)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(srs)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outBand.WriteArray(tif_matrix)
    outBand.FlushCache()

    outRaster = None
    outBand = None
    outRasterSRS = None
    driver = None
    
    
def zonal_statistics(shape_add,tif_add,out_add,ope='all'):
    gdf = gpd.read_file(shape_add)
    gdf = gdf.reset_index()
    gdf = gdf.rename(columns={"index":"chn_id"})
    gdf.to_file(shape_add.replace('.shp','code.shp'))
    
    
    mat = get_tif_array(tif_add)
    gt = get_tif_gt(tif_add)
    
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shape_add.replace('.shp','code.shp'), 1) #1 is read/write
    
    shp = ogr.Open(shape_add.replace('.shp','code.shp'))

    lyr = shp.GetLayer()
    
    mat2 = mat/mat
    mat[mat!=-9999]=0
    mat[mat==-9999]=-10000
    mat2 = mat2+mat
    mat2[mat2==1]=-1
    create_tif(tif_add.replace('.tif','_code.tif'),gt,mat2)
    data = gdal.Open(tif_add.replace('.tif','_code.tif'),gdal.GA_Update)
    
    gdal.RasterizeLayer(data, [1], lyr, None,options=['ATTRIBUTE=chn_id'])
    data = None
    lyr = None
    shp = None
    dataSource = None
    driver = None
    
    mat = get_tif_array(tif_add)
    mat_code = get_tif_array(tif_add.replace('.tif','_code.tif')).flatten()
    mat = mat.flatten()
    mat_code=mat_code[mat_code>=0]
    mat=mat[mat>=0]
    
    if ope == 'sum' or ope =='all':
        for n in gdf['chn_id']:
            op_array = mat[np.nonzero(mat_code == n)]
            gdf.loc[gdf.chn_id==n,'sum']=op_array.sum()
    if ope == 'mean' or ope=='all':
        for n in gdf['chn_id']:
            op_array = mat[np.nonzero(mat_code == n)]
            gdf.loc[gdf.chn_id==n,'mean']=op_array.mean()
    if ope == 'std' or ope=='all':
        for n in gdf['chn_id']:
            op_array = mat[np.nonzero(mat_code == n)]
            gdf.loc[gdf.chn_id==n,'std']=op_array.std()
    
    gdf = gdf.drop('chn_id')
    #df = gdf.drop('geometry',axis=1)
    gdf.to_file(out_add)
    
    file_path = tif_add.replace('.tif','_code.tif')
    os.remove(file_path)
    file_path = shape_add.replace('.shp','code.shp')
    file_path = file_path.replace('.shp','.*')
    file_path = glob.glob(file_path)
    for n in file_path:
        os.remove(n)
