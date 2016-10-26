# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>


import os
import shutil
import numpy as np
import multiprocessing

from python import tile_composer
from python import common
from config import cfg


def write_vrt_files(tiles_full_info):
    """
    Merges pieces of data into single VRT files : height map comprising the N pairs, height map for each signle pair, and err_rpc

    Args:
         tiles_full_info: a list of tile_info dictionaries
    """
    # VRT file : height map (N pairs)

    #-------------- tileComposerInfo --------------
    # tileComposerInfo : a simple list that gives the size of the tiles (after having removed the overlapping areas),
    # regarding their positions inside the ROI --> ULw , ULh, Lw , Lh, Uw ,
    # Uh, Mw , Mh
    x0, y0, tw, th = tiles_full_info[0]['coordinates']
    x, y, w, h = tiles_full_info[0]['roi_coordinates']
    ov = tiles_full_info[0]['overlap']
    nb_pairs = tiles_full_info[0]['number_of_pairs']

    ULw, ULh = tw - ov / 2, th - ov / 2  # Size of Corner tiles (UL UR BL BR)
    Lw, Lh = tw - ov / 2, th - ov  # Size of Left/Right tile (L R)
    Uw, Uh = tw - ov, th - ov / 2  # Size of Upper/Bottom tile (U B)
    Mw, Mh = tw - ov, th - ov  # Size of Middle tile

    rangey = np.arange(y, y + h - ov, th - ov)
    rangex = np.arange(x, x + w - ov, tw - ov)
    rowmin, rowmax = rangey[0], rangey[-1]
    colmin, colmax = rangex[0], rangex[-1]

    # Tile composer info (bis)
    imax = len(rangey) - 1
    jmax = len(rangex) - 1
    # height of the final height map after removing margins
    fh = 2 * ULh + (imax - 1) * Mh
    # width of the final height map after removing margins
    fw = 2 * ULw + (jmax - 1) * Mh
    #-------------- tileComposerInfo (end) --------------

    tileSizesAndPositions = {}
    for tile_info in tiles_full_info:
        col, row, tw, th = tile_info['coordinates']
        pos = tile_info['position_type']
        i, j = tile_info['index_in_roi']

        dicoPos = {}
        dicoPos['M'] = [ULw + (j - 1) * Mw, ULh + (i - 1) * Mh, Mw, Mh]
        dicoPos['UL'] = [0, 0, ULw, ULh]
        dicoPos['UR'] = [fw - ULw, 0, ULw, ULh]
        dicoPos['BR'] = [fw - ULw, fh - ULh, ULw, ULh]
        dicoPos['BL'] = [0, fh - ULh, ULw, ULh]
        dicoPos['L'] = [0, ULh + (i - 1) * Lh, Lw, Lh]
        dicoPos['R'] = [fw - ULw, ULh + (i - 1) * Lh, Lw, Lh]
        dicoPos['U'] = [ULw + (j - 1) * Uw, 0, Uw, Uh]
        dicoPos['B'] = [ULw + (j - 1) * Uw, fh - ULh, Uw, Uh]
        dicoPos['Single'] = [0, 0, tw, th]

        tile_reldir = 'tile_%d_%d_row_%d/col_%d/' % (tw, th, row, col)

        tileSizesAndPositions[tile_reldir] = dicoPos[pos]

    z = cfg['subsampling_factor']
    
    # VRT file : height map 
    height_map_path = os.path.join(cfg['out_dir'], 'height_map.vrt')
    if not (os.path.isfile(height_map_path) and cfg['skip_existing']):
        tile_composer.mosaic_stitch( height_map_path,
                               tileSizesAndPositions, 'height_map_crop.tif', fw, fh, 1, z)
                               
    # VRT file : rpc_err_rms_allsights
    rpc_err_path = os.path.join(cfg['out_dir'], 'rpc_err_rms_allsights.vrt')
    if not (os.path.isfile(rpc_err_path) and cfg['skip_existing']):
        tile_composer.mosaic_stitch( rpc_err_path,
                               tileSizesAndPositions, 'rpc_err_rms_allsights_crop.tif', fw, fh, 1, z)
    
    if cfg['full_vrt']:
        # VRT file : nb_sights
        nb_sights_path = os.path.join(cfg['out_dir'], 'nb_sights.vrt')
        if not (os.path.isfile(nb_sights_path) and cfg['skip_existing']):
            tile_composer.mosaic_stitch( nb_sights_path,
                               tileSizesAndPositions, 'nb_sights_crop.tif', fw, fh, 1, z)
        
        # selected_sight_i
        # rpc_err_norm_sight_i 
        # rpc_err_vec_sight_i   
        # rpc_err_rpjvec_sight_i                
        for img_id in xrange(1,len(cfg['images'])+1): 
            #selected sights
            selected_sighti = 'selected_sight_%d.vrt' % img_id
            selected_sighti_crop = 'selected_sight_%d_crop.tif' % img_id
            selected_sighti_path = os.path.join(cfg['out_dir'], selected_sighti)
            if not (os.path.isfile(selected_sighti_path) and cfg['skip_existing']):
                tile_composer.mosaic_stitch( selected_sighti_path,
                                   tileSizesAndPositions, selected_sighti_crop, fw, fh, 1, z)
             
            # err by sight                   
            rpc_err_sighti = 'rpc_err_norm_sight_%d.vrt' % img_id
            rpc_err_sighti_crop = 'rpc_err_norm_sight_%d_crop.tif' % img_id
            rpc_err_sighti_path = os.path.join(cfg['out_dir'], rpc_err_sighti)
            if not (os.path.isfile(rpc_err_sighti_path) and cfg['skip_existing']):
                tile_composer.mosaic_stitch( rpc_err_sighti_path,
                                   tileSizesAndPositions, rpc_err_sighti_crop, fw, fh, 1, z)
            # err vectors by sight
            rpc_err_veci = 'rpc_err_vec_sight_%d.vrt' % img_id
            rpc_err_veci_crop = 'rpc_err_vec_sight_%d_crop.tif' % img_id
            rpc_err_veci_path = os.path.join(cfg['out_dir'], rpc_err_veci)
            if not (os.path.isfile(rpc_err_veci_path) and cfg['skip_existing']):
                tile_composer.mosaic_stitch( rpc_err_veci_path,
                                   tileSizesAndPositions, rpc_err_veci_crop, fw, fh, 3, z)
                                   
            # reprojected err vectors by sight
            rpc_err_vec_rpji = 'rpc_err_rpjvec_sight_%d.vrt' % img_id
            rpc_err_vec_rpji_crop = 'rpc_err_rpjvec_sight_%d_crop.tif' % img_id
            rpc_err_vec_rpji_path = os.path.join(cfg['out_dir'], rpc_err_vec_rpji)
            if not (os.path.isfile(rpc_err_vec_rpji_path) and cfg['skip_existing']):
                tile_composer.mosaic_stitch( rpc_err_vec_rpji_path,
                                   tileSizesAndPositions, rpc_err_vec_rpji_crop, fw, fh, 3, z)
                                   
        # 2D disparities (if originaly computed in epipolar geometry)                     
        for pair_id in xrange(1,nb_pairs+1):                     
            disp2Di = 'disp2D_pair%d.vrt' % pair_id
            disp2Di_crop = 'pair_%d/disp2D_crop.tif' % pair_id
            disp2Di_path = os.path.join(cfg['out_dir'], disp2Di)
            if not (os.path.isfile(disp2Di_path) and cfg['skip_existing']):
                tile_composer.mosaic_stitch( disp2Di_path,
                                   tileSizesAndPositions, disp2Di_crop, fw, fh, 3, z)
                                   
        # Cleaning
        if(cfg['clean_intermediate'] and cfg['vrt_to_tiff']):
            
            for tile_info in tiles_full_info:
                tile_dir = tile_info['directory']
                nb_pairs = tile_info['number_of_pairs']
                
                for pair_id in xrange(1, nb_pairs + 1):
                    pair_dir = os.path.join(tile_dir,'pair_%d' %(pair_id+1))
                    shutil.rmtree( pair_dir )
                    
                # height maps
                common.remove_if_exists(os.path.join(tile_dir,'height_map_crop.tif'))
                # rpc_err_all
                common.remove_if_exists(os.path.join(tile_dir,'rpc_err_rms_allsights_crop.tif'))
                
                # other tif
                if cfg['full_vrt']:
                    for i in xrange(1,len(cfg['images'])+1):
                        common.remove_if_exists(os.path.join(tile_dir,'rpc_err_norm_sight_%d_crop.tif' % i ))
                        common.remove_if_exists(os.path.join(tile_dir,'rpc_err_rpjvec_sight_%d_crop.tif' % i ))
                        common.remove_if_exists(os.path.join(tile_dir,'rpc_err_vec_sight_%d_crop.tif' % i ))
                        common.remove_if_exists(os.path.join(tile_dir,'selected_sight_%d_crop.tif' % i ))
                    common.remove_if_exists(os.path.join(tile_dir,'nb_sights_crop.tif'))
                
            # height maps vrt
            common.remove_if_exists(os.path.join(cfg['out_dir'],'height_map.vrt'))
            # rpc_err_all vrt
            common.remove_if_exists(os.path.join(cfg['out_dir'],'rpc_err_rms_allsights.vrt'))
            
            # other vrt
            if cfg['full_vrt']:
                for i in xrange(1,len(cfg['images'])+1):
                    common.remove_if_exists(os.path.join(cfg['out_dir'],'rpc_err_norm_sight_%d.vrt' % i ))
                    common.remove_if_exists(os.path.join(cfg['out_dir'],'rpc_err_rpjvec_sight_%d.vrt' % i ))
                    common.remove_if_exists(os.path.join(cfg['out_dir'],'rpc_err_vec_sight_%d.vrt' % i ))
                    common.remove_if_exists(os.path.join(cfg['out_dir'],'selected_sight_%d.vrt' % i ))
                common.remove_if_exists(os.path.join(cfg['out_dir'],'nb_sights.vrt'))
                
                for pair_id in xrange(1, nb_pairs + 1):
                    disp2D_path = os.path.join(cfg['out_dir'],'disp2D_pair%d.vrt' %(pair_id+1))
                    common.remove_if_exists(disp2D_path)
            

def write_dsm():
    """
    Writes the DSM, from the ply files given by each tile.
    """
    dsm_pieces = os.path.join(cfg['out_dir'], 'dsm/dsm_*')
    final_dsm = os.path.join(cfg['out_dir'], 'dsm.vrt')
    if not (os.path.isfile(final_dsm) and cfg['skip_existing']):
        common.run("gdalbuildvrt %s %s" % (final_dsm, dsm_pieces))

    if cfg['vrt_to_tiff']:
        final_dsm_tif = os.path.splitext(final_dsm)[0]+".tif"
        if not (os.path.isfile(final_dsm_tif) and cfg['skip_existing']):
            common.run('gdal_translate %s %s' % (final_dsm, final_dsm_tif))
        if cfg['clean_intermediate']:
            shutil.rmtree(os.path.join(cfg['out_dir'],'dsm'))
            common.remove_if_exists(final_dsm);
                          

def lidar_preprocessor(output, input_plys):
    """
    Compute a multi-scale representation of a large point cloud.

    The output file can be viewed with LidarPreprocessor. This is useful for
    huge point clouds. The input is a list of ply files.

    Args:
        output: path to the output folder
        input_plys: list of paths to ply files
    """
    tmp = cfg['temporary_dir']
    nthreads = multiprocessing.cpu_count()
    plys = ' '.join(input_plys)
    common.run("LidarPreprocessor -to %s/LidarO -tp %s/LidarP -nt %d %s -o %s" % (
        tmp, tmp, nthreads, plys, output))
