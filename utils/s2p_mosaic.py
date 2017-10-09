#!/usr/bin/env python
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse, json, os
import sys

# This is needed to import from a sibling folder
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import s2p
from s2plib import common

def vrt_body_source(fname,band,data_type,x_size,y_size,block_x_size,block_y_size,src_x,src_y,src_w,src_h,dst_x,dst_y,dst_w,dst_h):
    """
    Generate a source section in vrt body.
    
    Args:
        fname: Relative path to the source image
        band: index of the band to use as source
        data_type: type of the source
        x_size, y_size: size of the source raster
        block_x_size, block_y_size, size of the blocks in source raster
        src_x, src_y, src_w, src_h: source window (cropped from source image)
        dst_x, dst_y, dst_w, dst_h: destination window (where crop will be pasted)
    """
    
    body ='\t\t<SimpleSource>\n'
    body+='\t\t\t<SourceFileName relativeToVRT=\'1\'>{}</SourceFileName>\n'.format(fname)
    body+='\t\t\t<SourceProperties RasterXSize=\'{}\' RasterYSize=\'{}\' DataType=\'{}\' BlockXSize=\'{}\' BlockYSize=\'{}\' />\n'.format(x_size,y_size,data_type,block_x_size,block_y_size)
    body+='\t\t\t<SourceBand>{}</SourceBand>\n'.format(band)
    body+='\t\t\t<SrcRect xOff=\'{}\' yOff=\'{}\''.format(src_x, src_y)
    body+='xSize=\'{}\' ySize=\'{}\'/>\n'.format(src_w, src_h)
    body+='\t\t\t<DstRect xOff=\'{}\' yOff=\'{}\''.format(dst_x, dst_y)
    body+='xSize=\'{}\' ySize=\'{}\'/>\n'.format(dst_w, dst_h)
    body+='\t\t</SimpleSource>\n'

    return body

def vrt_header(w,h):
    """
    Generate vrt general header

    Args:
    w,h: size of the corresponding raster
    """
    header = '<VRTDataset rasterXSize=\"{}\" rasterYSize=\"{}\">\n'.format(w,h)

    return header

def vrt_band_header(band_index=1,dataType='Float32'):
    """
    Generate vrt header for a band
    
    Args:
        band_index: index of the band to add
        dataType: Type of the raster (default is Float32)
    
    """
    header= '\t<VRTRasterBand dataType=\"{}\" band=\"{}\">\n'.format(dataType,band_index)
    header+= '\t\t<ColorInterp>Gray</ColorInterp>\n'

    return header

def vrt_footer():
    """
    Generate vrt footer
    """
    footer= '</VRTDataset>\n'

    return footer

def vrt_band_footer():
    """
    Generate vrt footer
    """
    footer = '\t</VRTRasterBand>\n'

    return footer


def get_img_properties(img):
    """
    Read an image with gdal and returns its properties
    
    Args:
       img: Path to the image file to read

    Returns:
       a dictionary of prorperties, with the following keys: nb_bands, data_type, x_size, y_size, block_x_size, block_y_size
    """
    from osgeo import gdal
    gdal.UseExceptions()
    print("Reading image properties from: {}".format(img))
    ds = gdal.Open(img, gdal.GA_ReadOnly)

    properties = {}
    
    properties['nb_bands'] = ds.RasterCount
    
    first_band = ds.GetRasterBand(1)
    
    properties['data_type'] = gdal.GetDataTypeName(first_band.DataType)
    properties['x_size'] = first_band.XSize
    properties['y_size'] = first_band.YSize

    (bx,by) = first_band.GetBlockSize()

    properties['block_x_size'] = bx
    properties['block_y_size'] = by

    return properties

def global_extent(tiles):
    """
    Compute the global raster extent from a list of tiles
    Args:
        tiles: list of config files loaded from json files
    Returns:
         (min_x,max_x,min_y,max_y) tuple
    """
    min_x = None
    max_x = None
    min_y = None
    max_y = None

    # First loop is to compute global extent 
    for tile in tiles:
        with open(tile,'r') as f:

            tile_cfg = json.load(f)

            x = tile_cfg['roi']['x']
            y = tile_cfg['roi']['y']
            w = tile_cfg['roi']['w']
            h = tile_cfg['roi']['h']
            
            if min_x is None or x < min_x:
                min_x = x
            if min_y is None or y < min_y:
                min_y = y
            if max_x is None or x + w > max_x:
                max_x = x + w
            if max_y is None or y + h > max_y:
                max_y = y + h
                
    return(min_x,max_x,min_y,max_y)


def write_main_vrt_plain(tiles,sub_img,vrt_name,min_x,max_x,min_y,max_y,properties):
    """
    Write the main vrt file
    
    Args:
        vrt_row: The vrt files dictionnary from write_row_vrts()
        vrt_name: The output vrt_name
        min_x,max_x,min_y,max_y: Extent of the raster 
    """
    vrt_basename = os.path.basename(vrt_name)
    vrt_dirname = os.path.dirname(vrt_name)
    
    with open(vrt_name,'w') as main_vrt_file:
        
        main_vrt_file.write(vrt_header(max_x-min_x,max_y-min_y))

        for band in range(properties['nb_bands']):
            main_vrt_file.write(vrt_band_header(band+1,properties['data_type']))
                                
            # Do not use items()/iteritems() here because of python 2 and 3 compat
            for tile in tiles:
                with open(tile,'r') as f:
                
                    tile_cfg = json.load(f)

                    x = tile_cfg['roi']['x']
                    y = tile_cfg['roi']['y']
                    w = tile_cfg['roi']['w']
                    h = tile_cfg['roi']['h']

                    tile_dir = os.path.dirname(tile)
                    tile_sub_img = os.path.join(tile_dir,sub_img)
                
                    vrt_body_src=vrt_body_source(tile_sub_img,band+1,properties['data_type'],
                                                 properties['x_size'],properties['y_size'],properties['block_x_size'],properties['block_y_size'],0,0,w,
                                                 h,x,y,w,h)
            
                    main_vrt_file.write(vrt_body_src)
            main_vrt_file.write(vrt_band_footer())
        main_vrt_file.write(vrt_footer())
            
def main(tiles_file,outfile,sub_img):

    outfile_basename = os.path.basename(outfile)
    outfile_dirname  = os.path.dirname(outfile)
    
    output_format = outfile_basename[-3:]

    print('Output format is '+output_format)

    # If output format is tif, we need to generate a temporary vrt
    # with the same name
    vrt_basename = outfile_basename
    
    if output_format == 'tif':
        vrt_basename = vrt_basename[:-3]+'vrt'
    elif output_format !='vrt':
        print('Error: only vrt or tif extension is allowed for output image.')
        return
    
    vrt_name = os.path.join(outfile_dirname,vrt_basename)
    
    # Read the tiles file
    tiles = s2p.read_tiles(tiles_file)

    print(str(len(tiles))+' tiles found')

    if len(tiles) == 0:
        print("Tiles file is empty.")
        sys.exit(1)

    first_img = os.path.join(os.path.dirname(tiles[0]),sub_img)

    properties = get_img_properties(first_img)

    print("Tiles properties: {}".format(properties))

    # Compute the global extent of the output image
    (min_x,max_x,min_y,max_y) = global_extent(tiles)
    
    print('Global extent: [%i,%i]x[%i,%i]'%(min_x,max_x,min_y,max_y))

    # Now, write all row vrts
    #print("Writing row vrt files "+vrt_basename)
    #vrt_row = write_row_vrts(tiles,sub_img,vrt_basename,min_x,max_x)
    
    # Finally, write main vrt
    print('Writing '+vrt_name)
    write_main_vrt_plain(tiles,sub_img,vrt_name,min_x,max_x,min_y,max_y,properties)

    # If Output format is tif, convert vrt file to tif
    if output_format == 'tif':
        print('Converting vrt to tif ...')
        common.run(('gdal_translate -ot %s -co TILED=YES -co'
                    ' BIGTIFF=IF_NEEDED %s %s'
                    %(properties['data_type'],common.shellquote(vrt_name),common.shellquote(outfile))))

        print('Removing temporary vrt files')

        try:
            os.remove(vrt_name)
        except OSError:
            pass
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=('S2P: mosaic tool'))
    
    parser.add_argument('tiles',metavar='tiles.txt',
                        help=('path to the tiles.txt file'))
    parser.add_argument('outfile',metavar='out.tif',
                        help=('path to the output file.'
                              ' File extension can be .tif or .vrt'))
    parser.add_argument('sub_img',metavar='pair_1/height_map.tif',
                        help=('path to the sub-image to mosaic.'
                              ' Can be (but not limited to) height_map.tif,'
                              ' pair_n/height_map.tif, pair_n/rpc_err.tif,'
                              ' cloud_water_image_domain_mask.png.'
                              ' Note that rectified_* files CAN NOT be used.'))
    args = parser.parse_args()

    main(args.tiles,args.outfile,args.sub_img)
