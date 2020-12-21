"""
Prepare LAC output product.

Author: Martino Ferrari (CS GROUP)
Email: martino.ferrari@c-s.fr
"""
import argparse
import os
import sys

from lxml import etree

__COMMON_BANDS__ = [('Rayleigh865',), ('RGB', 'jpg'), ('SolarAzimuth',), \
    ('SolarZenith',), ('ViewingAzimuth',), ('ViewingZenith',)]
__AERO_BANDS__ = [('Aerosol865',), ('Cloud60',), ('CloudMask', 'jpg')]

__COMMON_AUX__ = ['MeanO3.txt']
__AERO_AUX__ = ['Correc.txt', 'CloudAndEps_v3.txt']

#             <PRODUCT_URI>S2A_MSIL1C_20190712T085601_N0208_R007_T35TNH_20190712T110428.SAFE</PRODUCT_URI>

class S2Product:
    def __init__(self):
        self.name = None


def __args__():
    parser = argparse.ArgumentParser(description='Prepare LAC product.')
    parser.add_argument('-input', required=True, help='S2 l1c input product folder')
    parser.add_argument('-output', required=True, help='LAC output prodcut folder')
    parser.add_argument('-bands', required=True, help='bands to be used')
    parser.add_argument('-aero', choices=['y', 'n'], help='enable or disable areo')
    
    return parser.parse_args()

def __prepare__():
    args = __args__()
    prepare(args.input, args.output, args.bands, args.aero)


def prepare(input, output, bands, aero):
    """
    Prepare data structure and metadata for LAC product.
    """
    s2_metadata = os.path.join(input, 'MTD_MSIL1C.xml')
    
    if not os.path.exists(input) or not os.path.exists(s2_metadata):
        print(f'[ERROR] The input S2 prodcuts `{input}` does not exists')
        sys.exit(1)
    if not os.path.exists(output):
        print(f'[WARNING] The ouptut LAC directory `{output} does not exits')
        os.mkdir(output)
    elif len(os.listdir(output)) > 0:
        print(f'[ERROR] The output LAC directory `{output}` is not empty')
        sys.exit(1)

    gr_path = os.path.join(input, 'GRANULE')
    gr_path = os.listdir(os.path.join(gr_path, os.listdir(gr_path)[0], 'IMG_DATA'))[0]
    gr_name = gr_path.split('_')
    gr_name = '_'.join(gr_name[0:2])
    gr_date = gr_name.split('_')[1][:8]

    band_list = [b.upper() for b in bands.split('_')]

    lac_root = etree.Element('LAC-Product')
    lac_prod_info = etree.SubElement(lac_root, 'Product_Info')
    lac_proc = etree.SubElement(lac_prod_info, 'Processing')
    ptype = etree.SubElement(lac_prod_info, 'PRODUCT_FORMAT')
    ptype.text = 'LAC'

    pbands = etree.SubElement(lac_proc, "BANDS")
    pbands.text = bands
    pareo = etree.SubElement(lac_proc, "AREO")
    pareo.text = aero

    p_org = etree.SubElement(lac_prod_info,  'Product_Organisation')
    p_gran = etree.SubElement(p_org, 'Granule')

    for band in band_list:
        b = etree.SubElement(p_gran, 'IMAGE_FILE')
        b.text = f'GRANULE/C_{gr_name}_{band}.TIF'

    b = etree.SubElement(p_gran, 'IMAGE_FILE')
    b.text = f'GRANULE/O3-{gr_date}.TIF'
        
    for band in __COMMON_BANDS__:
        ext = 'TIF' if len(band) == 1 else band[1]
        b = etree.SubElement(p_gran, 'IMAGE_FILE')
        b.text = f'GRANULE/{band[0]}_{gr_name}.{ext}'

    
    p_aux = etree.SubElement(p_org, 'Aux_Data')
    for aux in __COMMON_AUX__:
        a = etree.SubElement(p_aux, 'AUX_FILE')
        a.text = f'AUX_DATA/{aux}'

    if aero == 'y':
        for band in __AERO_BANDS__:
            ext = 'TIF' if len(band) == 1 else band[1]
            b = etree.SubElement(p_gran, 'IMAGE_FILE')
            b.text = f'GRANULE/{band[0]}_{gr_name}.{ext}'
        for aux in __AERO_AUX__:
            a = etree.SubElement(p_aux, 'AUX_FILE')
            a.text = f'AUX_DATA/{aux}'
    else:
        a = etree.SubElement(p_aux, 'AUX_FILE')
        a.text = f'AUX_DATA/O3-{gr_date}.TIF.aux.xml'
        a = etree.SubElement(p_aux, 'AUX_FILE')
        a.text = f'AUX_DATA/Rayleigh865_{gr_name}.TIF.aux.xml'


    root = etree.parse(s2_metadata)
    for e in root.xpath("//Product_Info/Datatake"):
        lac_prod_info.append(e)
    for e in root.xpath("//Product_Info/PRODUCT_START_TIME"):
        lac_prod_info.append(e)
    for e in root.xpath("//Product_Info/PRODUCT_STOP_TIME"):
        lac_prod_info.append(e)
    for e in root.xpath("//Product_Info/PRODUCT_URI"):
        new_el = etree.SubElement(lac_prod_info, 'PRODUCT_URI')
        new_el.text = e.text[:-5] + '_LAC'
    # etree.indent(lac_root, space="  ")
    # print(etree.tostring(lac_root, encoding='utf-8').decode())
    # os.mkdir(output)
    with open(os.path.join(output, 'MTD_LAC.xml'), 'wb') as file:
        file.write(etree.tostring(lac_root, encoding='utf-8'))
    os.mkdir(os.path.join(output, 'GRANULE'))
    os.mkdir(os.path.join(output, 'AUX_DATA'))
    os.mkdir(os.path.join(output, 'TMP'))


if __name__ == "__main__":
    __prepare__()
