#!/usr/bin/env python
"""Read a THREDDS catalog [XML] and return a list of URLS for files of a given type (nc, csv.gz, png, json).
"""

# Module loading
import argparse
from lxml import etree
import re
import requests
import time as clock


__author__ = "Aaron Sweeney"
__copyright__ = ""
__credits__ = ["Aaron Sweeney"]
__license__ = "GNU General Public License, version 3.0"
__version__ = "1.0.1"
__maintainer__ = "Aaron Sweeney"
__email__ = "aaron.sweeney@colorado.edu"
__status__ = "Development"

#Modified from Aaron Sweeny code
def is_valid_ftype(string):

    ftypes = ['nc', 'csv.gz', 'json', 'png']
    if string not in ftypes:
        msg = 'Invalid file type: ', string
        raise argparse.ArgumentTypeError(msg)
    return string


def read_thredds_catalog(url_catalog, ftype):
    """Read a THREDDS catalog [XML] and return a list of URLS for files of a given type (nx, csv.gz, png, json)

    Args:
      url_catalog (str): URL of top-level of THREDDS catalog
      ftype (str): File type (one of 'nc', 'csv', 'png', or 'json')

    Returns:
      list of URLs for files of given type
    """

    try:
        response = requests.get(url_catalog)
        tree_data = etree.fromstring(response.content)
    except etree.ParseError:
        print('{}: XML parse error'.format(url_catalog))
        return

    station_nodes = tree_data.findall('.//catalogRef', namespaces=tree_data.nsmap)

    print('Number of stations: {}'.format(len(station_nodes)))

    urls_station_dataset = []

    ns = {k: v for k, v in tree_data.nsmap.items() if k is not None}

    # Loop over all station catalogs:
    for node in station_nodes:

        if node.xpath('.//@xlink:href', namespaces=ns) is not None:
            url_catalog_station = url_catalog.replace('catalog.xml', node.xpath('.//@xlink:href', namespaces=ns)[0])
        else:
            print('Error: Could not get xlink:href attribute of catalogReg element from catalog')
            return

        try:
            response = requests.get(url_catalog_station)
            station_tree_data = etree.fromstring(response.content)
        except etree.ParseError:
            print('{}: Parse error'.format(url_catalog_station))
            return

        station_datasets = station_tree_data.findall('.//dataset/dataset', namespaces=station_tree_data.nsmap)

        # Look for datasets that match given file type:
        for dataset in station_datasets:

            if re.match('^[\w]*.' + ftype + '$', dataset.get('name')):
                url_station_dataset = url_catalog_station.replace('/catalog/', '/fileServer/')
                url_station_dataset = url_station_dataset.replace('catalog.xml', dataset.get('name'))
                urls_station_dataset.append(url_station_dataset)

    return urls_station_dataset


if __name__ == "__main__":

    # Let's track how long this takes.
    program_start_time = clock.time()

    parser = argparse.ArgumentParser(
        description='Read THREDDS catalog [XML] and return list of files of given type (nc, csv.gz, json, png)')
    parser.add_argument('-u', '--url', metavar='url', type=str, required=True,
                        help='URL of THREDDS catalog [XML]')
    parser.add_argument('-t', '--file_type', metavar='file_type', type=is_valid_ftype,
                        required=True, help='Desired file type (one of nc, csv.gz, json, or png)')
    args = parser.parse_args()

    # print('Args: {}'.format(args))

    base_url = args.url
    ftype = args.file_type

    url_list = read_thredds_catalog(base_url, ftype)

    if url_list is not None:

        for url in url_list:

            print('{}'.format(url))
    else:

        print('No {} files found at {}'.format(ftype, base_url))

    # Capture the program's execution time.
    elapsed_time = clock.time() - program_start_time
    print('Program Execution Time [s]: {:0.3f}'.format(elapsed_time))

