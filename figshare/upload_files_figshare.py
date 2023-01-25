#! /usr/bin/env python

import argparse
import sys
import json
import figshare

def main(args):
    fconn = figshare.FigshareConnection(chunk_size = args.chunk_size)

    # load json
    if args.article_data is not None:
        f = open(args.article_data, 'r')
        article_data = json.load(f)
        # print(json.dumps(article_data, sort_keys=True, indent=4))
    else:
        article_data = {
            "title": args.title,
        }

    # if no article_id is specified, create a new entry
    if args.article_id is None:
        article_id = fconn.create_article(article_data)
    else:
        article_id = args.article_id

    # Then we upload the files.
    for filename in args.file_path:
        file_info = fconn.initiate_new_upload(article_id, filename)
        # Until here we used the figshare API; following lines use the figshare upload service API.
        fconn.upload_parts(file_info, filename)
        # We return to the figshare API to complete the file upload process.
        fconn.complete_upload(article_id, file_info['id'])
        print('Uploaded file {filename}'.format(filename = filename))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Script to upload files to Figshare using the API',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('file_path', nargs='+', metavar='Files',
        type=str, default=None, help='File(s) to upload [default: None]')
    parser.add_argument('--title', metavar='Str',
        type=str, default='"New Figshare entry"', 
        help='Title for entry')
    parser.add_argument('--article_id', metavar='Int', type=int, 
        default=None, help='Article id to add files to')
    parser.add_argument('--article_data', metavar='JSON FILE', type=str, 
        default=None, help='Name of JSON file of atricle data to use')
    parser.add_argument('--chunk_size', metavar='Int', type=int, 
        default=1048576, help='Size to chunk files by')
    # parser.add_argument('--log_level', action='store', default='INFO',
    #     choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
    #     help='Sets the logging level')
    args = parser.parse_args()
    main(args)