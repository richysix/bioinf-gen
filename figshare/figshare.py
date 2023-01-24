#! /usr/bin/env python

''' Module for interacting with the Figshare API 

This script is modified from

https://colab.research.google.com/drive/13CAM8mL1u7ZsqNhfZLv7bNb1rdhMI64d

See https://docs.figshare.com/#upload_files
'''

import hashlib
import json
import os
import sys
import requests
from requests.exceptions import HTTPError

class FigshareConnection:
    base_url = 'https://api.figshare.com/v2/{endpoint}'
    def __init__(self, chunk_size):
        self.chunk_size = chunk_size
        self.token = self.get_api_token()

    def get_api_token(self):
        API_TOKEN = os.getenv('FIGSHARE_TOKEN')
        if API_TOKEN is None:
            raise Exception("The required environment variable FIGSHARE_TOKEN is not set")
        else:
            return(API_TOKEN)

    def raw_issue_request(self, method, url, data=None, binary=False):
        headers = {'Authorization': 'token ' + self.token}
        if data is not None and not binary:
            data = json.dumps(data)
        response = requests.request(method, url, headers=headers, data=data)
        try:
            response.raise_for_status()
            try:
                data = json.loads(response.content)
            except ValueError:
                data = response.content
        except HTTPError as error:
            print('Caught an HTTPError: {}'.format(error.message))
            print('Body:\n', response.content)
            raise

        return data

    def issue_request(self, method, endpoint, *args, **kwargs):
        url = self.base_url.format(endpoint=endpoint)
        return self.raw_issue_request(method, url, *args, **kwargs)

    def list_articles(self):
        result = self.issue_request('GET', 'account/articles')
        print('Listing current articles:')
        if result:
            for item in result:
                print(u'  {url} - {title}'.format(**item))
        else:
            print('  No articles.')

    def list_files_of_article(self, article_id):
        result = self.issue_request('GET', 'account/articles/{}/files'.format(article_id))
        print('Listing files for article {}:'.format(article_id))
        if result:
            for item in result:
                print('  {id} - {name}'.format(**item))
        else:
            print('  No files.')

    def create_article(self, article_data):
        result = self.issue_request('POST', 'account/articles', data=article_data)
        print('Created article:', result['location'], '\n')

        result = self.raw_issue_request('GET', result['location'])
        # print(result)

        return result['id']

    def get_file_check_data(self, file_path):
        with open(file_path, 'rb') as fin:
            md5 = hashlib.md5()
            size = 0
            data = fin.read(self.chunk_size)
            while data:
                size += len(data)
                md5.update(data)
                data = fin.read(self.chunk_size)
            return md5.hexdigest(), size

    def initiate_new_upload(self, article_id, file_path):
        endpoint = 'account/articles/{}/files'
        endpoint = endpoint.format(article_id)

        md5, size = self.get_file_check_data(file_path)
        data = {'name': os.path.basename(file_path),
                'md5': md5,
                'size': size}

        result = self.issue_request('POST', endpoint, data=data)
        print('Initiated file upload:', result['location'], '\n')

        result = self.raw_issue_request('GET', result['location'])

        return result

    def upload_parts(self, file_info, file_path):
        url = '{upload_url}'.format(**file_info)
        result = self.raw_issue_request('GET', url)

        print('Uploading parts:')
        with open(file_path, 'rb') as fin:
            for part in result['parts']:
                self.upload_part(file_info, fin, part)

    def upload_part(self, file_info, stream, part):
        udata = file_info.copy()
        udata.update(part)
        url = '{upload_url}/{partNo}'.format(**udata)

        stream.seek(part['startOffset'])
        data = stream.read(part['endOffset'] - part['startOffset'] + 1)

        self.raw_issue_request('PUT', url, data=data, binary=True)
        print('  Uploaded part {partNo} from {startOffset} to {endOffset}'.format(**part))

    def complete_upload(self, article_id, file_id):
        self.issue_request('POST', 'account/articles/{}/files/{}'.format(article_id, file_id))

    def search_articles(self, search_string):
        kwargs = {
            'data': {
                'search_for': search_string
            }
        }
        result = self.issue_request('POST', 'account/articles/search', **kwargs)
        return result

    def fetch_article_by_id(self, article_id):
        endpoint = 'account/articles/{}'
        endpoint = endpoint.format(article_id)
        result = self.issue_request('GET', endpoint)
        return result
