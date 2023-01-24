#! /usr/bin/env python

''' Script for testing the Figshare API '''

import figshare
import json
import sys
def main():
    fconn = figshare.FigshareConnection(chunk_size = 10485760)
    print(fconn.chunk_size, fconn.token)

    article = fconn.fetch_article_by_id('11830380')
    print(json.dumps(article, sort_keys=True, indent=4))
    sys.exit(0)
    fconn.list_articles()

    fconn.search_articles(':title: "Entries from OMIM with Clinical Synopsis"')

    article_data = {
        "title": "Test New Entry",
        "keywords": [
            "tag1",
            "tag2"
        ]
    }
    article_id = fconn.create_article(article_data)
    # fconn.list_articles() #Prints the list again and you should see your new item at the top of the list
    # fconn.list_files_of_article(article_id)
    files = ['omim-clin-synopsis.log', 'omim-ph-map-key-3.json', 
        'requirements.txt']
    for filename in files:
        file_info = fconn.initiate_new_upload(article_id, filename)
        print(file_info)
        # Until here we used the figshare API; following lines use the figshare upload service API.
        fconn.upload_parts(file_info, filename)
        # We return to the figshare API to complete the file upload process.
        fconn.complete_upload(article_id, file_info['id'])
        print('Uploaded file {filename}'.format(filename = filename))

    # test adding file to existing article


if __name__ == '__main__':
    main()
