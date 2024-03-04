import datetime
from glob import glob
import json
import os
import pika
import psutil
import random
import redis
import shutil
import time


def main():
    print('============How to use============')
    print('You can request all images by simply clickng Enter.')
    print('If you want to specify the set of images, you need')
    print(' to set the parameters in command line in json format.')
    print('Json may have two keys: brand and colour, example: {"colour": "white"}')
    print('To stop the interaction type "exit".')
    print('===========/How to use============')
    print('Type your request:')
    connection = pika.BlockingConnection(pika.ConnectionParameters('localhost'))
    channel = connection.channel()
    channel.queue_declare(queue='server_client_interaction')
    
    while True:
        request = input()
        if request == 'exit':
            channel.basic_publish('', 'server_client_interaction', json.dumps({'exit': True}))
            break
        try:
            json.loads(request)
            channel.basic_publish('', 'server_client_interaction', request)
        except json.JSONDecodeError:
            print('Your input is incorrect! Watch \'How to Use\' and try again...')
    
    connection.close()


if __name__ == '__main__':
    main()