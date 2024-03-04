import datetime
from glob import glob
import json
import os
import pika
import random
import redis
import time

from PIL import Image

from logger import logger


RAW_DATA_FOLDER = './raw_data'


class TooLongLoadError(Exception):
    pass


def extract():
    result = {}
    files = glob(os.path.join(RAW_DATA_FOLDER, '*'))
    for fname in files:
        result[fname] = Image.open(fname)
    return result


def transform(extracted_data):
    for fname in extracted_data:
        extracted_data[fname] = extracted_data[fname].convert('RGB')
    return extracted_data


def get_new_name(new_folder_name: str, old_fname: str) -> str:
    file_name = old_fname.split('/')[-1]
    return os.path.join(new_folder_name, file_name)


def try_to_load_to_redis(db: redis.Redis, new_fname: str, load_start: datetime.datetime):
    key = new_fname.split('/')[-1].split('.')[0]
    while True:
        try:
            db.set(key, new_fname)
            break
        except redis.exceptions.ConnectionError:
            logger.error('[ETL] Failed to make set() request redis db...')
            now = datetime.datetime.now()
            if (now - load_start).seconds > 600:
                raise TooLongLoadError('Redis server does not respond for 600 seconds. Cannot load data')
            time.sleep(2)


def try_to_flush_redis(db: redis.Redis, load_start: datetime.datetime):
    while True:
        try:
            db.flushdb(True)
            break
        except redis.exceptions.ConnectionError:
            logger.error('[ETL] Failed to make flushdb() request redis db...')
            now = datetime.datetime.now()
            if (now - load_start).seconds > 10:
                raise TooLongLoadError('Redis server does not respond for 600 seconds. Cannot load data')
            time.sleep(2)


def load(extracted_data):
    load_start = datetime.datetime.now()
    
    redis_meta = json.load(open('redis_meta.json'))
    
    redis_1_meta = redis_meta[0]
    redis_folder_1 = redis_1_meta['file_folder']
    db_1 = redis.Redis(host='localhost', port=redis_1_meta['port'], db=0)
    try_to_flush_redis(db_1, load_start)
    redis_2_meta = redis_meta[1]
    redis_folder_2 = redis_2_meta['file_folder']
    db_2 = redis.Redis(host='localhost', port=redis_2_meta['port'], db=0)
    try_to_flush_redis(db_2, load_start)
    
    fnames = list(extracted_data.keys())
    for fname in fnames[:len(fnames) // 2]:
        new_fname = get_new_name(redis_folder_1, fname)
        extracted_data[fname].save(new_fname)
        try_to_load_to_redis(db_1, new_fname, load_start)
    for fname in fnames[len(fnames) // 2:]:
        new_fname = get_new_name(redis_folder_2, fname)
        extracted_data[fname].save(new_fname)
        try_to_load_to_redis(db_2, new_fname, load_start)


def send_time_to_parent(start: float):
    connection = pika.BlockingConnection(pika.ConnectionParameters('localhost'))
    channel = connection.channel()
    channel.queue_declare(queue='start_time_report')
    channel.basic_publish(exchange='', routing_key='start_time_report', body=f'{time.time() - start}')
    connection.close()


start = time.time()

data = extract()
logger.info('[ETL] Data has been extracted')
data = transform(data)
logger.info('[ETL] Data has been transformed')
load(data)
logger.info('[ETL] Data has been loaded')

send_time_to_parent(start)
                
