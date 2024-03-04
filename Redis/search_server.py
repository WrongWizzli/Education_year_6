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

from logger import logger


redis_meta = json.load(open('redis_meta.json'))
db_1 = redis.Redis(host='localhost', port=redis_meta[0]['port'], db=0)
db_2 = redis.Redis(host='localhost', port=redis_meta[1]['port'], db=0)


def consume_with_rabbitmq(queue: str, callback: object):
    connection = pika.BlockingConnection(pika.ConnectionParameters('localhost'))
    channel = connection.channel()
    channel.queue_declare(queue=queue)
    channel.basic_consume(queue=queue, on_message_callback=callback)
    logger.info(f'[SEARCH_SERVER] Ready to consume from queue {queue}')
    try:
        channel.start_consuming()
    except pika.exceptions.ConsumerCancelled:
        pass
    

pid_callback_result = None
def recieve_pids_callback(channel, method, properties, body):
    logger.info(f'[SEARCH_SERVER] Recieved body: {body.decode()}')
    if method is not None:
        global pid_callback_result
        pid_callback_result = json.loads(body.decode())
        channel.basic_ack(delivery_tag=method.delivery_tag)
        channel.basic_cancel(consumer_tag=method.consumer_tag)


def recieve_redis_pids() -> list:
    consume_with_rabbitmq('pid_info', recieve_pids_callback)
    global pid_callback_result
    return pid_callback_result


def get_redis_pattern(key: dict):
    redis_pattern = key.get('brand', '*')
    redis_pattern += '_'
    redis_pattern += key.get('colour', '*')
    redis_pattern += '_*'
    return redis_pattern


def copy_files_to_response(ipaths: list):
    folder = 'response'
    for ipath in ipaths:
        iname = ipath.split('/')[-1]
        shutil.copy(ipath, os.path.join(folder, iname))


def process_client_request_callback(channel, method, properties, body):
    logger.info(f'[SEARCH_SERVER] Recieved body: {body.decode()}')
    key_map = json.loads(body.decode())
    logger.info(f'[SEARCH_SERVER] Key_map for request: {key_map}')
    if 'exit' in key_map:
        logger.info(f'[SEARCH_SERVER] Stopping the server...')
        channel.basic_ack(delivery_tag=method.delivery_tag)
        channel.basic_cancel(consumer_tag=method.consumer_tag)
        return
    ipaths_to_response = []
    pattern = get_redis_pattern(key_map)
    keys = db_1.keys(pattern)
    if len(keys):
        ipaths_to_response.extend([elem.decode() for elem in db_1.mget(keys)])
    keys = db_2.keys(pattern)
    if len(keys):
        ipaths_to_response.extend([elem.decode() for elem in db_2.mget(keys)])
    logger.info(f'Images to load: {ipaths_to_response}')
    copy_files_to_response(ipaths_to_response)
    channel.basic_ack(delivery_tag = method.delivery_tag)
    

def run_server():
    redis_pids = recieve_redis_pids()
    logger.info(f'[SEARCH_SERVER] Redis pids: {redis_pids}')
    
    consume_with_rabbitmq('server_client_interaction', process_client_request_callback)

    for redis_pid in redis_pids:
        p = psutil.Process(int(redis_pid))
        p.terminate()

run_server()