from glob import glob
import json
import os
import pika
import time
import subprocess

from logger import logger


REDIS_FOLDER_1 = './redis_1'
REDIS_FOLDER_2 = './redis_2'
RESPONSE = './response'


def remove_files_from_last_launch(folder_name: str):
    old_files = glob(f'./{folder_name}/*')
    for fname in old_files:
        os.remove(fname)
    os.rmdir(folder_name)

def consume_with_rabbitmq(queue: str, callback: object):
    connection = pika.BlockingConnection(pika.ConnectionParameters('localhost'))
    channel = connection.channel()
    channel.queue_declare(queue=queue)
    channel.basic_consume(queue=queue, on_message_callback=callback)
    print(f'[SEARCH_SERVER] Ready to consume from queue {queue}')
    try:
        channel.start_consuming()
    except pika.exceptions.ConsumerCancelled:
        pass


time_callback_result = None
def recieve_etl_exec_time_callback(channel, method, properties, body):
    print(f'[MAIN] Recieved body: {body.decode()}')
    if method is None:
        return
    global time_callback_result
    time_callback_result = float(body.decode())
    channel.basic_ack(delivery_tag=method.delivery_tag)
    channel.basic_cancel(consumer_tag=method.consumer_tag)


def receive_execution_time_from_etl():
    consume_with_rabbitmq('start_time_report', recieve_etl_exec_time_callback)
    global time_callback_result
    return time_callback_result


def send_redis_pids_to_search_server(pids: list):
    connection = pika.BlockingConnection(pika.ConnectionParameters('localhost'))
    channel = connection.channel()
    channel.queue_declare(queue='pid_info')
    channel.basic_publish(exchange='', routing_key='pid_info', body=json.dumps(pids))
    connection.close()


def main():
    if os.path.exists(REDIS_FOLDER_1):
        remove_files_from_last_launch(REDIS_FOLDER_1)
    os.mkdir(REDIS_FOLDER_1)
    if os.path.exists(REDIS_FOLDER_2):
        remove_files_from_last_launch(REDIS_FOLDER_2)
    os.mkdir(REDIS_FOLDER_2)
    if os.path.exists(RESPONSE):
        remove_files_from_last_launch(RESPONSE)
    os.mkdir(RESPONSE)
    print('[MAIN] Cleaned redis image folders')
    
    redis_settings = json.load(open('redis_meta.json'))
    port_proc_1 = subprocess.Popen(['redis-server', '--port', str(redis_settings[0]['port'])])
    port_proc_2 = subprocess.Popen(['redis-server', '--port', str(redis_settings[1]['port'])])
    print(f'[MAIN] Ran redis servers. Pids: {port_proc_1.pid, port_proc_2.pid}')
    
    etl_proc = subprocess.Popen(['python3', 'etl.py'])
    etl_proc.wait()
    etl_exec_time = receive_execution_time_from_etl()
    print(f'[MAIN] ETL has finished. Elapsed time: {round(etl_exec_time, 5)} (s).')
    
    redis_pids = [port_proc_1.pid, port_proc_2.pid]
    print(redis_pids)

    ss = subprocess.Popen(['python3', 'search_server.py'])
    print(f'[MAIN] Search server pid: {ss.pid}')
    send_redis_pids_to_search_server(redis_pids)


if __name__ == '__main__':
    main()