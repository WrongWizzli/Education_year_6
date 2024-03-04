import logging

logging.basicConfig(
    level=logging.INFO,
    filename='log.txt',
    filemode='w',
    format='[%(levelname)s] %(asctime)s, %(funcName)s: %(message)s'
)

logger = logging.getLogger(__name__)