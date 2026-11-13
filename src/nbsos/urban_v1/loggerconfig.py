import logging
import shutil
import os



def config_logger(logger, mode='a'):
    logger.setLevel(logging.INFO)

    # define file handler and set formatter
    file_handler = logging.FileHandler('log_file.log', mode=mode)
    # formatter = logging.Formatter('[%(asctime)s] | %(levelname)-8s | %(name)-23s | %(message)s', "%Y-%m-%d %H:%M:%S")
    formatter = logging.Formatter('[%(asctime)s] %(levelname)-8s | %(filename)-25s:%(funcName)-40s() > %(message)s', "%Y-%m-%d %H:%M:%S")
    
    file_handler.setFormatter(formatter)
    logger.propagate = False

    # add file handler to logger
    logger.addHandler(file_handler)

    return logger

def save_local_logfile(local_root, country, city) -> str:
    src = 'log_file.log'
    logs_folder = os.path.join(local_root, country, city, 'logs')
    previous_logs = len(os.listdir(logs_folder))
    version = previous_logs + 1
    logfile_name = os.path.join(logs_folder, f'nbsos_log_V{version}.log')
    shutil.copyfile(src, logfile_name)
    return logfile_name

def save_logfile_to_cloud(logfile_name, country, city, bucket):
    try:
        cloud_name = f'nbsos/{country}/{city}/logs/{os.path.basename(logfile_name)}'
        blob = bucket.blob(cloud_name)
        blob.upload_from_filename(logfile_name)
    except Exception as e:
        logging.error(f'Exporting {logfile_name} to GCS: {cloud_name} failed: {e}', exc_info=True)

