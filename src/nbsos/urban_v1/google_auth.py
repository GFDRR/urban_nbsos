from . import loggerconfig

import json
import os
import ee
from google.cloud import storage
import logging

logger = logging.getLogger(__name__)
logger = loggerconfig.config_logger(logger)

class Authenticate:

    def __init__(self, credentials_json):
        """
		Authenticate with Google cloud to access data stored on cloud storage
		:param credentials_json:
			local json file with authentication details
		"""
        self.credentials_json = credentials_json
        logger.info('Opening credentials json file')
        try:
            with open(credentials_json) as cred:
                cred = json.load(cred)
        except Exception as e:
            logger.error(f'Could not open the provided credentials json file: {e}')
            raise e
        client_email = cred['client_email']
        project = cred['project_id']
        os.environ['SERVICE_ACCOUNT'] = client_email
        os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = credentials_json
        os.environ['CLOUD_PROJECT'] = project
        self.credentials = None
        self.client = None
        self.bucket = None

    def start_google_cloud_session(self, bucket):
        """
		Start a GCS session with specified bucket
		:param bucket: name of the bucket you want to access
		:return: GCS bucket object
		"""

        # get specified bucket from gcs
        logger.info(f'Connect to Google Cloud Bucket: {bucket}')
        try:
            self.client = storage.Client()
            self.bucket = self.client.get_bucket(bucket)
        except Exception as e:
            logging.error(f"Could not connect to {bucket}: {e}")
            raise e

        return self.bucket

    def initialize_ee(self):
        """
        Start an ee session with the account specified in the credentials.json
        :return: -
        """
        # set credetials for google earth engine
        logger.info(f"Sarting Google Earth Engine Session")
        try:
            self.credentials = ee.ServiceAccountCredentials(
                os.environ['SERVICE_ACCOUNT'],
                os.environ['GOOGLE_APPLICATION_CREDENTIALS']
            )

            # initialize google earth engine
            ee.Initialize(self.credentials)
        except Exception as e:
            logger.warning(f'Could not start Google Earth Engine Session: {e}')
            pass

class AuthenticateIndividual:
    def __init__(self, project):
        self.project = project

    def initialize_ee_indiviudal(self, project_name):
        logger.info(f"Sarting Google Earth Engine Session")
        try:
            ee.Authenticate()
            ee.Initialize(project=project_name
                          #opt_url='https://earthengine-highvolume.googleapis.com'
             )
        except Exception as e:
            logger.warning(f'Could not start Google Earth Engine Session: {e}')
            raise e









