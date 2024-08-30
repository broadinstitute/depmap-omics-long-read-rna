import logging

import functions_framework
import google.cloud.logging
from cloudevents.http import CloudEvent

from dogspa_long_reads.app import entrypoint


@functions_framework.cloud_event
def run(cloud_event: CloudEvent) -> None:
    """
    Wrapper around the primary `entrypoint` function (needed for remote execution
    inside a GCP Function).

    :param cloud_event: the pub/sub CloudEvent payload
    """

    client = google.cloud.logging.Client()
    client.setup_logging(log_level=logging.INFO)

    entrypoint(cloud_event)
