import base64
import json
import logging

import pandas as pd
from cloudevents.http import CloudEvent

from dogspa_long_reads.app import entrypoint

if __name__ == "__main__":
    # set this Pandas option now as a development best practice
    pd.set_option("display.max_columns", 30)
    pd.set_option("display.max_colwidth", 50)
    pd.set_option("display.max_info_columns", 30)
    pd.set_option("display.max_info_rows", 20)
    pd.set_option("display.max_rows", 20)
    pd.set_option("display.max_seq_items", None)
    pd.set_option("display.width", 200)
    pd.set_option("expand_frame_repr", True)
    pd.set_option("mode.chained_assignment", "warn")

    logger = logging.getLogger()
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.DEBUG)

    payload = {
        "workspace": {
            "namespace": "broad-firecloud-ccle",
            "name": "Long_Read_Omics",
        },
        "gcs_source": {
            "bucket": "fc-aaf4de93-c104-45c4-a01a-a036869119c6",
            "prefix": "tag2028/merge/",
        },
        "gcs_destination": {
            "bucket": "cclebams",
            "prefix": "rna_long_read/",
        },
        "gcp_project": "depmap-omics",
        "uuid_namespace": "00000000-0000-0000-0000-000000000000",
        "dry_run": True,
    }

    # encode to match pub/sub functionality
    ce_data = base64.b64encode(json.dumps(payload).encode())

    cloud_event = CloudEvent(
        attributes={
            "specversion": "1.0",
            "type": "com.github.pull_request.opened",
            "source": "https://github.com/cloudevents/spec/pull",
            "id": "A234-1234-1234",
            "time": "2018-04-05T17:31:00Z",
        },
        data={"message": {"data": ce_data}},
    )

    entrypoint(cloud_event)
