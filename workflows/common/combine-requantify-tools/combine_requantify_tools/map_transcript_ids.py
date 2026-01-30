import logging
import re
from pathlib import Path

from combine_requantify_tools.types import (
    Gtf,
    TrackingToProcess,
    TypedDataFrame,
)
from combine_requantify_tools.utils import (
    read_gtf_from_file,
    read_tracking_file,
    type_data_frame,
)


def do_map_transcript_ids(
    tracking_in: Path, gtf_in: str | Path, sample_ids: list[str]
) -> tuple[TypedDataFrame[TrackingToProcess], TypedDataFrame[Gtf]]:
    """
    Replace TCONS_* with annotated ENST* transcript IDs in tracking and GTF data.

    :param tracking_in: path to the tracking file
    :param gtf_in: path to the GTF file
    :param sample_ids: list of sample IDs (in the order expected in the tracking file)
    :return: updated data frames for tracking and GTF with updated transcript IDs
    """

    # read the input tracking file
    tracking = read_tracking_file(tracking_in, sample_ids)

    logging.info("Making transcript ID map")
    # split the transcript ID column into components
    tracking["transcript_id_parts"] = tracking["transcript_id"].str.split("|")

    # collect rows to be used for replacing IDs
    mapping_rows = tracking.loc[
        tracking["val"].eq("="), ["transcript_id_parts", "gene_id"]
    ].copy()

    # old_id = transcript_id before "|" (e.g. "TCONS_00002258")
    mapping_rows["old_id"] = mapping_rows["transcript_id_parts"].str.get(0)

    # new_id = gene_id split on "|", middle element (e.g. "ENST00000367067.8")
    mapping_rows["new_id"] = mapping_rows["gene_id"].str.split("|").str.get(1)

    id_map = dict(zip(mapping_rows["old_id"], mapping_rows["new_id"]))

    # rewrite transcript_id (e.g. "TCONS_00002258|9|3274"=>"ENST00000367067.8|9|3274")
    logging.info("Replacing transcript IDs in tracking")
    prefix = tracking["transcript_id_parts"].str.get(0)
    suffix = (
        tracking["transcript_id_parts"].str.get(1)
        + "|"
        + tracking["transcript_id_parts"].str.get(2)
    )

    new_prefix = prefix.map(id_map).fillna(prefix)

    tracking["transcript_id"] = new_prefix + suffix.fillna("|" + suffix)

    # finished processing the tracking data frame
    tracking = tracking.drop(columns="transcript_id_parts")
    tracking = type_data_frame(tracking, TrackingToProcess)

    # load the GTF
    gtf = read_gtf_from_file(gtf_in)

    # regex to find the transcript_id component and its quoted value
    tid_re = re.compile(r'transcript_id "([^"]+)"')

    def replace_attribute(attribute: str) -> str:
        """
        Replace transcript_id values inside a single attribute string.

        :param attribute: the "attribute" column in a single row of the GTF
        :return: the attribute string, potentially with the transcript_id replaced
        """

        def repl(match: re.Match) -> str:
            """
            Check for a transcrip_id and replace it in its regex match group

            :param match: a regex match object
            :return: the attribute string, potentially with the transcript_id replaced
            """

            tid = match.group(1)

            if tid in id_map:
                # replace using the mapping
                return f'transcript_id "{id_map[tid]}"'
            else:
                # just return original string
                return match.group(0)

        return tid_re.sub(repl, attribute)

    # replace the transcript_ids
    logging.info("Replacing transcript IDs in GTF")
    gtf["attribute"] = gtf["attribute"].apply(replace_attribute)

    return tracking, gtf
