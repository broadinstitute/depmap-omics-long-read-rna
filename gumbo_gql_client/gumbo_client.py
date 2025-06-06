from typing import Any, Dict, List

from .base_client import BaseClient
from .input_types import sequencing_alignment_insert_input
from .insert_sequencing_alignments import InsertSequencingAlignments
from .mapped_short_long_rna_samples import MappedShortLongRnaSamples


def gql(q: str) -> str:
    return q


class GumboClient(BaseClient):
    def mapped_short_long_rna_samples(self, **kwargs: Any) -> MappedShortLongRnaSamples:
        query = gql(
            """
            query MappedShortLongRnaSamples {
              records: omics_mapping(where: {datatype: {_in: ["rna", "long_read_rna"]}}) {
                model_id
                model_condition_id
                omics_profile_id
                omics_sequencing_id
                datatype
                priority
                model {
                  cell_line_name
                  stripped_cell_line_name
                }
                omics_sequencing {
                  sequencing_alignments {
                    sequencing_alignment_id: id
                    crai_bai_url: index_url
                    cram_bam_url: url
                    reference_genome
                    sequencing_alignment_source
                    size
                  }
                }
              }
            }
            """
        )
        variables: Dict[str, object] = {}
        response = self.execute(
            query=query,
            operation_name="MappedShortLongRnaSamples",
            variables=variables,
            **kwargs
        )
        data = self.get_data(response)
        return MappedShortLongRnaSamples.model_validate(data)

    def insert_sequencing_alignments(
        self,
        username: str,
        objects: List[sequencing_alignment_insert_input],
        **kwargs: Any
    ) -> InsertSequencingAlignments:
        query = gql(
            """
            mutation InsertSequencingAlignments($_username: String!, $objects: [sequencing_alignment_insert_input!]!) {
              set_username(args: {_username: $_username}) {
                username
              }
              insert_sequencing_alignment(objects: $objects) {
                affected_rows
              }
            }
            """
        )
        variables: Dict[str, object] = {"_username": username, "objects": objects}
        response = self.execute(
            query=query,
            operation_name="InsertSequencingAlignments",
            variables=variables,
            **kwargs
        )
        data = self.get_data(response)
        return InsertSequencingAlignments.model_validate(data)
