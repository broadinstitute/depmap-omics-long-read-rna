from typing import Any, Dict, List

from .base_client import BaseClient
from .get_models_and_children import GetModelsAndChildren
from .input_types import sequencing_alignment_insert_input
from .insert_sequencing_alignments import InsertSequencingAlignments
from .long_read_sequencing_alignments import LongReadSequencingAlignments


def gql(q: str) -> str:
    return q


class GumboClient(BaseClient):
    def get_models_and_children(self, **kwargs: Any) -> GetModelsAndChildren:
        query = gql(
            """
            query GetModelsAndChildren {
              records: model {
                model_id
                cell_line_name
                stripped_cell_line_name
                model_conditions(
                  where: {_and: [{_or: [{drug: {_is_null: true}}, {drug: {_eq: "DMSO"}}]}, {_or: [{expansion_team: {_is_null: true}}, {expansion_team: {_neq: "PRISM"}}]}]}
                ) {
                  model_condition_id
                  drug
                  expansion_team
                  source
                  omics_profiles(
                    where: {datatype: {_in: ["rna", "long_read_rna"]}, blacklist_omics: {_neq: true}}
                  ) {
                    omics_profile_id: profile_id
                    datatype
                    omics_sequencings(where: {blacklist: {_neq: true}}) {
                      omics_sequencing_id: sequencing_id
                      stranded
                      version
                      sequencing_alignments {
                        sequencing_alignment_id: id
                        crai_bai_url: index_url
                        cram_bam_url: url
                        reference_genome
                        sequencing_alignment_source
                      }
                    }
                  }
                }
              }
            }
            """
        )
        variables: Dict[str, object] = {}
        response = self.execute(
            query=query,
            operation_name="GetModelsAndChildren",
            variables=variables,
            **kwargs
        )
        data = self.get_data(response)
        return GetModelsAndChildren.model_validate(data)

    def long_read_sequencing_alignments(
        self, **kwargs: Any
    ) -> LongReadSequencingAlignments:
        query = gql(
            """
            query LongReadSequencingAlignments {
              records: omics_sequencing(where: {expected_type: {_eq: "long_read_rna"}}) {
                omics_sequencing_id: sequencing_id
                sequencing_alignments {
                  sequencing_alignment_source
                  size
                }
              }
            }
            """
        )
        variables: Dict[str, object] = {}
        response = self.execute(
            query=query,
            operation_name="LongReadSequencingAlignments",
            variables=variables,
            **kwargs
        )
        data = self.get_data(response)
        return LongReadSequencingAlignments.model_validate(data)

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
