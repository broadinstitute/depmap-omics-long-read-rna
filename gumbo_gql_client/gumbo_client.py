from typing import Any, Dict, List

from .base_client import BaseClient
from .get_models_and_children import GetModelsAndChildren
from .input_types import omics_sequencing_insert_input
from .insert_omics_sequencings import InsertOmicsSequencings


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
                model_conditions {
                  model_condition_id
                  omics_profiles {
                    profile_id
                    main_sequencing_id
                    datatype
                    blacklist_omics
                    omics_sequencings {
                      sequencing_id
                      legacy_size
                      sequencing_id
                      expected_type
                      bam_filepath
                      bai_filepath
                      blacklist
                      source
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

    def insert_omics_sequencings(
        self, username: str, objects: List[omics_sequencing_insert_input], **kwargs: Any
    ) -> InsertOmicsSequencings:
        query = gql(
            """
            mutation InsertOmicsSequencings($_username: String!, $objects: [omics_sequencing_insert_input!]!) {
              set_username(args: {_username: $_username}) {
                username
              }
              insert_omics_sequencing(objects: $objects) {
                affected_rows
              }
            }
            """
        )
        variables: Dict[str, object] = {"_username": username, "objects": objects}
        response = self.execute(
            query=query,
            operation_name="InsertOmicsSequencings",
            variables=variables,
            **kwargs
        )
        data = self.get_data(response)
        return InsertOmicsSequencings.model_validate(data)
