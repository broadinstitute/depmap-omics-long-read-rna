from typing import Any, Dict, List

from .base_client import BaseClient
from .get_models_and_children import GetModelsAndChildren
from .get_onboarding_workspace import GetOnboardingWorkspace
from .input_types import omics_sequencing_insert_input, onboarding_job_insert_input
from .insert_omics_sequencings import InsertOmicsSequencings
from .insert_onboarding_job import InsertOnboardingJob
from .set_status import SetStatus


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
                    datatype
                    blacklist_omics
                    omics_order_date
                    smid_ordered
                    smid_returned
                    omics_sequencings {
                      blacklist
                      expected_type
                      sequencing_id
                      source
                      version
                      sequencing_alignments {
                        id
                        url
                        index_url
                        size
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

    def set_status(
        self, username: str, profile_ids: List[str], status: str, **kwargs: Any
    ) -> SetStatus:
        query = gql(
            """
            mutation SetStatus($_username: String!, $profile_ids: [String!]!, $status: String!) {
              set_username(args: {_username: $_username}) {
                username
              }
              update_omics_profile_many(
                updates: {where: {profile_id: {_in: $profile_ids}}, _set: {status: $status}}
              ) {
                affected_rows
              }
            }
            """
        )
        variables: Dict[str, object] = {
            "_username": username,
            "profile_ids": profile_ids,
            "status": status,
        }
        response = self.execute(
            query=query, operation_name="SetStatus", variables=variables, **kwargs
        )
        data = self.get_data(response)
        return SetStatus.model_validate(data)

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

    def insert_onboarding_job(
        self, username: str, object: onboarding_job_insert_input, **kwargs: Any
    ) -> InsertOnboardingJob:
        query = gql(
            """
            mutation InsertOnboardingJob($_username: String!, $object: onboarding_job_insert_input!) {
              set_username(args: {_username: $_username}) {
                username
              }
              insert_onboarding_job_one(object: $object) {
                id
              }
            }
            """
        )
        variables: Dict[str, object] = {"_username": username, "object": object}
        response = self.execute(
            query=query,
            operation_name="InsertOnboardingJob",
            variables=variables,
            **kwargs
        )
        data = self.get_data(response)
        return InsertOnboardingJob.model_validate(data)

    def get_onboarding_workspace(self, **kwargs: Any) -> GetOnboardingWorkspace:
        query = gql(
            """
            query GetOnboardingWorkspace {
              records: onboarding_workspace(where: {expected_type: {_eq: "long_read_rna"}}) {
                id
                excluded_terra_sample_ids
                min_file_size
              }
            }
            """
        )
        variables: Dict[str, object] = {}
        response = self.execute(
            query=query,
            operation_name="GetOnboardingWorkspace",
            variables=variables,
            **kwargs
        )
        data = self.get_data(response)
        return GetOnboardingWorkspace.model_validate(data)
