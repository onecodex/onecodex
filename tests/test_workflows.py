from onecodex.models import Workflows


def test_query_for_workflows(ocx, api_data):
    sample_id = "002f55a11aae4b8f"
    workflow = ocx.Workflows.where(sample=sample_id)
    assert isinstance(workflow[0], Workflows)

    workflows = ocx.Workflows.all()
    assert len(workflows) == 1


def test_workflows_results(ocx, api_data):
    workflow_id = "abcde12345fedcba"
    workflow = Workflows.get(workflow_id)

    assert isinstance(workflow, Workflows)

    result = workflow.results()

    assert "workflowResults" in result
    assert result["workflowResults"]["summary"]["status"] == "complete"
    assert result["workflowResults"]["steps"][0]["name"] == "step_one"
