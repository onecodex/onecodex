# Run this against localhost logged in the ASM NGS user to get everything
import json
import requests
from requests.auth import HTTPBasicAuth

# asmngs@onecodex demo user account (test key!)
auth = HTTPBasicAuth("1eab4217d30d42849dbde0cd1bb94e39", "")
base_url = "http://localhost:3000"


def fetch(resource):
    return requests.get(base_url + resource, auth=auth).json()


schema = fetch("/api/v1/schema")
with open("schema.json", mode="w") as f:
    f.write(json.dumps(schema))


for resource in schema["properties"].keys():
    print("Fetching resource: {}".format(resource))
    instances = fetch("/api/v1/" + resource + "?per_page=500")
    with open("{}.json".format(resource), mode="w") as f:
        f.write(json.dumps(instances))

    schema = fetch("/api/v1/" + resource + "/schema")
    with open("schema_{}.json".format(resource), mode="w") as f:
        f.write(json.dumps(schema))
