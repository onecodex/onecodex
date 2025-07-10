from typing import ClassVar, Optional, List

from pydantic import BaseModel as PydanticBaseModel
from pydantic import ConfigDict, computed_field, Field


class ApiBaseModel(PydanticBaseModel):
    _api: ClassVar[Optional["Api"]] = (  # noqa: F821
        None  # Forward reference to avoid circular imports
    )
    _client: ClassVar[Optional["HTTPClient"]] = None  # noqa: F821
    _resource_path: ClassVar[str]  # Default resource path, subclasses override

    model_config = ConfigDict(
        # populate_by_name=True,
        populate_by_alias=True,
        extra="allow",
        arbitrary_types_allowed=True,
    )

    field_uri: str = Field(..., alias="$uri", exclude=True)

    @computed_field
    def id(self) -> str:
        return self.field_uri.split("/")[-1]

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.id}>"

    # def _repr_html_(self):
    #         return '''<table>
    #     <thead>
    #         <tr>
    #             <th colspan="2"><code>{cls}({id})</code></th>
    #         </tr>
    #     </thead>
    #     <tbody>{properties}</tbody>
    #     </table>'''.format(
    #         cls=self.__class__.__name__,
    #         id=escape(repr(self.id)),
    #         properties='\n'.join('<tr><td>{}</td><td><code>{}</code></td>'.format(escape(k),
    #                                                                               escape(pformat(v)))
    #                              for k, v in self._properties.items() if not k.startswith('$')))

    @classmethod
    def all(cls) -> List["ApiBaseModel"]:
        raise NotImplementedError("Subclasses must implement this method")

    @classmethod
    def get(cls, id: str) -> "ApiBaseModel":
        resp = cls._client.get(f"{cls._api._base_url}{cls._resource_path}/{id}?expand=all")
        return cls.model_validate(resp.json())

    @classmethod
    def where(cls, **kwargs) -> List["ApiBaseModel"]:
        raise NotImplementedError("Subclasses must implement this method")

    @classmethod
    def create(cls, **kwargs):
        raise NotImplementedError("Subclasses must implement this method")

    def update(self, id: str, **kwargs):
        raise NotImplementedError("Subclasses must implement this method")

    def delete(self, id: str):
        raise NotImplementedError("Subclasses must implement this method")

    def save(self):
        raise NotImplementedError("Subclasses must implement this method")
