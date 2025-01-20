import os
import yaml
import asyncio

from dotenv import load_dotenv
from notion_client import AsyncClient
from notion_client.helpers import async_collect_paginated_api
from dataclasses import dataclass, fields
from typing import Any, Type

YES = "Yes"
NO = "No"

load_dotenv()


notion = AsyncClient(auth=os.environ["NOTION_TOKEN"])
database_id = os.environ["NOTION_SERIES_DATABASE_ID"]


@dataclass
class SelectPropertyValue:
    id: str
    name: str
    color: str


@dataclass
class NotionDatabaseProperty:
    id: str
    type: str


@dataclass
class SelectProperty(NotionDatabaseProperty):
    select: SelectPropertyValue | None

    @staticmethod
    def from_dict(property: dict) -> "SelectProperty":
        return SelectProperty(
            id=property["id"],
            type=property["type"],
            select=(
                SelectPropertyValue(**property["select"])
                if property.get("select")
                else None
            ),
        )


@dataclass
class MultiSelectProperty(NotionDatabaseProperty):
    multi_select: list[SelectPropertyValue]

    @staticmethod
    def from_dict(property: dict) -> "SelectProperty":
        return MultiSelectProperty(
            id=property["id"],
            type=property["type"],
            multi_select=[
                SelectPropertyValue(**select)
                for select in property.get("multi_select", [])
            ],
        )



@dataclass
class LinkProperty(NotionDatabaseProperty):
    url: str

    @staticmethod
    def from_dict(property: dict) -> "LinkProperty":
        return LinkProperty(
            id=property["id"], type=property["type"], url=property["url"]
        )


@dataclass
class RichTextProperty(NotionDatabaseProperty):
    rich_text: list[dict]

    @staticmethod
    def from_dict(property: dict) -> "RichTextProperty":
        return RichTextProperty(
            id=property["id"],
            type=property["type"],
            rich_text=property["rich_text"],
        )


@dataclass
class NumberProperty(NotionDatabaseProperty):
    number: int | None

    @staticmethod
    def from_dict(property: dict) -> "NumberProperty":
        return NumberProperty(
            id=property["id"],
            type=property["type"],
            number=int(property["number"]) if property.get("number") is not None else None,
        )


@dataclass
class TitleProperty(NotionDatabaseProperty):
    title: list["NotionText"]

    @staticmethod
    def from_dict(property: dict) -> "TitleProperty":

        return TitleProperty(
            id=property["id"],
            type=property["type"],
            title=[NotionText.from_dict(text) for text in property["title"]],
        )


@dataclass
class NotionText:
    type: str
    text: dict
    plain_text: str
    href: str

    def from_dict(text: dict) -> "NotionText":
        return NotionText(
            type=text["type"],
            text=text["text"],
            plain_text=text["plain_text"],
            href=text["href"],
        )


@dataclass
class SeriesPageProperties:
    Tags: MultiSelectProperty
    AgePresent: SelectProperty
    SexPresent: SelectProperty
    Link: LinkProperty
    Platform: MultiSelectProperty
    Title: RichTextProperty
    Samples: NumberProperty
    Name: TitleProperty


def dict_to_series_page_property(property: dict) -> NotionDatabaseProperty:
    if property["type"] == "select":
        return SelectProperty.from_dict(property)
    if property["type"] == "multi_select":
        return MultiSelectProperty.from_dict(property)
    if property["type"] == "url":
        return LinkProperty.from_dict(property)
    if property["type"] == "rich_text":
        return RichTextProperty.from_dict(property)
    if property["type"] == "number":
        return NumberProperty.from_dict(property)
    if property["type"] == "title":
        return TitleProperty.from_dict(property)

    raise ValueError(f"Unsupported property type: {property['type']}")


def create_series_page_properties(properties: dict) -> SeriesPageProperties:
    keys = fields(SeriesPageProperties)
    data = {
        key.name: dict_to_series_page_property(properties[key.name])
        for key in keys
    }
    return SeriesPageProperties(**data)


@dataclass
class SeriesPage:
    page_id: str
    created_time: str
    last_edited_time: str
    parent: dict
    archived: bool
    in_trash: bool
    url: str
    properties: SeriesPageProperties


@dataclass
class SeriesItem:
    id: str
    title: str
    samples: int
    parser: dict
    age_present: bool = False
    sex_present: bool = False
    query: str = None
    platforms: list[dict] = None
    tags: list[str] = None


async def get_series_dataset_from_notion() -> list[SeriesPage]:
    pages = await async_collect_paginated_api(
        notion.databases.query, database_id=database_id
    )
    result = [map_repsonse_series_page(page) for page in pages]
    return [page for page in result if len(page.properties.Name.title) > 0]


def get_data_library_file(name: str) -> str:
    file = f"../biolearn/data/{name}"
    base_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(base_dir, file))


def merge_libraries(target_library: list[SeriesItem], source_library: list[SeriesItem]) -> list[SeriesItem]:
    """Merge two lists of SeriesItem with source taking precedence on conflicts.
    
    Args:
        target_library: The base list of SeriesItem
        source_library: The list of SeriesItem whose items take precedence
        
    Returns:
        list[SeriesItem]: Merged list of SeriesItem
    """
    merged_items = {item.id: item for item in target_library}
    
    for item in source_library:
        merged_items[item.id] = item
        
    return list(merged_items.values())


def read_series_dataset_from_library() -> list[SeriesItem]:
    autoscan_library = get_data_library_file("geo_autoscan_library.yaml")
    library = get_data_library_file("library.yaml")
    
    with open(autoscan_library, "r") as file:
        autoscan_library_data = yaml.safe_load(file)
        autoscan_library_items = [dict_to_dataclass(SeriesItem, item) for item in autoscan_library_data["items"]]

        # Check if age and sex metadata keys are present in the parser and set corresponding flags
        for item in autoscan_library_items:
            if "age" in item.parser["metadata_keys_parse"]:
                item.age_present = True
            if "sex" in item.parser["metadata_keys_parse"]:
                item.sex_present = True
    
    with open(library, "r") as file:
        library_data = yaml.safe_load(file)
        library_items = [dict_to_dataclass(SeriesItem, item) for item in library_data["items"]]
        
        # Set age and sex presence fields for library items
        for item in library_items:
            item.age_present = True
            item.sex_present = True

        # Add GEO query URL if not specified
        for item in library_items:
            if item.query is None:
                item.query = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={item.id}"

    return merge_libraries(autoscan_library_items, library_items)


def dict_to_dataclass(cls: Type[Any], data: dict) -> Any:
    field_names = {f.name for f in fields(cls)}
    filtered_data = {k: v for k, v in data.items() if k in field_names}
    return cls(**filtered_data)


def map_repsonse_series_page(page: dict) -> SeriesPage:
    page_id = page["id"]
    created_time = page["created_time"]
    last_edited_time = page["last_edited_time"]
    parent = page["parent"]
    archived = page["archived"]
    in_trash = page["in_trash"]
    url = page["url"]
    properties = create_series_page_properties(page["properties"])

    return SeriesPage(
        page_id=page_id,
        created_time=created_time,
        last_edited_time=last_edited_time,
        parent=parent,
        archived=archived,
        in_trash=in_trash,
        url=url,
        properties=properties,
    )


def get_data_to_add_to_notion(
    local_dataset: list[SeriesItem], remote_dataset: list[SeriesPage]
) -> list[dict]:

    remote_series = [
        page.properties.Name.title[0].plain_text for page in remote_dataset
    ]

    list_to_add = []
    for item in local_dataset:
        if item.id not in remote_series:
            list_to_add.append(create_notion_page_creation(item))

    return list_to_add


def create_notion_page_creation(local_item: SeriesItem) -> dict:
    creation = {"parent": {"database_id": database_id}, "properties": {}}

    creation["properties"]["Name"] = {
        "title": [{"text": {"content": local_item.id}}]
    }
    creation["properties"]["Link"] = {"url": local_item.query}
    creation["properties"]["Title"] = {
        "rich_text": [{"text": {"content": local_item.title}}]
    }
    platforms = [
        {"name": platform["name"]} for platform in (local_item.platforms or [])
    ]
    creation["properties"]["Platform"] = {"multi_select": platforms}
    creation["properties"]["Samples"] = {"number": local_item.samples}

    if local_item.tags:
        creation["properties"]["Tags"] = {
            "multi_select": [{"name": tag} for tag in local_item.tags]
        }

    local_age_select = YES if local_item.age_present else NO
    creation["properties"]["AgePresent"] = {
        "select": {"name": local_age_select}
    }

    local_sex_select = YES if local_item.sex_present else NO
    creation["properties"]["SexPresent"] = {
        "select": {"name": local_sex_select}
    }


    return creation


def get_data_to_update_in_notion(
    local_dataset: list[SeriesItem], remote_dataset: list[SeriesPage]
) -> list[dict]:

    updates = []
    for item in local_dataset:
        remote_item = find_series_in_remote_dataset(item.id, remote_dataset)
        if remote_item:
            update = create_notion_page_update(item, remote_item)
            if update:
                updates.append(update)

    return updates


def create_notion_page_update(
    local_item: SeriesItem, remote_item: SeriesPage
) -> dict:

    if local_item.id != remote_item.properties.Name.title[0].plain_text:
        raise ValueError(f"Local item id {local_item.id} does not match remote item id {remote_item.properties.Name.title[0].plain_text}")

    update = {}
    if local_item.query != remote_item.properties.Link.url:
        update["Link"] = {"url": local_item.query}

    local_age_select = YES if local_item.age_present else NO
    remote_age_select = (
        remote_item.properties.AgePresent.select.name 
        if remote_item.properties.AgePresent.select is not None 
        else None
    )
    if local_age_select != remote_age_select:
        update["AgePresent"] = {"select": {"name": local_age_select}}

    local_sex_select = YES if local_item.sex_present else NO
    remote_sex_select = (
        remote_item.properties.SexPresent.select.name 
        if remote_item.properties.SexPresent.select is not None 
        else None
    )
    if local_sex_select != remote_sex_select:
        update["SexPresent"] = {"select": {"name": local_sex_select}}

    notion_platforms = [
        select.name for select in remote_item.properties.Platform.multi_select
    ]
    local_platforms = [platform["name"] for platform in (local_item.platforms or [])]
    if set(local_platforms) != set(notion_platforms):
        update["Platform"] = {
            "multi_select": [
                {"name": platform} for platform in local_platforms
            ]
        }

    if local_item.samples != remote_item.properties.Samples.number:
        update["Samples"] = {"number": local_item.samples}

    remote_title = (remote_item.properties.Title.rich_text[0]["text"]["content"] 
                   if remote_item.properties.Title.rich_text 
                   else None)
    if local_item.title != remote_title:
        update["Title"] = {
            "rich_text": [{"text": {"content": local_item.title}}]
        }

    if local_item.tags:
        notion_tags = [
            select.name for select in remote_item.properties.Tags.multi_select
        ]
        if set(local_item.tags) != set(notion_tags):
            update["Tags"] = {
                "multi_select": [{"name": tag} for tag in local_item.tags]
            }

    if update:
        return {
            "page_id": remote_item.page_id,
            "params": {"properties": update},
        }

    return None


def find_series_in_remote_dataset(
    series: str, remote_dataset: list[SeriesPage]
) -> SeriesPage | None:
    for item in remote_dataset:
        if item.properties.Name.title[0].plain_text == series:
            return item


def get_data_to_delete_in_notion(
    local_dataset: list[SeriesItem], remote_dataset: list[SeriesPage]
) -> list[SeriesPage]:

    local_ids = [item.id for item in local_dataset]

    list_to_delete = []
    for item in remote_dataset:
        series = item.properties.Name.title[0].plain_text
        if series not in local_ids:
            list_to_delete.append(
                {"page_id": item.page_id, "params": {"archived": True}}
            )

    return list_to_delete


async def limited_concurrent_run(tasks, max_concurrent_tasks):
    semaphore = asyncio.Semaphore(max_concurrent_tasks)

    async def sem_task(task):
        async with semaphore:
            return await task

    # Wrap tasks with semaphore control
    wrapped_tasks = [sem_task(task) for task in tasks]
    await asyncio.gather(*wrapped_tasks)


async def sync_library_entry_to_notion():

    series_items = read_series_dataset_from_library()
    series_pages = await get_series_dataset_from_notion()

    list_to_add = get_data_to_add_to_notion(series_items, series_pages)
    creations = [notion.pages.create(**item) for item in list_to_add]
    list_to_update = get_data_to_update_in_notion(series_items, series_pages)
    updates = [
        notion.pages.update(item["page_id"], **item["params"])
        for item in list_to_update
    ]
    list_to_delete = get_data_to_delete_in_notion(series_items, series_pages)
    deletes = [
        notion.pages.update(item["page_id"], **item["params"])
        for item in list_to_delete
    ]

    print(
        f"Has {len(creations)} to add, {len(updates)} to update and {len(deletes)} to delete"
    )

    all_tasks = creations + updates + deletes
    await limited_concurrent_run(all_tasks, 3)


asyncio.run(sync_library_entry_to_notion())
