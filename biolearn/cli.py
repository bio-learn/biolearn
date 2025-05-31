"""Command-line interface for biolearn."""

import click

from biolearn.metadata import search_metadata


@click.group()
def cli():
    """Biolearn command-line interface."""


@cli.command("search-metadata")
@click.option("--field", required=True, help="Metadata key (e.g., sex)")
@click.option("--value", help="Exact match for --field")
@click.option(
    "--min-age", type=float, default=None, help="Minimum age in years"
)
def _search_metadata(field, value, min_age):
    """
    Query library.yaml without downloading big matrices.

    Examples
    --------
    biolearn search-metadata --field sex --value male
    biolearn search-metadata --field sex --value female --min-age 70
    """
    if field == "age" and min_age:
        df = search_metadata(min_age=min_age)
    else:
        df = search_metadata(**{field: value})
    click.echo(df.to_string(index=False))
