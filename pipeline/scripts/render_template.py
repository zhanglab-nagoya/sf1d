from jinja2 import Environment, FileSystemLoader
from pathlib import Path

# Take template path from snakemake input
template_path = Path(snakemake.input[0])

template_dir = template_path.parent
template_name = template_path.name

# Set up Jinja environment
env = Environment(
    loader=FileSystemLoader(str(template_dir)),
    autoescape=False,
)

template = env.get_template(template_name)
rendered = template.render(
    wildcards=snakemake.wildcards,
    params=snakemake.params,
    input=snakemake.input,
    output=snakemake.output,
    config=snakemake.config,
)

with open(snakemake.output[0], "w") as f:
    f.write(rendered)
