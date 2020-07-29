from box import Box

ACTIVATE_POREC = "set +u; source ~/miniconda3/etc/profile.d/conda.sh ; conda activate ; conda activate poreC; "

DASK_SETTINGS = "--dask-scheduler-port 0"


def create_path_accessor(prefix: Path) -> Box:
    """Create a Box to provide '.' access to hierarchy of paths"""
    data = yaml.load(Path("file_layout.yaml").open(), Loader=yaml.SafeLoader)
    paths = {}
    for directory in data.keys():
        paths[directory] = {}
        for file_alias, file_name in data[directory].items():
            p = str(prefix / directory / file_name)
            if "{refgenome_id}" in p:
                p = p.replace("{refgenome_id}", config["refgenome"]["refgenome_id"])
            paths[directory][file_alias] = str(p)
    return Box(paths, frozen_box=True)


def to_log(path: str) -> str:
    """Log file location based on output file"""
    return str(outdir / "logs" / path) + ".log"


def to_benchmark(path: str) -> str:
    """Log file location based on output file"""
    return str(outdir / "benchmarks" / path) + ".bench.txt"


def to_prefix(path: str, components=2) -> str:
    """Strip trailing extensions to create an output prefix"""
    return path.rsplit(".", components)[0]


def expand_rows(path: str, df: pd.DataFrame):
    """Expand a templated path string with values from a dataframe"""
    res = df.apply(lambda x: path.format(**x.to_dict()), axis=1)
    return list(res)


def lookup_value(column, df):
    """Use wildcards to 'lookup' a value in a dataframe. The wildcard keys must
    match the dataframe index.
    """
    index_names = tuple(df.index.names)
    assert column in df.columns

    def _inner(wildcards):
        row = df.xs(tuple(wildcards[k] for k in index_names), level=index_names, drop_level=True)
        assert len(row) == 1
        return row[column].values[0]

    return _inner
