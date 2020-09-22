from box import Box

ACTIVATE_POREC = "set +u; source ~/miniconda3/etc/profile.d/conda.sh ; conda activate ; conda activate poreC; "

DASK_SETTINGS = "--dask-scheduler-port 0"
WRAPPER_PREFIX = f"file:{workflow.basedir}/wrappers"
OUTDIR = Path(config['output_dir'])



if config["pore_c_version"] == "rel":
    PORE_C_CONDA_FILE = "../envs/pore_c.yml"
else:
    # if running in development mode then we expect pore_c to be installed in the
    # environment that contains the snakemake command
    assert config["pore_c_version"] == "dev"
    PORE_C_CONDA_FILE = ""




def create_config_dataframes():
    # one row per reference genome
    reference_df = pd.read_table(config["references"], comment="#").set_index(["refgenome_id"], drop=False)

    # one row per input pore-c fastq
    basecall_df = pd.read_table(config["basecalls"], comment="#").set_index(["run_id", "enzyme"], drop=False)

    #easier for dataentry to have list of refgenomes on basecall table, here we stack so that there's one row per (basecall, regenome) pair
    run_to_ref_df = basecall_df['refgenome_ids'].str.split(',', expand=True).stack().rename("refgenome_id").to_frame().droplevel(-1, axis=0)
    # now join this to the refgenome table
    mapping_df = basecall_df.join(run_to_ref_df).drop(["refgenome_ids"], axis=1)

    if Path(config['phased_vcfs']).exists():
        # if we have a phased VCF for some biospecimen, regenome combination we add it to the mapping table
        phased_vcf_df = pd.read_table(config['phased_vcfs'], comment="#").set_index(["refgenome_id", "biospecimen"])
        mapping_df = (
            pd.merge(mapping_df, phased_vcf_df, left_on=["refgenome_id", "biospecimen"], right_index=True, how="left")
            .fillna({"phase_set_id": "unphased", "vcf_path": ""})
        )
    else:
        mapping_df['phase_set_id'] = 'unphased'
        mapping_df['vcf_path'] = ''

    mapping_df = mapping_df.set_index(['refgenome_id', 'phase_set_id'], append=True, drop=False)

    # can subset the mapping
    if config['mapping_query']:
        mapping_df = mapping_df.query(config['mapping_query'])
    return basecall_df, reference_df, mapping_df


def create_path_accessor(prefix: Path = OUTDIR) -> Box:
    """Create a Box to provide '.' access to hierarchy of paths"""
    data = yaml.load(Path(config['file_layout']).open(), Loader=yaml.SafeLoader)
    paths = {}
    for directory in data.keys():
        paths[directory] = {}
        for file_alias, file_name in data[directory].items():
            p = str(prefix / directory / file_name)
            paths[directory][file_alias] = str(p)
    return Box(paths, frozen_box=True)


def to_log(path: str) -> str:
    """Log file location based on output file"""
    return str(OUTDIR / "logs" / path) + ".log"


def to_benchmark(path: str) -> str:
    """Log file location based on output file"""
    return str(OUTDIR / "benchmarks" / path) + ".bench.txt"


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
    assert column in df.columns, column

    def _inner(wildcards):
        #print(df)
        if df.index.nlevels == 1:
            res = df.loc[wildcards[index_names[0]], :]
            return res[column]
        else:
            row = df.xs(tuple(wildcards[k] for k in index_names), level=index_names, drop_level=True)
            assert len(row) == 1, row
            return row[column].values[0]

    return _inner

def lookup_json(columns, df):
    """Use wildcards to 'lookup' a value in a dataframe. The wildcard keys must
    match the dataframe index. Return the row as single-line json string
    """
    index_names = tuple(df.index.names)

    for column in df.columns:
        assert column in df.columns, column

    def _inner(wildcards):
        #print(df)
        if df.index.nlevels == 1:
            row = df.loc[wildcards[index_names[0]], :]
        else:
            row = df.xs(tuple(wildcards[k] for k in index_names), level=index_names, drop_level=True)
        assert len(row) == 1, row
        return row.iloc[0,:][columns].to_json()
    return _inner
