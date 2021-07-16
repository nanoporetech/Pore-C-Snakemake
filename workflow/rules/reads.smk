checkpoint import_basecalls:
    output:
        paths.basecall.catalog,
        paths.basecall.summary,
        paths.basecall.read_metadata,
    params:
        fname=lookup_value("fastq_path", basecall_df),
        prefix=to_prefix(paths.basecall.catalog),
        dask_settings=config["software"]["dask"]["settings"],
        max_read_length=config["max_read_length"],
        reads_per_batch=config["reads_per_batch"],
    log:
        to_log(paths.basecall.catalog),
    benchmark:
        to_benchmark(paths.basecall.catalog)
    conda:
        PORE_C_CONDA_FILE
    threads: 1
    shell:
        "pore_c {params.dask_settings} --dask-num-workers {threads} "
        "reads prepare {params.fname} {params.prefix} --max-read-length {params.max_read_length} "
        " --batch-size {params.reads_per_batch} 2> {log}"


rule import_fast5s:
    # Makes symlinks to the fast5 files, which can then be renamed for
    # for f5c index without modifying the original fast5s.
    output:
        directory(paths.fast5.fast5),
    params:
        fname=lookup_value("fast5_directory", basecall_df),
    run:
        expand_tilde = os.path.expanduser(params.fname)
        path_is_absolute = os.path.isabs(expand_tilde)
        # cp -rs requires an absolute path for the source
        if path_is_absolute:
            shell("cp -rs {params.fname} {output}")
        else:
            shell('cp -rs "$(pwd)/{params.fname}" {output}')


rule import_sequencing_summary:
    output:
        paths.fast5.seq_summary,
    params:
        fname=lookup_value("sequencing_summary_path", basecall_df),
    shell:
        """
        cp {params.fname} {output}
        """


rule adjust_pass_fail:
    # Addresses 3 cases re: pass and fail fast5s
    # 1) Older sequencing summary files do not specify pass or fail
    # and there are fast5s with the same name in both the pass and fail folders.
    # This changes the summary to specify ".fail.fast5" or ".pass.fast5" and
    # renames the symlinks according to the folder
    # 2) Newer sequencing summary files have "_pass_" or "_fail_" in their names
    # and will not be changed.
    # 3) Users may have disabled splitting pass/fail fast5s into separate folders.
    # The fast5 directory will not have 'fast5_pass' or 'fast5_fail', and the
    # filenames are already unique, and will not be changed.
    input:
        fast5=ancient(paths.fast5.fast5),
        summary=paths.fast5.seq_summary,
    output:
        modified_summary=paths.fast5.seq_summary_pf,
    run:
        def fast5_are_not_split(directory):
            folders = []
            for split_type in ["*fast5_pass*", "*fast5_fail*", "*fast5_skip*"]:
                folders.extend(list(Path(directory).rglob(split_type)))
            return len(folders) == 0


        seq_sum = pd.read_table(input.summary)
        print(seq_sum.head(1).T)
        if "filename_fast5" in seq_sum.columns:
            fast5_column = "filename_fast5"
        elif "filename" in seq_sum.columns:
            fast5_column = "filename"
        else:
            raise ValueError("Neither 'filename' nor 'filename_fast5' columns found in summary file")

        example_fn = seq_sum[fast5_column][0]
        if any(tag in example_fn for tag in ["_pass_", "_fail_", "_skip_"]):
            print("Pass/fail files already uniquely named; No adjustments needed")
            shell("cp {input.summary} {output.modified_summary}")
        elif fast5_are_not_split(input.fast5):
            print("Fast5 files were not split into pass/fail folders; " + "No adjustments needed")
            shell("cp {input.summary} {output.modified_summary}")
        else:
            seq_sum.loc[seq_sum.passes_filtering == True, fast5_column] = seq_sum[fast5_column].replace(
                "\.fast5$", ".pass.fast5", regex=True
            )
            seq_sum.loc[seq_sum.passes_filtering == False, fast5_column] = seq_sum[fast5_column].replace(
                "\.fast5$", ".fail.fast5", regex=True
            )
            seq_sum.to_csv(output.modified_summary, sep="\t", index=False)

            def adjust_filenames_command(passes_filtering):
                command = (
                    "find "
                    + input.fast5
                    + " -type l -name '*.fast5' "
                    + "| grep fast5_{status} "
                    + "| grep -v -e '\.{status}\.fast5$' -e '_{status}_' "
                    + "| xargs rename 's!\.fast5$!\.{status}\.fast5!' ;"
                ).format(status=passes_filtering)
                return command

            shell(adjust_filenames_command("pass"))
            shell(adjust_filenames_command("fail"))
