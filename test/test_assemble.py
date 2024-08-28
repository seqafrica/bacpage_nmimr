import shutil
from pathlib import Path

import pandas as pd
import pytest
from snakemake import WorkflowError
from snakemake.utils import validate

from bacpage.src import assemble, common_funcs


def test_error_if_project_not_found():
    with pytest.raises( AssertionError ):
        assemble.run_assemble(
            project_directory="/foo/bar",
            configfile=".",
            sample_data=".",
            denovo=False,
            qc=True,
            threads=1
        )


def test_exit_if_zero_threads_allocated():
    with pytest.raises( SystemExit ) as excinfo:
        assemble.run_assemble(
            project_directory="test/test_pipeline",
            configfile=".",
            sample_data=".",
            denovo=False,
            qc=True,
            threads=0
        )


def test_error_if_config_not_automatically_found():
    with pytest.raises( AssertionError ):
        common_funcs.load_configfile( ".", Path( "test/test_fastqs" ).absolute() )


def test_automatically_find_config_file():
    config = common_funcs.load_configfile( ".", Path( "test/test_pipeline" ).absolute() )
    assert config, "An object should be returned but nothing was."


def test_config_is_valid_yaml():
    config = common_funcs.load_configfile( ".", Path( "test/test_pipeline" ).absolute() )
    assert isinstance( config, dict ), f"Returned config is a {type( config )}, expected <class 'dict'>."
    validate( config, "bacpage/schemas/Illumina_config.schema.yaml" )


def test_config_finds_local_parameters():
    config = common_funcs.load_configfile(
        "test/configs/local_parameters.yaml", Path( "test/test_pipeline" ).absolute(),
    )
    reference = Path( config["reference"] )
    assert reference.is_absolute() and reference.exists(), f"Local reference in resources directory was not correctly found ({reference})."


def test_abridged_config_accepted():
    search_directory = "test/test_tree_project_directory"
    config = "test/configs/assemble_only.yaml"
    config, snakefile = assemble.run_assemble(
        project_directory=search_directory,
        configfile=config,
        sample_data=".",
        denovo=False,
        qc=False,
        threads=1,
        dryrun=True
    )

    assert config is not None


def test_error_if_sample_data_not_automatically_found():
    with pytest.raises( SystemExit ) as excinfo:
        assemble.load_sampledata( ".", Path( "test/test_fastqs" ).absolute() )
    assert excinfo.value.code < 0


def test_automatically_find_sample_data():
    sample_data, skipped_samples = assemble.load_sampledata( ".", Path( "test/test_pipeline" ).absolute() )
    assert sample_data is not None, "An object should be returned but nothing was."


def test_sample_data_is_valid_dataframe():
    sample_data, skipped_samples = assemble.load_sampledata( ".", Path( "test/test_pipeline" ).absolute() )
    assert isinstance(
        sample_data, pd.DataFrame
    ), f"Returned sample data is {type( sample_data )}, expected <class 'pd.DataFrame'>."
    validate( sample_data, "bacpage/schemas/Illumina_metadata.schema.yaml" )


def test_duplicate_samples_causes_exit():
    with pytest.raises( SystemExit ) as excinfo:
        sample_data, skipped_samples = assemble.load_sampledata(
            "test/sample_datas/duplicate_names.csv", Path( "." )
        )
    assert excinfo.value.code == -2


def test_number_name_is_accepted():
    sample_data, skipped_samples = assemble.load_sampledata( "test/sample_datas/number_name.csv", Path( "." ) )


def test_error_when_sample_name_contains_illegal_characters():
    with pytest.raises( WorkflowError ) as excinfo:
        sample_data, skipped_samples = assemble.load_sampledata(
            "test/sample_datas/illegal_characters.csv", Path( "." )
        )


@pytest.mark.slow
def test_assemble_snakemake_runs_correctly():
    project_directory = Path( "test/test_pipeline/" )
    expected_output = ["results/consensus/test.consensus.fasta", "results/reports/depth/test.depth.pdf",
                       "results/reports/qc_report.html"]
    expected_output = [project_directory / file for file in expected_output]

    assemble.run_assemble(
        project_directory=str( project_directory ),
        configfile=".",
        sample_data=".",
        denovo=False,
        qc=True,
        threads=-1,
        verbose=False
    )
    results = dict()
    for file in expected_output:
        results[file] = file.exists()
    assert all( results.values() ), f"Not all expected files generated. Expected {results}"

    if (project_directory / "results").exists():
        shutil.rmtree( project_directory / "results" )
    if (project_directory / "intermediates").exists():
        shutil.rmtree( project_directory / "intermediates" )


def test_assemble_postamble():
    project_directory = Path( "test/test_pipeline/" ).absolute()
    assemble.postamble( False, project_directory )


@pytest.mark.slow
def test_denono_assembly_snakemake_runs_correctly():
    project_directory = Path( "test/test_pipeline/" )
    expected_output = ["results/assembly/test.annotated.gff", "results/reports/qc_report.html"]
    expected_output = [project_directory / file for file in expected_output]

    assemble.run_assemble(
        project_directory=str( project_directory ),
        configfile=".",
        sample_data=".",
        denovo=True,
        qc=True,
        threads=8
    )
    results = dict()
    for file in expected_output:
        results[file] = file.exists()
    assert all( results.values() ), f"Not all expected files generated. Expected {results}"

    if (project_directory / "results").exists():
        shutil.rmtree( project_directory / "results" )
    if (project_directory / "intermediates").exists():
        shutil.rmtree( project_directory / "intermediates" )

# def test_find_all_genes():
#    genes_found = list( assemble.get_genes( "resources/cholera_ref_genes/" ).keys() )
#    genes_expected = ["ctxA", "tcpA_classical", "tcpA_eltor", "toxR", "wbeO1", "wbfO139"]
#    assert sorted( genes_found ) == sorted(
#        genes_expected
#    ), f"Incorrect genes found. Expected {genes_expected}, got {genes_found}."
#
#
# def test_error_on_empty_directory():
#    with pytest.raises( AssertionError ):
#        genes_found = list( assemble.get_genes( "test" ).keys() )
#
