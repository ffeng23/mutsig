"""
Author: Feng, Feng 

Run SigProfilerMatrixGenerator/SigProfilerExtractor/SigProfilerAssignment on a set of VCF files,
SigProfilerMatrixGenerator to generate mutation matrices, 
SigProfilerExtractor to de-novo extract mutation signatures.
SigProfilerAssignment to assign COSMIC signatures to the samples.

Usage:
    # Use the current directory as the VCF folder (default):
    python run_sigature_generate.py

    # Specify the VCF folder explicitly:
    python run_sigature_generate.py --vcf-dir HN_HPV/HNhpv_vcf
    python run_sigature_generate.py --vcf-dir /absolute/path/to/vcfs

    # Override genome, project name, or signature extraction range:
    python run_sigature_generate.py --vcf-dir HN_HPV/HNhpv_vcf --genome GRCh38 --project MyProject
    python run_sigature_generate.py --vcf-dir HN_HPV/HNhpv_vcf --min-sigs 2 --max-sigs 10 --nmf-reps 50

    # Run matrix generation only (skip extractor):
    python run_sigature_generate.py --vcf-dir HN_HPV/HNhpv_vcf --skip-extraction
"""

import argparse
import logging
from pathlib import Path

from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sigExt
from SigProfilerAssignment import Analyzer as spa

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DEFAULT_GENOME = "GRCh37"   # hg19 — change to GRCh38 / mm10 etc. if needed

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def sanitize_vcfs(vcf_dir: Path) -> None:
    """
    Strip Windows carriage returns (\r) from all VCF files in *vcf_dir*.

    VCF files exported from Windows or transferred without conversion often
    carry CRLF (\r\n) line endings. The \r bleeds into the last field (ALT),
    turning single-base substitutions like 'C' into 'C\r'. This makes every
    SNV appear to have a 2-character ALT, which prevents SigProfilerMatrixGenerator
    from classifying any mutation correctly and suppresses the ID/SBS output.
    """
    vcf_files = list(vcf_dir.glob("*.vcf"))
    fixed = 0
    for vcf in vcf_files:
        text = vcf.read_bytes()
        if b"\r" in text:
            vcf.write_bytes(text.replace(b"\r", b""))
            fixed += 1
    if fixed:
        log.info("Stripped Windows line endings (CRLF) from %d VCF file(s).", fixed)
    else:
        log.info("VCF line endings are clean (no CRLF found).")


def install_genome(genome: str) -> None:
    """Download the reference genome data if it has not been installed yet."""
    log.info("Checking / installing reference genome: %s", genome)
    genInstall.install(genome, rsync=False, bash=True)
    log.info("Reference genome ready.")


def run_matrix_generator(
    project: str,
    input_dir: Path,
    output_dir: Path,
    genome: str,
) -> None:
    """
    Call SigProfilerMatrixGenerator to build SBS / DBS / ID matrices.

    Parameters
    ----------
    project    : short name used for output file prefixes
    input_dir  : directory that contains the VCF files
    output_dir : root directory for all output files
    genome     : reference genome identifier (e.g. 'GRCh37')
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    log.info("Input  VCF directory : %s", input_dir)
    log.info("Output directory     : %s", output_dir)
    log.info("Project name         : %s", project)
    log.info("Reference genome     : %s", genome)

    # SigProfilerMatrixGenerator expects string paths.
    # NOTE: 'vcf', 'tsv_format', and 'project_dir' were removed in v1.3+.
    # Use 'output_directory' to control where results are written.
    # tsb_stat=False  → skip the transcriptional strand-bias (TSB/24/384) matrices.
    # seqInfo=False   → skip per-mutation sequence-context files (keeps output smaller).
    matrices = matGen.SigProfilerMatrixGeneratorFunc(
        project,                      # project name / output file prefix
        genome,                       # reference genome identifier
        str(input_dir),               # directory that contains the VCF files
        exome=False,                  # whole-genome (not exome)
        bed_file=None,                # no region filter
        chrom_based=False,            # genome-wide matrices (not per-chromosome)
        plot=True,                    # write per-sample spectrum plots
        #tsb_stat=True,               # skip TSB (tsb_stat=True triggers a pandas KeyError in v1.3.x)
        #seqInfo=True,                 # write per-mutation seqInfo files
        output_directory=str(output_dir),  # where to save all results
    )

    log.info("Matrix generation complete.")
    log.info("Generated matrices:")
    for matrix_type, df in matrices.items():
        log.info("  %-6s  shape=%s", matrix_type, df.shape)
    return matrices


def run_extractor(
    project: str,
    matrix_output_dir: Path,
    extractor_output_dir: Path,
    genome: str,
    minimum_signatures: int = 1,
    maximum_signatures: int = 10,
    nmf_replicates: int = 100,
) -> None:
    """
    Run SigProfilerExtractor on the SBS96 matrix produced by the matrix generator.

    The extractor performs de-novo NMF decomposition to identify the operative
    mutational signatures and map them to COSMIC reference signatures.

    Parameters
    ----------
    project              : project name (used for display only)
    matrix_output_dir   : directory where SigProfilerMatrixGenerator wrote its results
                          (contains the SBS/ subfolder)
    extractor_output_dir : directory where SigProfilerExtractor will write results
    genome               : reference genome identifier (e.g. 'GRCh37')
    minimum_signatures   : minimum number of NMF signatures to try
    maximum_signatures   : maximum number of NMF signatures to try
    nmf_replicates       : number of NMF replicates per signature count
    """
    # Locate the SBS96 matrix file produced by the matrix generator
    sbs96_file = matrix_output_dir / "SBS" / f"{project}.SBS96.all"
    if not sbs96_file.exists():
        log.warning(
            "SBS96 matrix file not found at %s — skipping signature extraction. "
            "(This happens when the VCF files contain no SNVs.)",
            sbs96_file,
        )
        return

    extractor_output_dir.mkdir(parents=True, exist_ok=True)

    log.info("--- SigProfilerExtractor ---")
    log.info("Input matrix   : %s", sbs96_file)
    log.info("Output dir     : %s", extractor_output_dir)
    log.info("Genome         : %s", genome)
    log.info("Signatures     : %d – %d", minimum_signatures, maximum_signatures)
    log.info("NMF replicates : %d", nmf_replicates)

    sigExt.sigProfilerExtractor(
        input_type="matrix",            # pre-built mutation matrix (not raw VCFs)
        output=str(extractor_output_dir),
        input_data=str(sbs96_file),
        reference_genome=genome,
        opportunity_genome=genome,
        context_type="SBS96",           # extract SBS96 signatures
        exome=False,
        minimum_signatures=minimum_signatures,
        maximum_signatures=maximum_signatures,
        nmf_replicates=nmf_replicates,
        cpu=20,                         # use all available CPUs
        gpu=False,
        make_decomposition_plots=True,  # plot decomposition into COSMIC sigs
        export_probabilities=True,      # write per-mutation signature probabilities
    )

    log.info("Signature extraction complete. Results: %s", extractor_output_dir)


def run_assignment(
    project: str,
    matrix_output_dir: Path,
    assignment_output_dir: Path,
    genome: str,
) -> None:
    """
    Run SigProfilerAssignment on the SBS96 matrix produced by the matrix generator.

    This maps the samples' mutation catalogs directly to known COSMIC signatures.
    """
    # Locate the SBS96 matrix file produced by the matrix generator
    sbs96_file = matrix_output_dir / "SBS" / f"{project}.SBS96.all"
    if not sbs96_file.exists():
        log.warning(
            "SBS96 matrix file not found at %s — skipping signature assignment. "
            "(This happens when the VCF files contain no SNVs.)",
            sbs96_file,
        )
        return

    assignment_output_dir.mkdir(parents=True, exist_ok=True)

    log.info("--- SigProfilerAssignment ---")
    log.info("Input matrix : %s", sbs96_file)
    log.info("Output dir   : %s", assignment_output_dir)
    log.info("Genome       : %s", genome)

    spa.cosmic_fit(
        samples=str(sbs96_file),
        output=str(assignment_output_dir),
        input_type="matrix",
        context_type="96",
        genome_build=genome,
        cosmic_version=3.5,
        make_plots=True,
        collapse_to_SBS96=True,
        export_probabilities=True,
    )

    log.info("Signature assignment complete. Results: %s", assignment_output_dir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate mutation matrices from VCF files with SigProfilerMatrixGenerator, "
            "then extract de-novo signatures with SigProfilerExtractor."
        )
    )
    parser.add_argument(
        "--vcf-dir",
        type=Path,
        default=Path("."),
        metavar="DIR",
        help=(
            "Path to the folder containing VCF files. "
            "Defaults to the current working directory."
        ),
    )
    parser.add_argument(
        "--genome",
        default=DEFAULT_GENOME,
        metavar="GENOME",
        help=(
            "Reference genome identifier, e.g. GRCh37, GRCh38, mm10. "
            f"Defaults to '{DEFAULT_GENOME}'."
        ),
    )
    parser.add_argument(
        "--project",
        default=None,
        metavar="NAME",
        help=(
            "Project name used for output file prefixes. "
            "Defaults to the name of the parent directory of --vcf-dir."
        ),
    )
    parser.add_argument(
        "--min-sigs",
        type=int,
        default=1,
        metavar="N",
        help="Minimum number of signatures for NMF extraction. Default: 1.",
    )
    parser.add_argument(
        "--max-sigs",
        type=int,
        default=10,
        metavar="N",
        help="Maximum number of signatures for NMF extraction. Default: 10.",
    )
    parser.add_argument(
        "--nmf-reps",
        type=int,
        default=100,
        metavar="N",
        help="Number of NMF replicates per signature count. Default: 100.",
    )
    parser.add_argument(
        "--skip-extraction",
        action="store_true",
        help="Skip SigProfilerExtractor (de-novo extraction).",
    )
    parser.add_argument(
        "--skip-assignment",
        action="store_true",
        help="Skip SigProfilerAssignment (known signature fitting).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Resolve the VCF directory to an absolute path
    input_dir = args.vcf_dir.resolve()
    if not input_dir.is_dir():
        raise SystemExit(f"ERROR: VCF directory not found: {input_dir}")

    # Output directory sits next to the VCF folder, suffixed with '_output'
    output_dir = input_dir.parent / (input_dir.name + "_output")

    # Extractor output sits next to the VCF folder, suffixed with '_signatures'
    extractor_dir = input_dir.parent / (input_dir.name + "_signatures")

    # Assignment output sits next to the VCF folder, suffixed with '_assignment'
    assignment_dir = input_dir.parent / (input_dir.name + "_assignment")

    # Project name defaults to the parent folder of the VCF directory
    project = args.project or input_dir.parent.name

    log.info("=== SigProfiler Pipeline — %s ===", project)

    # Step 1: sanitize VCF line endings (CRLF -> LF)
    sanitize_vcfs(input_dir)

    # Step 2: ensure the reference genome is available
    install_genome(args.genome)

    # Step 3: generate the mutation matrices
    run_matrix_generator(
        project=project,
        input_dir=input_dir,
        output_dir=output_dir,
        genome=args.genome,
    )
    log.info("Matrix results saved to: %s", output_dir)

    # Step 4: de-novo signature extraction (optional)
    if args.skip_extraction:
        log.info("--skip-extraction flag set; skipping SigProfilerExtractor.")
    else:
        run_extractor(
            project=project,
            matrix_output_dir=output_dir,
            extractor_output_dir=extractor_dir,
            genome=args.genome,
            minimum_signatures=args.min_sigs,
            maximum_signatures=args.max_sigs,
            nmf_replicates=args.nmf_reps,
        )

    # Step 5: signature assignment (optional)
    if args.skip_assignment:
        log.info("--skip-assignment flag set; skipping SigProfilerAssignment.")
    else:
        run_assignment(
            project=project,
            matrix_output_dir=output_dir,
            assignment_output_dir=assignment_dir,
            genome=args.genome,
        )

    log.info("=== Pipeline complete ===")


if __name__ == "__main__":
    main()
