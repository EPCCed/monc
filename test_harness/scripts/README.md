# Working with the test_harness Scripts.

The Scripts require a Python environment.  On ARCHER2, this was configured using a local `uv` venv:

```
curl -LsSf https://astral.sh/uv/install.sh | sh
```

This was installed into `$HOME/.local/bin`, which is added to the $PATH.

A uv environment was created as: `uv venv mth --python 3.13`  (mth for MONC test_harness...), and activated: `source mth/bin/activate`.

Necessary libraries were installed into the active env: `uv pip install netCDF4 matplotlib xarray[complete]`


---

## Test Harness Concept:

1. Select a set of simulation cases
2. Run these for two processor number specifications (and ultimately MONC grid decompositions, given the same domain): 36 and 72 cores
3. Compare to a previous execution of the same simulation cases to allow inspection of and reports on:
  - bitwise comparability (where this is expected across decomposition--not where domain means are being considered)
  - physical property consistency

Within the top-level `test_harness` directory, you will find:
- test suite configurations - collections of MONC `.mcf` specification cases
  - `casim_aerosol_processing`
  - `monc_casim_socrates`
  - `monc_ecse`
  - `monc_main`
- `submonc_template.pbs` and `submonc_template.sb` - scheduler-dependent submission script templates
- `continuation.sh` - a test_harness-specific copy of the continuation script to enable automatic sub-cycling between checkpoints

---

## Operation

Tend to begin by making the scripts visible on the path and activating the env:
```
    export PATH=$PATH:..../monc/test_harness/scripts
    source ...../mth/bin/activate
```

And think about the two main tools:
1. `monc_test_harness.py` - configure `test_harness` suite, build model, run `test_harness` jobs
2. `ck_progress.sh` - inspect the status of submitted work (also can submit the jobs)
3. `monc_kgo_bit_compare.py` - bitwise comparison of two `test_harness` runs


### `monc_test_harness.py`

**Execute from the MONC top-level directory.**

This script will create new directories under `test_harness/`, named to match your specifications, beginning with the name of the system on which the `test_harness` is being run.  This will contain these directories:
- `bit_comparison` - to receive output from `monc_kgo_bit_compare.py`, comparing this execution to a previous execution
- `build_dir` - to store the corresponding executable, preprocessed code, and output of the build steps  (**when built with this script**)
- `chk_point_dump` - to store checkpoints created during execution
- `monc_stdout` - to store stdout created during execution
- `ncfiles` - to store diagnostic output files created during execution
- `pbs_dir` - to store `*.o*` files created for each job (name is historical)
- `tmp_config_dir` - to store case configuration files
- `tmp_qsub_dir` - to store submission scripts (name is historical)

**Note that loaded modules are only configured for the listed known sets.**  The system would need to be modified to work directly with others.  This is similarly true for use on other systems.
- Bypassing this, one could build the model separately, and this system will copy the existing build information and executable already present in the MONC main directories.  However, automatic naming conventions could be a problem, since this only expects certain values.


```
Options:
    -t or --test_suite,    DEFAULT: main_component
        the suite type, i.e. standard, main_component, ecse, casim_socrates, or casim_aerosol_processing
        N.B.: The standard suite uses existing directories in the local branch, while
              main_component and casim_socrates use new directories under test_harness.
              These have not yet been configured to run under Slurm (ARCHER2).

    -b or --build,         DEFAULT: False
        whether to build MONC, True or False

    -r or --run,           DEFAULT: False
        whether to run the test harness, True or False

    -m or --modules,       DEFAULT: default
        the modules to load, i.e. default or cdt, where cdt is latest modules from
        Cray Development tools (Monsoon ONLY). Also, on ARCHER2: md1.0, md2.0, 2.0-cce20

    -c or --compiler,      DEFAULT: cray
        the compiler, i.e. cray or gnu

    -o or --opt_level,     DEFAULT: high
        the level of optimisation and compiler flags, i.e. debug, safe, high

    -s or --system,        DEFAULT: MO-XC40
        the system to run the test suite on, MO-XC40 or ARCHER2

    -a or --account,       DEFAULT: ""
        REQUIRED for ARCHER2, the name of the project account to charge

    -l or --label,         DEFAULT: ""
        Completely optional label appended to test_harness run directory.  Useful for repeated tests with
        the same test_harness options, but with potentially different model options or simple retests.

    -i or --input_option,  DEFAULT: ""
        A list of option/value pairs to be added to the configuration files for all cases in this
        test_suite.  These are given as a single string, beginning with and delimited by the ? symbol.
        For example, we can ensure that threading is disabled in the IO server and that the flux_budget
        component is disabled in all cases by supplying:
            -i  ?l_thoff=.true.?flux_budget_enabled=.false.
        Any number of these can be supplied.  They are not checked for validity until used by the model.
        Quotation marks are not retained when passed in this way.

default settings =
    monc_test_harness.py -t main_component -b False -r False -m default -c cray -o high -s MO-XC40
```

### `ck_progress.sh`

This is a tool used to check on all (or a subset) existing `test_harness` configurations.  It will check the running status and stdout files for a few kinds of errors and **optionally** submits jobs that have not been submitted or resubmits those that have encountered errors and could benefit from restart.

This can be run repeatedly for existing configurations until you have success on all of them.

```
#  Input (all optional):
#    -a  Perform submission and restart actions.
#          DEFAULT: Not enabled. Only checks on status and logs success.
#    -c  Provide a single text pattern that will be globbed as:
#              test_harness/*${pattern}*
#        to select configuration directories for evaluation.
#          DEFAULT: All existing configuration directories are evaluated
#    -r  Number of restarts to allow for each job.
#          DEFAULT: 2
#    -o  Don't ask to confirm settings.
#          DEFAULT: Confirm settings with user.
#
#
#  Output:
#    stdout - print to screen job status ['still running', 'ATTENTION', or 'Success']
#           - tail of most recent cycle stdout
#    copies of stdout from failed jobs as case.rerun in test_harness/.../monc_stdout
#    empty case.success files for completed cases as case.success in test_harness/.../monc_stdout
#
#  Actions:
#    Does nothing for running jobs
#    Marks successful jobs with a .success file in the monc_stdout directory
#      - .success files will contain any 'Miss match or un-enabled' warnings from the
#        final submission cycle
#    Attempts to restart failed jobs (two attpemts only, unless changed with -r)
#      - on second attempt, the previous checkpoint is deleted to restart further back
#      - If you wish to try more times, you can clear the contents of the offending .rerun file,
#        being careful to clear checkpoints as needed, too.
#
#  Example usage:
#    To run on all cases without submission actions:
#      ck_progress.sh
#    To run with submission and restart actions enabled for a very specific configuration
#    directory pattern:
#      ck_progress.sh -a -c ARCHER2_main_component_gnu_high
```

### `monc_kgo_bit_compare.py`

This script compares two (probably) completed executions of the test_harness: `kgo` vs. `new` .  Typically, we are comparing two versions of code that have been compiled under the same settings, but this is not strictly required or at all checked.  That is, compiler and optimisation settings are requested for reporting and identifying the **new** simulation data.

The system needs to know the correct test suite cases because specific handling varies between test cases (e.g., times, variables inspected).  

Otherwise, point it at the `test_harness/< test harness configuration >/ncfiles` directory for both `kgo` and `new` data.  Where the paths are not provided, the user will be prompted to enter this information anyway.

Typically, we like to run the comparison test from the `test_harness/<master_dir>/bit_comparison` directory for the `new` data.

```
Options:
    -t or --test_suite,  DEFAULT: ""
        the suite type, i.e. standard, main_component, ecse, casim_socrates, casim_aerosol_processing

    -c or --compiler     DEFAULT: ""
        the compiler, i.e. cray or gnu

    -o or --opt_level    DEFAULT: ""
        optimisation level, i.e. high, safe, debug

    -k or --kgo_path,    DEFAULT: undefined
        absolute or relative path to KGO data directory, e.g., ../ncfiles
        if not provided, user will be prompted to enter this information

    -n or --new_path,    DEFAULT: undefined
        absolute or relative path to New data directory, e.g., ../ncfiles
        if not provided, user will be prompted to enter this information

NOTE: the compiler and opt_level options are only printed to the output file. They are not checked.
```

---

## Comparison Output

In all cases, `bit_compare_results.txt` will be produced.  This is a textual summary of the comparison results.   We will be looking for successes and failures in bit comparison.  Comparisons are made between:
1. 36_pe vs. 72_pe (primarily on the `new` case, but also with reference to the `kgo`)
2. `new` vs. `kgo` for both processor decompositions

Specific MONC prognostic fields (u, v, w, q_vapour, q_cloud_liquid_mass) are inspected for each relevant simulation case.
E.g.:

```
****************Bit comparison results from  main_component  component testing**************
opt_level: safe
compiler:  cray
Known good MONC output = /mnt/lustre/a2fs-work4/work/ecsegb12/ecsegb12/toddjg/EPCCed/monc/test_harness/ARCHER2_main_component_cray_safe_vn1/ncfiles/
New MONC output = /mnt/lustre/a2fs-work4/work/ecsegb12/ecsegb12/toddjg/EPCCed/monc/test_harness/ARCHER2_main_component_cray_md2.0_safe_vn1/ncfiles/
********************************************************************************************
******************************************bubble******************************************
==========================================================================================
  Case:  WarmPw_  [decomp success expected]
Field = u
Decomposition comparison:
Field = u - New 36_pe vs 72_pe WarmPw_ SUCCESS (KGO also Passed)
New vs KGO comparison:
Field = u - New 36_pe vs kgo 36_pe WarmPw_ SUCCESS
Field = u - New 72_pe vs kgo 72_pe WarmPw_ SUCCESS
..........................................................................................
Field = v
Decomposition comparison:
Field = v - New 36_pe vs 72_pe WarmPw_ SUCCESS (KGO also Passed)
New vs KGO comparison:
Field = v - New 36_pe vs kgo 36_pe WarmPw_ SUCCESS
Field = v - New 72_pe vs kgo 72_pe WarmPw_ SUCCESS
..........................................................................................
Field = w
Decomposition comparison:
Field = w - New 36_pe vs 72_pe WarmPw_ SUCCESS (KGO also Passed)
New vs KGO comparison:
Field = w - New 36_pe vs kgo 36_pe WarmPw_ SUCCESS
Field = w - New 72_pe vs kgo 72_pe WarmPw_ SUCCESS
..........................................................................................
```

Where unexpected differences are found, `.png` figure files will be generated for later manual inspection.  Inspection will help to identify the severity of the discovered differences and whether these have meaningfully changed the simulation from a physical perspective. 

Alternatively, there is an unadvertised option to produce figures for all possible comparisons: `-a True` to "plot all".


