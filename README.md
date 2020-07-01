# annotateVcf
| Master                                              | Develop                                               |
| --------------------------------------------------- | ----------------------------------------------------- |
| [![Master Badge][travis-master-badge]][travis-repo] | [![Develop Badge][travis-develop-badge]][travis-repo] |

This project hosts scripts to annotate VCF files using user defined driver genes and mutations

<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Design](#design)
- [Tools](#tools)
	- [annotateVcf](#annotateVcf)
	- [inputFormat](#inputformat)
	- [outputFormat](#outputformat)
- [INSTALL](#install)
	- [Package Dependencies](#package-dependencies)
- [Development environment](#development-environment)
	- [Development Dependencies](#development-dependencies)
		- [Setup VirtualEnv](#setup-virtualenv)
- [Cutting a release](#cutting-a-release)
	- [Install via `.whl` (wheel)](#install-via-whl-wheel)
- [Reference](#reference)

<!-- /TOC -->

## Design
Uses bcftools and tabix from htslib

## Tools

`annotateVcf` has multiple commands, listed with `annotateVcf --help`.

### annotateVcf
Takes vcf file as input along with driver gene information and outputs VCF with updated INFO field.

Various exceptions can occur for malformed input files.

### inputFormat

 * ```VCF file```  snp or indel vcf file annoated using [VAGrENT]
 * ```Driver annotation and reference files``` please refer the files in annotate/config folder
 * ```lof_genes.txt``` list of known loss of function [LoF] genes
 * ```driver_mutations_grch38.tsv.gz``` tab separated driver mutations along with effect type
 * ```genome_grch38.tab.gz``` tab separated chromosme lenght file created from genome index file
 * ```lof_consequences.txt``` consequence types to be considerd under LoF category
 * ```vcf_info.txt``` vcf header and INFO tag

### outputFormat

 * ```input.drv.vcf.gz ``` output vcf file with DRV info tag 


## INSTALL
Installing via `pip install`. Simply execute with the path to the compiled 'whl' found on the [release page][annotateVcf-releases]:

```bash
pip install annotateVcf.X.X.X-py3-none-any.whl
```

Release `.whl` files are generated as part of the release process and can be found on the [release page][annotateVcf-releases]

### Package Dependancies

`pip` will install the relevant dependancies, listed here for convenience, please refer requirements.txt for versions:

## Development environment

This project uses git pre-commit hooks.  As these will execute on your system it
is entirely up to you if you activate them.

If you want tests, coverage reports and lint-ing to automatically execute before
a commit you can activate them by running:

```
git config core.hooksPath git-hooks
```

Only a test failure will block a commit, lint-ing is not enforced (but please consider
following the guidance).

You can run the same checks manually without a commit by executing the following
in the base of the clone:

```bash
./run_tests.sh
```

### Development Dependencies

#### Setup VirtualEnv

```
cd $PROJECTROOT
hash virtualenv || pip3 install virtualenv
virtualenv -p python3 env
source env/bin/activate
python setup.py develop # so bin scripts can find module
```

For testing/coverage (`./run_tests.sh`)

```
source env/bin/activate # if not already in env
pip install pytest
pip install pytest-cov
```

__Also see__ [Package Dependancies](#package-dependancies)

### Cutting a release

__Make sure the version is incremented__ in `./setup.py`

### Install via `.whl` (wheel)

Generate `.whl`

```bash
source env/bin/activate # if not already
python setup.py bdist_wheel -d dist
```

Install .whl

```bash
# this creates an wheel archive which can be copied to a deployment location, e.g.
scp dist/annotateVcf.X.X.X-py3-none-any.whl user@host:~/wheels
# on host
pip install --find-links=~/wheels annotateVcf
```

### Reference
<!--refs-->
 [VAGrENT]: https://github.com/cancerit/VAGrENT 
 [travis-master-badge]: https://travis-ci.org/cancerit/annotateVcf.svg?branch=master
 [travis-develop-badge]: https://travis-ci.org/cancerit/annotateVcf.svg?branch=develop
 [travis-repo]: https://travis-ci.org/cancerit/annotateVcf
 [annotateVcf-releases]: https://github.com/cancerit/annotateVcf/releases
