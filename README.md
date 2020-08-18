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

Uses [bcftools], [tabix] and [bgzip] in user's path , these are part of [htslib] or can be installed separately

## Tools

`annotateVcf` has multiple command line options, listed with `annotateVcf --help`.

### annotateVcf
Takes vcf file as input along with driver gene information, and optional unmatched normal panel vcf and outputs VCF with added  DRV INFO field.

Various exceptions can occur for malformed input files.

### inputFormat

 * ```input_vcf.gz```  snv or indel vcf file annotated using [VAGrENT]
 * ```normal_panel.vcf.gz```  normal panel to tag germline variants [VAGrENT]
 * ```lof_genes.txt ``` list of known loss of function [LoF] genes along with previous gene symbols ( to make sure all gene synonyms were matched with input vcf)
 * ```consequence types``` lof consequence types to restrict driver gene annotations  
 * ```driver_mutations.tsv.gz``` tab separated driver mutations along with consequence type 
 * ```info.header``` vcf header INFO line 

### outputFormat

 * ```<input>.drv.vcf.gz ``` output vcf file with DRV info field and consequence type if known, LoF in case annotated using LoF gene list.

## INSTALL
Installing via `pip install`. Simply execute with the path to the compiled 'whl' found on the [release page][annotateVcf-releases]:

```bash
pip install annotateVcf.X.X.X-py3-none-any.whl
```

Release `.whl` files are generated as part of the release process and can be found on the [release page][annotateVcf-releases]

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

pytest
radon
pytest-cov

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
pip install radon
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
 [htslib]: https://github.com/samtools/htslib
 [bcftools]: https://github.com/samtools/bcftools
 [tabix]: https://github.com/samtools/tabix
 [VAGrENT]: https://github.com/cancerit/VAGrENT 
 [travis-master-badge]: https://travis-ci.org/cancerit/annotateVCF.svg?branch=master
 [travis-develop-badge]: https://travis-ci.org/cancerit/annotateVCF.svg?branch=develop
 [travis-repo]: https://travis-ci.org/cancerit/annotateVCF
 [annotateVcf-releases]: https://github.com/cancerit/annotateVCF/releases
