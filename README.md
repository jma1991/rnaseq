# rnaseq <img align="right" width="200" src="images/shield.png">

A workflow for RNA-seq analysis in Snakemake.

[![Snakemake][shield-snakemake]](https://snakemake.readthedocs.io)
[![MIT license][shield-license]](https://choosealicense.com/licenses/mit)

Table of Contents
-----------------

  * [Introduction](#introduction)
  * [Requirements](#requirements)
  * [Usage](#usage)
  * [Contributing](#contributing)
  * [Support and Migration](#support-and-migration)
  * [Thanks](#thanks)
  * [License](#license)

Introduction
------------

The Ripple workflow is a bioinformatics analysis pipeline for RNA sequencing data. The workflow is built using [Snakemake - a scalabale bioinformatics workflow engine](https://doi.org/10.1093/bioinformatics/bts480)


Requirements
------------

This workflow requires the following software to run:

  * [Snakemake][snakemake] 0.10+
  * [Conda][code] (normally comes with Node.js)

Usage
-----


Ripple is easiest to use when installed with [npm][npm]:


Clone workflow into working directory:

```sh
git clone https://github.com/jma1991/rnaseq.git
```

Execute workflow and deploy software dependencies via conda:

```sh
snakemake --use-conda
```

Configuration
-------------

Configure the workflow by editing the files in the `config` folder:

- The `config.yaml` is the basic configuration file for the workflow. 

- The `samples.csv` file contains sample metadata, with 1 row per sample.

- The `units.csv` file

Contributing
------------

To contribute to the workflow, clone this repository locally and commit your code on a separate branch. Please generate unit tests for your code, and run the linter before opening a pull-request:

```sh
snakemake --generate-unit-tests # generate unit tests
snakemake --lint # run the linter
```

You can find more detail in our [Contributing Guide](#). Participation in this open source project is subject to a [Code of Conduct](#).


Support and Migration
---------------------

Ripple major versions are normally supported for 6 months after their last minor release. This means that patch-level changes will be added and bugs will be fixed. The table below outlines the end-of-support dates for major versions, and the last minor release for that version.

We also maintain a [migration guide](#) to help you migrate.

| :grey_question: | Major Version | Last Minor Release | Support End Date |
| :-------------- | :------------ | :----------------- | :--------------- |
| :heart:         | 3             | N/A                | N/A              |
| :hourglass:     | 2             | 2.1                | 2016-07-04       |
| :no_entry_sign: | 1             | 1.4                | 2015-01-26       |

If you're opening issues related to these, please mention the version that the issue relates to.


Thanks
------

I would like to thank Johannes Köster for developing the Snakemake workflow engine, Istvan Albert for, and 




License
-------

This workflow is licensed under the [MIT](#) license.  
Copyright &copy; 2020, James Ashmore


[shield-snakemake]: https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg
[shield-license]: https://img.shields.io/badge/license-MIT-blue.svg
