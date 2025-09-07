# metagenome-pipeline

This repository stores the codes to reproduce the analysis in Fujimoto et al. 2024. `software_version.list.txt` lists the version of software used in the analysis. The analysis was conducted on SHIROKANE OS5.

The pipeline codes are written by Dr. Yasumasa Kimura.

## Usage

The codes are intended to be run on [SHIROKANE Supercomputer](https://gc.hgc.jp/en/) environment. Each script submits the necessary jobs by `qsub` in each step per sample. `inifile` package must be installed by `gem install inifile`. The `dir_home_` variable in each script must be set to the pipeline directory holding the `scripts` directory and a user-defined `profile` file specifying the environmental variables while running the script. Also, the software directory (e.g. `cutadapt_dir`) should be set to the path to the software.

Two configuration files must be made to specify the sample and project information as below.

- `info.raw_data.cfg`

```shell
[run1]
project = (Project title)
dir_project_ = (Absolute path to the project directory)
dir_run_ = (Absolute path to the sequence directory [e.g. directory produced after HiSeq CASAVA software])
samples = (Comma-separated sample names)
```

```shell
[run1]
project = Project1
dir_project_ = /home/user/Project1
dir_run_ = /home/user/Project1/fastq
samples = SAMPLE001,SAMPLE002,SAMPLE003
```

- `setting.cfg`

```shell
[project1]
project = (Project title)
type = (virome OR bacteriome)
dir_project_ = (Absolute path to the project directory)
samples = (Comma-separated sample names)
```

```shell
[project1]
project = Project1
type = bacteriome
dir_project_ = /home/user/Project1
samples = SAMPLE001,SAMPLE002,SAMPLE003
```

The database files must be downloaded or compiled, and the custom path must be specified beforehand for the scripts below.

- `08_bacterial_taxonomic_profile.rb`
	- `db` variable
	- Metaphlan database: http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/
- `09_bacterial_contig_taxonomy.rb`
	- `ppsp_dir` variable
	- PhyloPythiaS+ software: https://github.com/algbioi/ppsp/wiki
- `14_pfam_search.split.rb`
	- `db` variable
	- Pfam HMM files: https://www.ebi.ac.uk/interpro/download/Pfam/
- `14_kegg_search.split.rb`
	- [GHOST-MP](https://www.bi.cs.titech.ac.jp/ghostmp/index.html) database must be compiled beforehand by using KEGG GENES fasta file (e.g. `kegg.#{db_date}/genes/fasta/prokaryotes.pep`), and be specified to the `db` path. Accordingly, the `db_date` parameter in the `h_params` should be changed to reflect the version of KEGG database used.
- `scripts/*`
	- The path specification beginning with `/home/user` must be replaced with the path pointing to the databases (e.g. `ko_genes_ = "/home/user/kegg.20180916/genes/ko/ko_genes.list"` in `assign_kegg.blast.py` must be replaced with the path pointing to KEGG database starting with `kegg.#{db_date}`).
    - Note that the files `ko_genes.list`, `ko_module.list`, and `ko_pathway.list` in the scripts directory is the mapping file from `links.tar.gz` or `genes_ko.list.gz` in the database, with the modification that KO identifiers as the first column and linked identifiers as the second.

After the configuration, each step can be run by:

```bash
ruby 01_trim_qc.HiSeq.rb
```
