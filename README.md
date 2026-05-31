# nf-ucsc-hub

> **⚠️ Work in progress** — this pipeline is under active development.
 
Nextflow pipeline for building UCSC track hubs from genome assemblies.
 
## Requirements
 
UCSC tools on `$PATH`: `faToTwoBit`, `twoBitInfo`, `gtfToGenePred`, `genePredToBed`, `bedToBigBed`, `ixIxx`
 
## Usage
 
```bash
nextflow run main.nf \
    --genome genome.fa \
    --gtf genes.gtf \
    --trf trf_output.dat \
    --outdir ./hub
```
 
## Parameters
 
| Parameter    | Description                        | Default                     |
|-------------|------------------------------------|-----------------------------|
| `--genome`  | Genome FASTA                       | required                    |
| `--gtf`     | Gene annotation GTF                | required                    |
| `--trf`     | TRF `.dat` output                  | required                    |
| `--assembly`| Assembly name for the hub          | derived from genome filename|
| `--outdir`  | Output hub directory               | `./hub`                     |
| `--email`   | Contact email shown in hub         | `your@email.cz`             |
 
## Output
 
```
hub/
├── hub.txt
├── genomes.txt
└── <assembly>/
    ├── <assembly>.2bit
    ├── trackDb.txt
    ├── genes.bb
    ├── genes.ix
    ├── genes.ixx
    └── trf.bb
```
 
## Tracks
 
- **Genes** — BED12 BigBed from GTF, with searchable gene name index
- **TRF** — BED3 BigBed from Tandem Repeat Finder `.dat` output