# hubmap.build.reference

Builds HuBMAP references from .h5ad files
## Example Workflow

1. **Add cell type markers of interest to a file (e.g, modify `misc/markers.txt`).**

    FeaturePlots of these markers will be automatically generated for use in annotation.

2. **Run `./prepare_integrated.R`** 

    This will make Seurat objects from .h5ad files, filter cells, integrate, process the integrated object (dimensional reduction, neighbor-finding, clustering), and generate plots for use in annotation. 

    Example command (run in download directory): `./prepare_integrated.R --files $(pwd)/data/ed8a4dbbb1554a5e3227d6dfb2368828.h5ad,$(pwd)/data/7fd04d1aba61c35843dd2eb6a19d2545.h5ad --i ed8a4,7fd04 --ensembl $(pwd)/misc/mart_export.txt --markers $(pwd)/misc/markers.txt --o $(pwd)/output` 

3. **View plots and add annotations to a file (e.g, modify `misc/annotations.txt`)**

    Plots are saved in the designated output directory. The objects themselves are also saved in output directory if you want to take a closer look in R/Seurat.

4. **Run `./prune.R`** 

   This will add annotations to the integrated object, prune 'border' cells, and construct final pruned reference with labels + scores. The pruned object is also saved in the designated output directory. 

   Example command (run in download directory): ` ./prune.R --files $(pwd)/data/ed8a4dbbb1554a5e3227d6dfb2368828.h5ad,$(pwd)/data/7fd04d1aba61c35843dd2eb6a19d2545.h5ad --i ed8a4,7fd04 --annotations $(pwd)/misc/annotations.txt --o $(pwd)/output`
 
 ## Documentation
 1. **`prepare_integrated.R: Prepare integrated HuBMAP object.`**

```
Usage:
  prepare_integrated.R [--files <h5ad-files> --i <ids> --ensembl <path-to-ensembl-file> --markers <path-to-markers-file> --o <output-directory> --verbose <verbosity>]

Options:
  --files <h5ad-files>                Absolute paths to h5ad files to be integrated (unspaced, comma-separated).
  --i <ids>                           ID names of objects in h5ad files (unspaced, comma-separated).
  --ensembl <path-to-ensembl-file>    Absolute path to file with ensemble ID dataframe.
  --markers <path-to-markers-file>    Absolute path to file with cell type markers.
  --o <output-directory>              Output directory.
  --verbose <verbosity>               Verbosity [default: TRUE].
  ```

 2. **`prune.R: Annotate and prune integrated HuBMAP object.`**

```
Usage:
  
  prune.R [--files <h5ad-files> --i <ids> --annotations <path-to-annotation-file> --o <output-directory> --verbose <verbosity> --threshold <threshold> --nsample <nsample>]

Options:
  --files <h5ad-files>                        Absolute paths to h5ad files to be integrated (unspaced, comma-separated).
  --i <ids>                                   ID names of objects in h5ad files (unspaced, comma-separated).
  --annotations <path-to-annotation-file>     Absolute path to file with cell type markers.
  --o <output-directory>                      Output directory.
  --threshold <threshold>                     Threshold for pruning [default: 0.75].
  --nsample <nsample>                         Number of subsamples for pre- and post-pruned objects [default: 3].
  --verbose <verbosity>                       Verbosity [default: TRUE].
  ```
