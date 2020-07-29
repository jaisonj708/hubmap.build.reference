# hubmap.build.reference

Builds HuBMAP references from .h5ad files
## Dummy Workflow

1. **Add cell type markers of interest to a file (e.g, modify `misc/markers.txt`).**

    FeaturePlots of these markers will be automatically generated for use in annotation.

2. **Run `./prepare_integrated.R`** 

    This will make Seurat objects from .h5ad files, filter cells, integrate, process the integrated object (dimensional reduction, neighbor-finding, clustering), and generate plots for use in annotation. 

    Example command (run in download directory): `./prepare_integrated.R --files $(pwd)/data/ed8a4dbbb1554a5e3227d6dfb2368828.h5ad,$(pwd)/data/7fd04d1aba61c35843dd2eb6a19d2545.h5ad --i ed8a4,7fd04 --ensembl $(pwd)/misc/mart_export.txt --markers $(pwd)/misc/markers.txt --o $(pwd)/output` 

3. **View plots and add annotations to a file (e.g, modify `misc/annotations.txt`)**

    Plots are saved in the designated output directory. The objects themselves are also saved in output directory if you want to take a closer look in R/Seurat.

     Example command (run in download directory): # ./prune.R --files `$(pwd)/data/ed8a4dbbb1554a5e3227d6dfb2368828.h5ad,$(pwd)/data/7fd04d1aba61c35843dd2eb6a19d2545.h5ad --i ed8a4,7fd04 --annotations $(pwd)/misc/annotations.txt --o $(pwd)/output`

4. **Run `./prune.R`** 

   This will add annotations to the integrated object, prune 'border' cells, and construct final pruned reference with labels + scores. The pruned object is also saved in the designated output directory. 
 
