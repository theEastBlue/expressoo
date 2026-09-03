# expressoo

## Features

### 1. Metadata Proportions & Dynamic Filtering
Filtering based on specific cellular hierarchies (L1, L2, L3) among other metrics. The interface includes stacked barplots to visualize population distributions, such as cell type proportions grouped by Sex.

![Metadata Proportions & Filtering](/docs/assets/cell-type.png)

### 2. QC Metrics (Boxplots)
Evaluate quality control metrics across specific cell populations. Users can see transcript counts (`nCount_RNA`) for selected cell types, like CD8 T cells, split by metadata variables.

![QC Metrics Boxplot](/docs/assets/box-plot.png)

### 3. Module Score & Gene Expression (Violin Plots)
Analyse specific expression signatures across your cohorts (e.g., `ifnk_rna1`). Users can toggle between raw counts and `Log(1+x)` normalisation on the fly to visualise expression distributions across selected cell types.

![Gene Expression Violin Plot](/docs/assets/violin.png)

### 4. High-Performance WebGL UMAPs
The Pool UMAP tab uses Plotly WebGL to render thousands of cells smoothly (showing 17,428 sampled points here). It toggles dataset modes and adjust marker size and opacity to handle overplotting.

![Pool UMAP View](/docs/assets/pool-umap-cell-types.png)

### 5. Mapping Score
Evaluate reference mapping quality with a UMAP controlled by a mapping score range slider. This allows users to dynamically filter cells based on their confidence scores.

![Mapping Score](/docs/assets/pool-umap-mapping-score.png)