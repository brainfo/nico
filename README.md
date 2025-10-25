# NiCo - Neighborhood informed Cell-type cOmposition analysis

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**NiCo** is a Python package for spatial transcriptomics analysis that performs neighborhood-informed cell-type composition analysis and niche interaction prediction.

## Features

- **Cell Type Annotation**: Transfer labels from scRNA-seq reference data to spatial transcriptomics data using mutual nearest neighbors (MNN) in PCA space
- **Niche Analysis**: Reconstruct spatial niche interaction patterns using logistic regression
- **Gene Covariation**: Analyze gene covariation patterns and ligand-receptor interactions between cell types
- **Visualization**: Comprehensive visualization tools for spatial data, UMAP projections, and interaction networks

## Installation

### Install from source (editable mode)

```bash
# Clone or navigate to the repository
cd /path/to/nico

# Install in editable mode
pip install -e .
```

### Install with development dependencies

```bash
pip install -e ".[dev]"
```

## Quick Start

### Import Styles

NiCo supports multiple import styles to suit your preference:

```python
# Style 1: Import modules with aliases
from nico import Annotations as sann
from nico import Interactions as sint
from nico import Covariations as scov

# Style 2: Import specific functions directly
from nico import find_anchor_cells_between_ref_and_query
from nico import spatial_neighborhood_analysis

# Style 3: Import the whole package
import nico
```

### Example Workflow

```python
import nico
import scanpy as sc

# Load your data
ref_adata = sc.read_h5ad('reference_scrnaseq.h5ad')
query_adata = sc.read_h5ad('spatial_data.h5ad')

# Find anchor cells between reference and query
anchors = nico.find_anchor_cells_between_ref_and_query(
    ref_adata=ref_adata,
    query_adata=query_adata,
    neigh=50,
    no_of_pc=50
)

# Perform NiCo-based annotation
annotations = nico.nico_based_annotation(
    previous=anchors,
    ref_cluster_tag='celltype',
    across_spatial_clusters_dispersion_cutoff=0.15
)

# Save annotations
nico.save_annotations_in_spatial_object(annotations)

# Visualize results
nico.visualize_umap_and_cell_coordinates_with_all_celltypes()

# Perform spatial neighborhood analysis
niche_results = nico.spatial_neighborhood_analysis(
    spatial_cluster_tag='nico_ct',
    Radius=0
)

# Plot niche interactions
nico.plot_niche_interactions_with_edge_weight(niche_results)
```

## Main Functions

### Annotation Module
- `find_anchor_cells_between_ref_and_query()` - Find anchor cells between reference and query datasets
- `nico_based_annotation()` - Perform cell type annotation using NiCo method
- `visualize_spatial_anchored_cell_mapped_to_scRNAseq()` - Visualize anchor mappings
- `save_annotations_in_spatial_object()` - Save annotations to AnnData object
- `visualize_umap_and_cell_coordinates_with_all_celltypes()` - Visualize all cell types
- `visualize_umap_and_cell_coordinates_with_selected_celltypes()` - Visualize selected cell types

### Interactions Module
- `spatial_neighborhood_analysis()` - Perform spatial neighborhood niche analysis
- `find_interacting_cell_types()` - Display regression coefficients for interacting cell types
- `plot_confusion_matrix()` - Plot confusion matrix from niche prediction
- `plot_coefficient_matrix()` - Plot coefficient matrix from niche prediction
- `plot_roc_results()` - Plot ROC curves for cell type predictions
- `plot_niche_interactions_with_edge_weight()` - Plot niche interaction network with edge weights
- `plot_niche_interactions_without_edge_weight()` - Plot niche interaction network without edge weights
- `plot_evaluation_scores()` - Plot evaluation metrics

### Covariations Module
- `gene_covariation_analysis()` - Analyze gene covariation patterns
- `find_LR_interactions_in_interacting_cell_types()` - Find ligand-receptor interactions
- `plot_significant_regression_covariations_as_circleplot()` - Plot covariations as circle plot
- `plot_significant_regression_covariations_as_heatmap()` - Plot covariations as heatmap
- `visualize_factors_in_scRNAseq_umap()` - Visualize factors in scRNA-seq UMAP
- `visualize_factors_in_spatial_umap()` - Visualize factors in spatial UMAP

## Package Structure

```
nico/
├── pyproject.toml          # Package configuration and dependencies
├── README.md               # This file
├── LICENSE                 # MIT License
├── .gitignore             # Git ignore file
└── src/
    └── nico/              # Main package
        ├── __init__.py    # Package initialization with all exports
        ├── Annotations.py # Cell type annotation functions
        ├── Interactions.py # Niche interaction analysis
        ├── Covariations.py # Gene covariation analysis
        └── utils/         # Utility functions
            ├── __init__.py
            ├── SCTransform.py
            └── pyliger_utilities.py
```

## Requirements

- Python ≥ 3.8
- scanpy ≥ 1.9.0
- anndata ≥ 0.8.0
- numpy ≥ 1.20.0
- pandas ≥ 1.3.0
- scipy ≥ 1.7.0
- matplotlib ≥ 3.4.0
- seaborn ≥ 0.11.0
- networkx ≥ 2.6.0
- scikit-learn ≥ 1.0.0
- statsmodels ≥ 0.13.0
- KDEpy ≥ 1.1.0

All dependencies are automatically installed when you install the package.

## Development Tips

### Reloading Modules During Development

When working in editable mode (`pip install -e .`), you can reload modules after making changes:

```python
import importlib
from nico import Annotations as sann

# After editing Annotations.py
importlib.reload(nico.Annotations)
from nico import Annotations as sann  # Re-import to update reference
```

For more details, see [RELOADING.md](RELOADING.md).

### Testing Installation

After installing, you can test all import styles work correctly:
```bash
python test_import.py
```

## Documentation

For detailed documentation and tutorials, please visit the [NiCo documentation](https://github.com/ankitbioinfo/nico).

## Citation

If you use NiCo in your research, please cite:

```
[Citation information will be added]
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

For questions and issues, please open an issue on the [GitHub repository](https://github.com/ankitbioinfo/nico/issues).
