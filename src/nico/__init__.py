# Import modules so they can be accessed as: from nico import Annotations as sann
# Using relative imports to ensure proper module registration for reload
from . import Annotations
from . import Interactions
from . import Covariations

# Import all functions for direct access: from nico import function_name
from .Annotations import create_directory
from .Annotations import find_index
from .Annotations import find_match_index_in_dist
from .Annotations import find_mutual_nn
from .Annotations import sct_return_sc_sp_in_shared_common_PC_space
from .Annotations import find_annotation_index
from .Annotations import find_commnon_MNN
from .Annotations import visualize_spatial_anchored_cell_mapped_to_scRNAseq
from .Annotations import find_anchor_cells_between_ref_and_query
from .Annotations import nico_based_annotation
from .Annotations import read_dist_and_nodes_as_graph
from .Annotations import return_singlecells
from .Annotations import findSpatialCells
from .Annotations import find_all_the_spatial_cells_mapped_to_single_cells
from .Annotations import write_annotation
from .Annotations import find_unmapped_cells_and_deg
from .Annotations import resolved_confused_and_unmapped_mapping_of_cells_with_weighted_average_of_inverse_distance_in_neighbors
from .Annotations import resolved_confused_and_unmapped_mapping_of_cells_with_majority_vote
from .Annotations import visualize_umap_and_cell_coordinates_with_all_celltypes
from .Annotations import visualize_umap_and_cell_coordinates_with_selected_celltypes
from .Annotations import remove_extra_character_from_name
from .Annotations import plot_all_ct
from .Annotations import plot_specific_ct
from .Annotations import delete_files
from .Covariations import plot_significant_regression_covariations_as_circleplot
from .Covariations import plot_significant_regression_covariations_as_heatmap
from .Covariations import save_LR_interactions_in_excelsheet_and_regression_summary_in_textfile_for_interacting_cell_types
from .Covariations import extract_and_plot_top_genes_from_chosen_factor_in_celltype
from .Covariations import plot_top_genes_for_pair_of_celltypes_from_two_chosen_factors
from .Covariations import plot_top_genes_for_a_given_celltype_from_all_factors
from .Covariations import find_interest_of_genes
from .Annotations import save_annotations_in_spatial_object
from .Covariations import gene_covariation_analysis
from .Covariations import plot_cosine_and_spearman_correlation_to_factors
from .Covariations import plot_feature_matrices
from .Covariations import find_LR_interactions_in_interacting_cell_types
from .Covariations import make_excel_sheet_for_gene_correlation
from .Covariations import pathway_analysis
from .Covariations import create_directory
from .Covariations import find_index
from .Covariations import read_spatial_data
from .Covariations import find_correlation_bw_genes_and_PC_component_in_singlecell
from .Covariations import find_correlation_bw_genes_and_PC_component_in_singlecell_cosine
from .Covariations import top_genes_in_correlation_list_without
from .Covariations import alignment_score
from .Covariations import multiplicative_method
from .Covariations import remove_extra_character_from_name
from .Covariations import find_PC_of_invidualCluster_in_SC
from .Covariations import makePCneighboorhoodFeatureMatrix
from .Covariations import compute_PC_space
from .Covariations import model_linear_regression
from .Covariations import run_ridge_regression
from .Covariations import find_logistic_regression_interacting_score
from .Covariations import findXYZC
from .Covariations import create_subtitle
from .Covariations import find_fold_change
from .Covariations import sorting_of_factors_for_showing_the_value_in_excelsheet
from .Covariations import triangulation_for_triheatmap
from .Covariations import  plot_ligand_receptor_in_interacting_celltypes
from .Covariations import visualize_factors_in_scRNAseq_umap
from .Covariations import plot_all_ct
from .Covariations import visualize_factors_in_spatial_umap
from .Covariations import read_LigRecDb
from .Covariations import sort_index_in_right_order
from .Interactions import create_directory
from .Interactions import findNeighbors_in_given_radius
from .Interactions import find_neighbors
from .Interactions import create_spatial_CT_feature_matrix
from .Interactions import euclidean_dist
from .Interactions import reading_data
from .Interactions import plot_multiclass_roc
from .Interactions import plot_roc_results
from .Interactions import plot_confusion_matrix
from .Interactions import plot_coefficient_matrix
from .Interactions import plot_predicted_probabilities
from .Interactions import plot_niche_interactions_with_edge_weight
from .Interactions import plot_niche_interactions_without_edge_weight
from .Interactions import read_processed_data
from .Interactions import model_log_regression
from .Interactions import find_interacting_cell_types
from .Interactions import remove_extra_character_from_name
from .Interactions import spatial_neighborhood_analysis
from .Interactions import plot_evaluation_scores

__version__ = '1.5.0'

# Export modules and main functions
__all__ = [
    # Modules
    'Annotations',
    'Interactions', 
    'Covariations',
    # Main functions - Annotations
    'find_anchor_cells_between_ref_and_query',
    'nico_based_annotation',
    'visualize_spatial_anchored_cell_mapped_to_scRNAseq',
    'save_annotations_in_spatial_object',
    'visualize_umap_and_cell_coordinates_with_all_celltypes',
    'visualize_umap_and_cell_coordinates_with_selected_celltypes',
    'delete_files',
    # Main functions - Interactions
    'spatial_neighborhood_analysis',
    'find_interacting_cell_types',
    'plot_confusion_matrix',
    'plot_coefficient_matrix',
    'plot_predicted_probabilities',
    'plot_roc_results',
    'plot_niche_interactions_with_edge_weight',
    'plot_niche_interactions_without_edge_weight',
    'plot_evaluation_scores',
    # Main functions - Covariations
    'gene_covariation_analysis',
    'find_LR_interactions_in_interacting_cell_types',
    'plot_significant_regression_covariations_as_circleplot',
    'plot_significant_regression_covariations_as_heatmap',
    'visualize_factors_in_scRNAseq_umap',
    'visualize_factors_in_spatial_umap',
    'plot_ligand_receptor_in_interacting_celltypes',
    'pathway_analysis',
]
