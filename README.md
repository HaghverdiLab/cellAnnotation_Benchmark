# Adjustments to the reference dataset design improves cell type label transfer

This GitHub repositiory includes the scripts used in the analysis of referece dataset design for cell typle label transfer methods. The pre-print is available at XXX.

The transfer of cell type labels from prior annotated (reference) to newly collected data is an important task in single-cell data analysis. As the number of publicly available annotated datasets which can be used as a reference,
as well as the number of computational methods for cell type label transfer are constantly growing, rationals to understand and decide which reference design and which method to use for a particular query dataset is needed.
Here, we benchmark a set of five popular cell type annotation methods (Seurat, SingleR, CellID, SingleCellNet, ItClust), study the performance on different celltypes and highlight the importance of the design of the reference data (number of cell samples for each cell type,
inclusion of multiple datasets in one reference, gene set selection, etc.) for more reliable predictions.  We introduce a weighted bootstrapping-based approach to imprive the prediction of less represented cell types in the reference data.
