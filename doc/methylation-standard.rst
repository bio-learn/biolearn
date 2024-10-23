===============================
DNA Methylation Array Data Standard V-2410
===============================

This standard describes a consistent file format for methylation matrix data and associated metadata. 
It focuses on file formatting for storage and transmission without mandating any specific form of data processing.

Why is this standard needed?
============================

GEO Series Matrix files are a very commonly used standard for representing methylation data and sample metadata. However, they have two problems:

- While metadata field inclusion follows a standard format, GEO does not standardize metadata values even for common data fields like sex and age.
- As methylation data has grown larger, the matrix can no longer be included in the Series Matrix files and is instead added as a supplementary file; however, these files do not have a standard format.

Preprocessing
=============

This standard assumes that the array data has been pre-processed, including quality control, batch correction, and normalization. 
We recommend using SeSAMe for preprocessing, as recommended by Illumina.

Methylation Matrix File Standard
================================

File Format
----------

- Files must be formatted as CSV (Comma-Separated Values).
- Each row in the file must contain the same number of entries.
- **Rows**: Each row corresponds to measurements for a CpG site, with the CpG ID in the first column.
- **Columns**: Each column corresponds to a sample, with the sample ID as the header.
    - **Header Row**: The first row must contain sample IDs.
    - If published on GEO, the sample IDs must match GEO IDs exactly.
- **Length**: No more than 1000 samples should be included in a single file.

Data Values
----------

- **Beta Values**: All data values must be beta values (methylation levels).
- **Range**: Beta values must be decimal numbers between 0 and 1 or “NaN”.
- **Precision**: Beta values must be rounded to no more than three decimal places.
- **Handling Multiple Values**: If multiple values exist for the same CpG site, use the average of the available values.

Metadata File Standard
======================

File Format
----------

- Files must be formatted as CSV (Comma-Separated Values).
- Each row represents a sample, and each column represents a metadata field.

Standard Metadata Fields
------------------------

- **Sex**: Integer value with 0 for female and 1 for male. “NaN” (no values) is permitted for unknown values.
- **Age**: Integer or decimal value representing age in years.

Additional Metadata Fields
--------------------------

- Other metadata fields are permitted.

Example Files
=============

Example DNAm File::

    cpgSite     GSM3074480    GSM3074481    GSM3074482    GSM3074483
    cg00000029    0.395         0.499         0.42         0.423
    cg00000109    0.821         0.868         0.877        0.841
    cg00000155    0.878         0.89          0.894        0.822

Example Metadata File::

    SampleID    Sex    Age     Disease State
    GSM3074480    0    25.2    None
    GSM3074481    0    29.5    COVID
    GSM3074482    1    39.7    None

References
==========

- `Data standards – ENCODE <https://www.encodeproject.org/data-standards/>`_
- `SOFT file format and content - GEO - NCBI <https://www.ncbi.nlm.nih.gov/geo/info/soft.html>`_
- `Methylation Array Data Analysis Tips - Illumina <https://www.illumina.com/techniques/microarrays/methylation-arrays/methylation-array-data-analysis-tips.html>`_
