# rename_vcf_headers

## Description
This module renames the headers in VCF files to include the correct sample names. It ensures compatibility with downstream tools that require properly formatted VCF headers.

## Inputs
- **meta**: Metadata containing the sample name.
- **vcf**: Input VCF file.

## Outputs
- **renamed_vcf**: VCF file with updated headers.

## Usage
Include this module in your Nextflow workflow and provide the required inputs. The module uses `awk` to modify the VCF headers.

## Authors
- Manas Sehgal

## Requirements
- `awk`
