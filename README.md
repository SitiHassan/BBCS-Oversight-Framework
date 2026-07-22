# NHS Oversight Framework Metric Engine

This git repository contains a reusable metric engine for processing [NHS Oversight Framework](https://www.england.nhs.uk/nhs-oversight-framework/) metrics.

The engine was developed using the existng [Outcomes Framework](https://github.com/BBCS-PHI/2_BSOL_Outcomes_Framework) metric engine as a foundation and adapted to meet the Oversight Framework's specific data structures, calculation methods and technical requirements. 

The engine is metadata-driven, using a central control file to define the downstream processing and calculation logic for each metric.

## Purpose
The project provides a consistent process for:
* importing metric data
* applying metric-specific calculations
* producing standardised outputs
* carrying out data quality checks
* preparing results for reporting and downstream analysis

# Supported calculations
The engine can process a range of metric types, including:
* counts
* percentages
* proportions
* directly age-standardised rates
* crude rates
* ratios
* percentage changes
* percentage point differences
  
# Running the project
1. Clone the repository
2. Open the R project in RStudio
3. Update the required input paths and configuration
4. Run the main processing script
5. Review the generated outputs and DQ checks
   
This repository is dual licensed under the [Open Government v3]([https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/) & MIT. All code and outputs are subject to Crown Copyright.
