# pk_semaglutide

## Overview

This project provides an implementation of a two-compartment pharmacokinetic (PK) model for semaglutide, a medication used in the treatment of Type 2 diabetes. The model simulates the concentration-time profile of semaglutide in a patient's body following multiple injections and produces a plot of the drug concentration over time.

The PK model is based on the research presented in the article:
> Overgaard, R. V., Delff, P. H., Petri, K. C. C., Anderson, T. W., Flint, A., & Ingwersen, S. H. (2019). Population Pharmacokinetics of Semaglutide for Type 2 Diabetes. Diabetes therapy : research, treatment and education of diabetes and related disorders, 10(2), 649–662. https://doi.org/10.1007/s13300-019-0581-y

## Model Description

The two-compartment model describes the disposition and absorption of semaglutide following subcutaneous or oral administration. It takes into account patient-specific characteristics such as weight, sex, race, age, and injection site. The model is also capable of simulating recurring injection events, allowing users to specify the number of days between each injection, the number of recurring injection events, and the injection concentration.

Covariate effects, as well as reference parameter values, are defined based on the findings of the research article. The model calculates the drug concentration in the central and peripheral compartments over time, using ordinary differential equations.

## Caveats

- The model may not have accurately translated the equations from the original paper over to code; this is a hobby project, and while I have done my best accurately translate the original model over to this code, there is a possibility that I made mistakes.  Please reach out to me with any suggested modificaitons and if I can verify that I made a mistake, I will promply correct.
- The model was developed using clinical data from specific trials and may not be generalizable to all patient populations or dosing conditions.
- The model includes assumptions and parameter estimates from the referenced article. Users should consult the original article for full details on the model development process and limitations.
- The model allows for flexibility in specifying injection schedules, but the accuracy of predictions may vary based on the input data and assumptions.

## Usage

The project contains the following main files:
- `two_compartment_model.py`: Contains the implementation of the `TwoCompartmentModel` class, including methods for calculating patient parameters, simulating drug concentration-time profiles, and plotting results.
- `app.py`: Contains the Streamlit app code for interactive simulation and visualization of the drug concentration-time profile based on user inputs.

To run the Streamlit app locally:
1. Install the required Python packages: `numpy`, `scipy`, `matplotlib`, and `streamlit`.
2. Run the following command in the terminal: `streamlit run app.py`.
3. Open the provided URL in your web browser to access the app and enter the required inputs.

## Citation
Overgaard, R. V., Delff, P. H., Petri, K. C. C., Anderson, T. W., Flint, A., & Ingwersen, S. H. (2019). Population Pharmacokinetics of Semaglutide for Type 2 Diabetes. Diabetes therapy : research, treatment and education of diabetes and related disorders, 10(2), 649–662. https://doi.org/10.1007/s13300-019-0581-y

Please cite the following article when using this project:

## License

Please refer to the LICENSE file for licensing information.