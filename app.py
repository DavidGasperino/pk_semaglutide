# Streamlit app for the Two-Compartment Model
import streamlit as st

# Import the TwoCompartmentModel class
from two_compartment_model import TwoCompartmentModel

# Define patient characteristics input fields
st.title("Two-Compartment Model")
st.header("Patient Characteristics")
weight = st.number_input("Body weight (kg)", min_value=0.0, step=0.1)
sex = st.radio("Sex", options=["female", "male"])
race = st.radio("Race", options=["white", "black", "asian"])
ethnicity = st.radio("Ethnicity", options=["non-Hispanic", "hispanic"])
glycaemic_status = st.radio("Glycaemic status", options=["normoglycemia", "T2D"])
age_group = st.radio("Age group", options=["<= 65", "> 65"])
hepatic_impairment = st.checkbox("Hepatic impairment")
injection_site = st.radio("Injection site", options=["abdomen", "thigh"])
drug_product_strength = st.number_input(
    "Drug product strength (mg)", min_value=0.25, max_value=10.0, step=0.01
)
T2D = st.checkbox("Type 2 diabetes")

# Get patient parameters
patient_params = TwoCompartmentModel.get_patient_params(
    weight=weight,
    sex=sex,
    race=race,
    ethnicity=ethnicity,
    glycaemic_status=glycaemic_status,
    age_group=age_group,
    hepatic_impairment=hepatic_impairment,
    injection_site=injection_site,
    drug_product_strength=drug_product_strength,
    T2D=T2D,
)

# Define recurring injection event input fields
st.header("Recurring Injection Events")
days_between_injections = st.number_input(
    "Days between each injection", min_value=1, step=1
)
num_injections = st.number_input(
    "Number of recurring injection events", min_value=1, step=1
)
injection_concentration = st.number_input(
    "Injection concentration (mg)", min_value=0.25, max_value=10.0, step=0.01
)
days_before_first_injection = st.number_input(
    "Days before first injection", min_value=0, step=1
)

# Create an instance of the TwoCompartmentModel class
model = TwoCompartmentModel(patient_params, {})

# Create recurring injection schedule
injections = model.create_recurring_injection_schedule(
    days_between_injections,
    num_injections,
    injection_concentration,
    days_before_first_injection,
)

# Define end time input field
st.header("Simulation")
t_end = st.number_input("End time for the simulation (days)", min_value=1, step=1)

# Simulate the drug concentration-time profile
time, drug_concentration = model.simulate_concentration_time_profile(
    t_end * 24
)  # Convert t_end to hours

# Plot the time series
model.plot_time_series(time, drug_concentration)
