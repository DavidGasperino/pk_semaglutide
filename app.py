import streamlit as st
from two_compartment_model import TwoCompartmentModel

# Add a title across the top of the app
st.title("Semaglutide Pharmacokinetic Model")

# Input fields for patient parameters in the left sidebar
st.sidebar.header("Patient Parameters")
weight = st.sidebar.number_input("Weight (kg)", min_value=0.0, value=85.0)
sex = st.sidebar.selectbox("Sex", options=["male", "female"])
race = st.sidebar.selectbox("Race", options=["white", "black", "asian", "other"])
ethnicity = st.sidebar.selectbox("Ethnicity", options=["non-hispanic", "hispanic"])
glycaemic_status = st.sidebar.selectbox("Glycaemic Status", options=["T2D", "non-T2D"])
age_group = st.sidebar.selectbox("Age Group", options=["<=65", ">65"])
hepatic_impairment = st.sidebar.checkbox("Hepatic Impairment")
injection_site = st.sidebar.selectbox("Injection Site", options=["abdomen", "thigh"])

# Create two columns (dosing schedule column, plot column) with relative widths
middle_col, right_col = st.columns(2)

# Dosing schedule input and management (in the middle column)
middle_col.header("Dosing Schedule")
dosing_schedules = st.session_state.get("dosing_schedules", [])

if "add_schedule" not in st.session_state:
    st.session_state.add_schedule = False

if st.button("Add Dosing Schedule"):
    st.session_state.add_schedule = True

# Calculate the default value for the "Start Day" input
default_start_day = 0
if dosing_schedules:
    last_schedule = dosing_schedules[-1]
    last_start_day = last_schedule["start_day"]
    last_interval_days = last_schedule["interval_days"]
    last_recurring_events = last_schedule["recurring_events"]
    last_end_day = last_start_day + (last_interval_days * (last_recurring_events - 1))
    default_start_day = last_end_day  # Start one day after the last schedule ends

# Add a "Clear" button to clear all dosing schedules
if st.button("Clear Dosing Schedules"):
    dosing_schedules.clear()
    st.session_state.dosing_schedules = dosing_schedules

if st.session_state.add_schedule:
    dose_mg = st.number_input(
        "Dosage Concentration (mg)",
        min_value=0.25,  # Minimum allowed value
        max_value=10.0,  # Maximum allowed value
        value=0.25,  # Initial value
        step=0.01,  # Increment step
        format="%.2f",  # Number format
    )
    start_day = st.number_input(
        "Start Day",
        min_value=0,
        value=default_start_day,  # Use the calculated default value
        step=1,
        format="%i",
    )
    interval_days = st.number_input(
        "Interval Between Injections (days)", min_value=0, value=7, step=1, format="%i"
    )
    recurring_events = st.number_input(
        "Number of Recurring Injection Events",
        min_value=0,
        value=4,
        step=1,
        format="%i",
    )
    if st.button("Confirm Dosing Schedule"):
        dosing_schedules.append(
            {
                "dose_mg": dose_mg,
                "start_day": start_day,
                "interval_days": interval_days,
                "recurring_events": recurring_events,
            }
        )
        st.session_state.dosing_schedules = dosing_schedules
        st.session_state.add_schedule = False

# Display existing dosing schedules (in the middle column)
for idx, schedule in enumerate(dosing_schedules):
    middle_col.write(
        f"Dosing Schedule {idx + 1}: {schedule['dose_mg']} mg every {schedule['interval_days']} days, starting on day {schedule['start_day']}, with {schedule['recurring_events']} recurring events."
    )

# Convert dosing schedules into the injections dictionary format required by the TwoCompartmentModel
# and calculate ka values for each dosing schedule
injections = {}
ka_values = {}
for schedule in dosing_schedules:
    dose_mg = schedule["dose_mg"]
    ka_value = TwoCompartmentModel.interpolate_ka(dose_mg)
    for i in range(schedule["recurring_events"]):
        injections[schedule["start_day"] + i * schedule["interval_days"]] = dose_mg
        ka_values[schedule["start_day"] + i * schedule["interval_days"]] = ka_value

# Simulate and plot the PK profile if dosing schedules have been added
if dosing_schedules:
    # Get patient parameters (excluding "ka" attribute as it's not needed)
    patient_params = TwoCompartmentModel.get_patient_params(
        weight=weight,
        sex=sex,
        race=race,
        ethnicity=ethnicity,
        glycaemic_status=glycaemic_status,
        age_group=age_group,
        hepatic_impairment=hepatic_impairment,
        injection_site=injection_site,
    )

    # Create a model instance with the patient parameters and empty injections
    model = TwoCompartmentModel(patient_params, {}, [])

    # Convert dosing schedules into the injections dictionary format required by the TwoCompartmentModel
    # and calculate ka values for each dosing schedule
    injections = {}
    ka_values = []
    for schedule in dosing_schedules:
        dose_mg = schedule["dose_mg"]
        interval_days = schedule["interval_days"]
        start_day = schedule["start_day"]
        recurring_events = schedule["recurring_events"]
        ka_value = TwoCompartmentModel.interpolate_ka(dose_mg)
        schedule_injections = model.create_recurring_injection_schedule(
            interval_days, recurring_events, dose_mg, start_day
        )
        for t, dose in schedule_injections.items():
            injections[t] = dose
            ka_values.append(ka_value)  # Append the ka value for each injection event

    # Update the model instance with the new injections and ka_values
    model.injections = model.convert_to_hours(injections)
    model.ka_values = ka_values

    # Calculate the end day of the last dosing schedule
    last_schedule = dosing_schedules[-1]
    last_start_day = last_schedule["start_day"]
    last_interval_days = last_schedule["interval_days"]
    last_recurring_events = last_schedule["recurring_events"]
    last_end_day = last_start_day + (last_interval_days * (last_recurring_events - 1))

    # Calculate t_end based on the end day of the last dosing schedule
    # and extend the simulation by 7 days (24*7 hours) after the last dose
    t_end = (last_end_day + 7) * 24
    time_points, concentrations = model.simulate_concentration_time_profile(t_end)

    # Plot the concentration-time profile (in the right column)
    fig = model.plot_time_series(time_points, concentrations)
    right_col.pyplot(fig)
else:
    middle_col.warning("Please add at least one dosing schedule before simulating.")

# Add a hyperlink to the original document at the bottom of the sidebar
st.sidebar.markdown(
    "Link to the article this model is based on: "
    "[Semaglutide Pharmacokinetic Modeling](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6437231/?report=classic)"
)

# Add a hyperlink to the GitHub page for the project.
st.sidebar.markdown(
    "Code can be found here: "
    "[pk_semaglutide](https://github.com/DavidGasperino/pk_semaglutide)"
)

st.markdown("## Disclaimer")
st.markdown(
    "This application is for educational and informational purposes only. "
    "It is not intended to provide medical advice or to replace consultation "
    "with a qualified healthcare provider. The results of the simulation are "
    "based on a mathematical model and should not be used for clinical decision-making."
)
