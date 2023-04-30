import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint


class TwoCompartmentModel:
    def __init__(self, patient_params, injections, ka_values):
        self.patient_params = patient_params
        self.injections = self.convert_to_hours(injections)
        self.ka_values = ka_values  # Store the list of ka values

    @staticmethod
    def get_patient_params(
        weight,
        sex,
        race,
        ethnicity,
        glycaemic_status,
        age_group,
        hepatic_impairment,
        injection_site,
        weight_ref=85,
    ):
        # Covariate effects
        theta_weight_CL = 1.01
        theta_weight_V = 0.923
        theta_sex_male = 1  # No effect found in the study
        theta_black_race = 1  # No effect found in the study
        theta_asian_race = 1  # No effect found in the study
        theta_hispanic_ethnicity = 1  # No effect found in the study
        theta_T2D_glycaemic_status = 1.12
        theta_age_group_over_65 = 1  # No effect found in the study
        theta_hepatic_impairment = 1  # No effect found in the study
        theta_thigh_F = 0.883

        # Reference parameter values
        CL_ref = 0.0348  # L/h
        Q_ref = 0.304  # L/h
        Vc_ref = 3.59  # L
        Vp_ref = 4.10  # L
        ka_ref = 0.0253  # h^-1
        F_ref = 0.847  # unitless

        # Covariate adjustment and parameter calculation
        E_weight_CL = (weight / weight_ref) ** theta_weight_CL
        E_weight_V = (weight / weight_ref) ** theta_weight_V
        E_sex = theta_sex_male if sex == "male" else 1
        E_race = (
            theta_black_race
            if race == "black"
            else theta_asian_race
            if race == "asian"
            else 1
        )
        E_ethnicity = theta_hispanic_ethnicity if ethnicity == "hispanic" else 1
        E_glycaemic_status = (
            theta_T2D_glycaemic_status if glycaemic_status == "T2D" else 1
        )
        E_age_group = theta_age_group_over_65 if age_group == ">65" else 1
        E_hepatic_impairment = theta_hepatic_impairment if hepatic_impairment else 1
        E_injection_site_F = theta_thigh_F if injection_site == "thigh" else 1

        # Final parameter values for the individual
        CL = (
            CL_ref
            * E_weight_CL
            * E_sex
            * E_race
            * E_ethnicity
            * E_glycaemic_status
            * E_age_group
            * E_hepatic_impairment
        )
        Q = Q_ref * E_weight_CL  # Q is scaled with the same factor as CL
        Vc = Vc_ref * E_weight_V
        Vp = Vp_ref * E_weight_V  # Vp is scaled with the same factor as Vc
        F = F_ref * E_injection_site_F

        return {"CL": CL, "Q": Q, "Vc": Vc, "Vp": Vp, "F": F}

    @staticmethod
    def interpolate_ka(dosage):
        # Define the ka values for specific dosages
        dosages_ka = {
            1: 0.0346,
            3: 0.0526,
            10: 0.139,
        }

        # Ensure dosage is within the range of 0.25 mg to 10 mg
        if dosage < 0.25 or dosage > 10:
            raise ValueError("Dosage must be between 0.25 mg and 10 mg (inclusive).")

        # If the dosage has a corresponding ka value, return it
        if dosage in dosages_ka:
            return dosages_ka[dosage]

        # Otherwise, interpolate ka value based on the closest dosages
        sorted_dosages = sorted(dosages_ka.keys())
        lower_dosages = list(filter(lambda d: d < dosage, sorted_dosages))
        lower_dosage = max(lower_dosages) if lower_dosages else min(sorted_dosages)
        upper_dosages = list(filter(lambda d: d > dosage, sorted_dosages))
        upper_dosage = min(upper_dosages) if upper_dosages else max(sorted_dosages)

        lower_ka = dosages_ka[lower_dosage]
        upper_ka = dosages_ka[upper_dosage]

        # Check if upper_dosage and lower_dosage are the same
        if upper_dosage == lower_dosage:
            return lower_ka

        # Perform interpolation
        ka = lower_ka + (upper_ka - lower_ka) * (dosage - lower_dosage) / (
            upper_dosage - lower_dosage
        )
        return ka

    def create_recurring_injection_schedule(
        self,
        days_between_injections,
        num_injections,
        injection_concentration,
        days_before_first_injection,
    ):
        # Initialize an empty injection schedule dictionary
        injection_schedule = {}

        # Calculate the time of each injection event
        for i in range(num_injections):
            injection_time_days = (
                days_before_first_injection + i * days_between_injections
            )
            injection_schedule[injection_time_days] = injection_concentration

        return injection_schedule

    @staticmethod
    def convert_to_hours(injections):
        # Convert the times in the injections dictionary from days to hours
        injections_hours = {
            t_days * 24: dose_mg for t_days, dose_mg in injections.items()
        }
        return injections_hours

    def simulate_concentration_time_profile(self, t_end):
        def pk_model(y, t, CL, Q, Vc, Vp, ka, F):
            A_c, A_p, A_g = y
            dA_cdt = Q * (A_p / Vp - A_c / Vc) - CL * A_c / Vc + ka * A_g
            dA_pdt = Q * (A_c / Vc - A_p / Vp)
            dA_gdt = -ka * A_g
            return [dA_cdt, dA_pdt, dA_gdt]

        # Automatically calculate n_points by adding 200 points for every 7 days
        n_points = int((t_end / 7) * 200)

        # Time array for simulation
        t_eval = np.linspace(0, t_end, n_points)
        y_init = [
            0,
            0,
            0,
        ]  # Initial amounts in central, peripheral, and gut compartments

        # Initialize lists to store time and concentration data
        time_data = []
        concentration_data = []

        # Sort the injection schedule by time
        sorted_injections = sorted(self.injections.items(), key=lambda x: x[0])
        next_dose_index = 0  # Index of the next scheduled dose

        # Set default ka value to the value for the first dosing schedule
        ka = self.ka_values[0] if self.ka_values else 0.0253  # Default ka value

        for i, t in enumerate(t_eval):
            # Check if there is a dose administration at the current time point
            if (
                next_dose_index < len(sorted_injections)
                and t >= sorted_injections[next_dose_index][0]
            ):
                dose_time, dose_mg = sorted_injections[next_dose_index]
                y_init[2] += (
                    dose_mg / self.patient_params["F"]
                )  # Add dose to gut compartment
                ka = self.ka_values[
                    next_dose_index
                ]  # Update ka value for the current dosing schedule
                next_dose_index += 1  # Update the index of the next scheduled dose

            # Solve the ODEs for the PK model
            sol = odeint(
                pk_model,
                y_init,
                [t, t_eval[i + 1] if i + 1 < len(t_eval) else t_end],
                args=(
                    self.patient_params["CL"],
                    self.patient_params["Q"],
                    self.patient_params["Vc"],
                    self.patient_params["Vp"],
                    ka,
                    self.patient_params["F"],
                ),
            )
            y_init = sol[-1]  # Update initial conditions for the next step

            # Append the current time and concentration in the central compartment to the data lists
            time_data.append(t)
            concentration_data.append(y_init[0] / self.patient_params["Vc"])

        return time_data, concentration_data

    def plot_time_series(self, time, drug_concentration):
        # Create a figure and axis using plt.subplots()
        fig, ax = plt.subplots(figsize=(10, 6))

        # Molecular weight of semaglutide (g/mol)
        molecular_weight_semaglutide = 4113.57

        # Convert drug concentration from mg/L to nmol/L
        drug_concentration_nm = (
            np.array(drug_concentration) * 1e6
        ) / molecular_weight_semaglutide

        # Convert time from hours to weeks
        time_weeks = np.array(time) / (24 * 7)

        # Plot the data on the axis
        ax.plot(time_weeks, drug_concentration_nm, label="Semaglutide Concentration")

        # Plot the injection points on the plot line
        for t_hours, dose_mg in self.injections.items():
            # Convert time of injection from hours to weeks
            t_weeks = t_hours / (24 * 7)

            # Find the index of the time point closest to the injection time
            idx = np.argmin(np.abs(time_weeks - t_weeks))

            # Get the corresponding drug concentration value at the injection time
            concentration_at_t = drug_concentration_nm[idx]

            # Plot a marker at the injection point on the plot line
            ax.plot(t_weeks, concentration_at_t, "ko", markersize=5)

            # Label the marker with the dose in mg
            ax.text(
                t_weeks, concentration_at_t, f"{dose_mg} mg", ha="left", va="bottom"
            )

        # Set up x-axis ticks
        max_weeks = np.ceil(max(time_weeks))  # Get the maximum number of weeks
        major_ticks = np.arange(0, max_weeks + 1, 1)  # Major ticks at each week
        minor_ticks = np.arange(
            0, max_weeks + 1, 1 / 7
        )  # Minor ticks with 6 divisions between each week

        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor=True)

        # Customize vertical grid lines
        ax.grid(
            which="major", axis="x", linestyle="-", color="grey"
        )  # Darker lines for major ticks (weeks)
        ax.grid(
            which="minor", axis="x", linestyle="-", color="lightgrey"
        )  # Lighter lines for minor ticks

        # Only show labels for the weeks (integer values)
        ax.set_xticklabels(
            [int(tick) if tick.is_integer() else "" for tick in ax.get_xticks()]
        )

        ax.set_xlabel("Time (weeks)")
        ax.set_ylabel("Semaglutide Concentration (nmol/L)")
        ax.set_title("Two-Compartment Model: Semaglutide Concentration vs. Time")
        ax.legend()
        ax.grid(True)

        # Return the figure object
        return fig
