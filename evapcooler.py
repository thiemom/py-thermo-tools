# Evaporative cooler calculations - 2024 Atmospheric Composition

from math import log

class EvaporativeCooler:
    """
    Class for performing evaporative air cooling calculations.
    
    This class provides methods to calculate various thermodynamic properties 
    of moist air during an evaporative cooling process, including humidity ratio, 
    enthalpy, dew point temperature, and other related parameters.
    
    Attributes:
        molar_mass_air (float): Molar mass of air in g/mol
        R_air (float): Specific gas constant for dry air in J/(kg·K)
        cp_air (float): Specific heat capacity of dry air in J/(kg·K)
        latent_heat_evaporation (float): Latent heat of water evaporation in kJ/kg at ~25°C
    """

    def __init__(self):
        """
        Initialize the EvaporativeCooler with standard air composition constants.
        
        Attributes:
        - molar_mass_air: Molar mass of air in g/mol
        - R_air: Specific gas constant for dry air in J/(kg·K)
        - cp_air: Specific heat capacity of dry air in J/(kg·K)
        - latent_heat_evaporation: Latent heat of water evaporation in kJ/kg at ~25°C
        """
        # Constants for 2024 air composition
        self.molar_mass_air = 28.966  # g/mol
        self.R_air = 287.042  # J/(kg·K), specific gas constant for dry air
        self.cp_air = 1006  # J/(kg·K), specific heat capacity of dry air
        self.latent_heat_evaporation = 2442.3  # kJ/kg, at ~25°C (adjustable)

    def saturation_pressure(self, T_K):
        """Calculate the saturation pressure of water vapor using the Magnus formula."""
        T_C = T_K - 273.15  # Convert to Celsius
        A = 17.625
        B = 243.04
        return 6.1094 * 10 ** ((A * T_C) / (T_C + B)) * 100  # Convert to Pa

    def humidity_ratio(self, RH, P, T):
        """Calculate the humidity ratio (W) given RH, pressure, and temperature."""
        P_sat = self.saturation_pressure(T)  # Pa
        P_vapor = RH * P_sat  # Partial pressure of water vapor
        return 0.622 * P_vapor / (P - P_vapor)  # kg water vapor / kg dry air

    def relative_humidity(self, W, P, T):
        """Calculate relative humidity from humidity ratio W, pressure, and temperature."""
        P_sat = self.saturation_pressure(T)  # Pa
        P_vapor = W * P / (0.622 + W)  # Partial pressure of water vapor
        return P_vapor / P_sat  # Fractional RH

    def dew_point_temperature(self, W, P):
        """Calculate the dew point temperature (temperature at RH=1)."""
        P_vapor = W * P / (0.622 + W)  # Partial pressure of water vapor
        A = 17.625
        B = 243.04
        alpha = (P_vapor / 100) / 6.1094  # Convert Pa to hPa for Magnus formula
        T_C = (B * (log(alpha) / log(10))) / (A - (log(alpha) / log(10)))  # Celsius
        return T_C + 273.15  # Convert back to Kelvin

    def enthalpy(self, T_K, W):
        """Calculate enthalpy of moist air."""
        return self.cp_air * T_K + W * (self.latent_heat_evaporation * 1000)

    def temperature(self, h, W):
        """Calculate temperature (T) of moist air from enthalpy."""
        return (h - W * (self.latent_heat_evaporation * 1000)) / self.cp_air # K

    def calculate(self, T_in, P_in, RH_in, T_out):
        """Main calculation function for outlet conditions and cooling energy."""
        warnings = []
        if T_out > T_in:
            T_out = T_in
            warnings.append(f"Outlet temperature adjusted to inlet temperature {T_out:.2f}°C to avoid non-physical results.")

        # Convert temperatures to Kelvin
        T_in_K = T_in + 273.15
        T_out_K = T_out + 273.15

        # Inlet properties
        W_in = self.humidity_ratio(RH_in, P_in, T_in_K)
        H_in = self.enthalpy(T_in_K, W_in)

        # Outlet properties
        W_out = W_in + (self.cp_air * (T_in_K - T_out_K)) / (self.latent_heat_evaporation * 1000)
        H_out = self.enthalpy(T_out_K, W_out)

        # Calculate dew point temperature (temperature at RH=1)
        T_dew_K = self.dew_point_temperature(W_out, P_in)

        # If outlet temperature is below dew point, return dew point temperature and raise a warning
        if T_out_K < T_dew_K:
            T_out_K = T_dew_K
            T_out = T_dew_K - 273.15
            warnings.append(f"Outlet temperature adjusted to dew point temperature {T_out_K - 273.15:.2f}°C to avoid non-physical results.")

        # Calculate Relative Humidity
        RH_out = self.relative_humidity(W_out, P_in, T_out_K)

        # Check if warnings were issued
        if len(warnings) == 0:
            warnings = None

        return {
            "W_in": W_in,
            "W_out": W_out,
            "H_in": H_in,
            "H_out": H_out,
            "T_in": T_in,
            "T_out": T_out,
            "RH_in": RH_in,
            "RH_out": RH_out,
            "warnings": warnings,
        }

if __name__ == "__main__":
    
    # Example data and units
    T_in = 33.0  # °C, inlet air temperature
    P_in = 1.013 * 1e5  # Pa, atmospheric pressure
    RH_in = 0.3  # Inlet relative humidity
    T_out = 28.0  # °C, outlet air temperature

    cooler = EvaporativeCooler()
    results = cooler.calculate(T_in, P_in, RH_in, T_out)

    print("Evaporative Cooler Results:")
    print(f"Inlet Humidity Ratio (W_in): {results['W_in']:.6f} kg water vapor / kg dry air")
    print(f"Outlet Humidity Ratio (W_out): {results['W_out']:.6f} kg water vapor / kg dry air")
    print(f"Inlet Enthalpy (H_in): {results['H_in']:.2f} J/kg dry air")
    print(f"Outlet Enthalpy (H_out): {results['H_out']:.2f} J/kg dry air")
    print(f"Outlet Relative Humidity {results['RH_out']:.2f}")
    print(f"Outlet Temperature {results['T_out']:.2f}")
    if results['warnings']: 
        print('\n'.join(results['warnings']))

