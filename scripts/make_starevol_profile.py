"""Convert Spiroid's Starevol table (examples/data/star/evolution/savgol_10.csv) into a Posidonius evolution file
input/Starevol/M_10.dat read by EvolutionType::Starevol(1.0):  age[yr]  radius[Rsun]  radius_of_gyration_2[total, (I_rad+I_conv)/(M R^2)]
Spiroid stores the moments of inertia in units of M R^2 (star_csv.rs::convert_units), so rg2 = I_rad + I_conv as tabulated."""
import numpy as np, os
src = os.path.expanduser("~/Documents/GitHub/spiroid/examples/data/star/evolution/savgol_10.csv")
dst = os.path.expanduser("~/Documents/GitHub/posidonius/input/Starevol/M_10.dat")
t = np.genfromtxt(src, delimiter=",", names=True)
rg2 = t["radiative_moment_of_inertia"] + t["convective_moment_of_inertia"]
out = np.column_stack([t["age"], t["radius"], rg2, t["convective_moment_of_inertia"], t["mass"]])
np.savetxt(dst, out, fmt="%.16e", header="Starevol 1 Msun (Spiroid savgol_10.csv). columns: age[yr] radius[Rsun] radius_of_gyration_2_total convective_moment_of_inertia/(M R^2) mass[Msun]", comments="# ")
i = np.searchsorted(t["age"], 5e6)
print("rows", len(out), "age range %.3e..%.3e yr" % (t["age"][0], t["age"][-1]), "| at 5 Myr: R = %.6f Rsun, rg2 = %.6f (I_conv/MR2 = %.6f)" % (np.interp(5e6, t["age"], t["radius"]), np.interp(5e6, t["age"], rg2), np.interp(5e6, t["age"], t["convective_moment_of_inertia"])))
print("monotonic age:", np.all(np.diff(t["age"]) > 0))
