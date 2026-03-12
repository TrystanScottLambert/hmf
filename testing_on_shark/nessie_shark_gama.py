"""
Script to run Nessie on a GAMA-like mock version of Shark.
"""

import polars as pl
import numpy as np
from nessie import FlatCosmology
from nessie.helper_funcs import create_density_function
from nessie import RedshiftCatalog
from nessie.optimizer import optimize_nm

cosmo = FlatCosmology(h=0.67, omega_matter=0.3)

# READ in Shark Data and make it gama like
MOCK_DIR = "/Users/00115372/Desktop/masking_mock_cat/"
df_galaxies = pl.read_parquet(f"{MOCK_DIR}galaxies_shark.parquet")
df_groups = pl.read_parquet(f"{MOCK_DIR}groups_shark.parquet")

df_groups = df_groups.filter(pl.col("dec") > 0)
df_galaxies = df_galaxies.filter(pl.col("dec") > 0)
# df_galxies = df_galaxies.filter(pl.col("id_fof") != -1)
df_galaxies = df_galaxies.filter(pl.col("mass_stellar_total") > 1e8)
df_galaxies = df_galaxies.filter(pl.col("mag_r_SDSS") < 19.8)
id_counts = df_galaxies.group_by("id_fof").agg(pl.len().alias("count"))
singleton_ids = id_counts.filter(
    (pl.col("count") == 1) & (pl.col("id_fof") != -1)
).select("id_fof")
df_galaxies = df_galaxies.with_columns(
    pl.when(pl.col("id_fof").is_in(singleton_ids["id_fof"]))
    .then(pl.lit(-1))
    .otherwise(pl.col("id_fof"))
    .alias("id_fof")
)

# Prepare the inputs to Nessie
redshifts = np.array(df_galaxies["redshift_observed"])
ra = np.array(df_galaxies["ra"])
dec = np.array(df_galaxies["dec"])
mock_ids = np.array(df_galaxies["id_fof"])
abs_mags = np.array(df_galaxies["mag_r_SDSS"]) - cosmo.dist_mod(redshifts)
total_counts = len(redshifts)
SURVEY_AREA = 7.362531e-03
running_density = create_density_function(redshifts, total_counts, SURVEY_AREA, cosmo)

# Tune and run Nessie
red_cat = RedshiftCatalog(ra, dec, redshifts, running_density, cosmo)
red_cat.set_completeness()
red_cat.mock_group_ids = mock_ids
b0, r0 = optimize_nm(red_cat, 5)
print(f"b0 = {b0}")
print(f"r0 = {r0}")
red_cat.run_fof(b0, r0)

# Write finished data products

df_groups = pl.DataFrame(
    red_cat.calculate_group_table(abs_mags, np.repeat(50, len(redshifts)))
)
df_groups.write_parquet("nessie_groups.parquet")
