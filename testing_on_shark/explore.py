import polars as pl
import numpy as np
import pylab as plt

mock_directory = "/Users/00115372/Desktop/masking_mock_cat/"
group_file_name = "groups_shark.parquet"
galaxy_file_name = "galaxies_shark.parquet"

df_groups = pl.read_parquet(f"{mock_directory}{group_file_name}")
# df_groups = df_groups.filter(pl.col("mass_virial") > 1e12)
df_groups = df_groups.filter(pl.col("dec") > 0)
df_galaxies = pl.read_parquet(f"{mock_directory}{galaxy_file_name}")
df_galaxies = df_galaxies.filter(pl.col("dec") > 0)
df_galxies = df_galaxies.filter(pl.col("id_fof") != -1)
df_galaxies = df_galaxies.filter(pl.col("mass_stellar_total") > 1e8)
df_galaxies = df_galaxies.filter(pl.col("mag_r_SDSS") < 19.8)

df_brightest = df_galaxies.group_by("id_group_sky").agg(
    pl.col("mass_stellar_total").max()
)
df_plot = df_brightest.join(df_groups, on="id_group_sky")

plt.scatter(
    df_plot["redshift_cosmological"], np.log10(df_plot["mass_stellar_total"]), s=0.1
)
plt.xlim(0, 0.4)
plt.ylim(4, 11)
plt.show()
