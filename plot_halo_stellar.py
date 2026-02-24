import matplotlib.pyplot as plt
from astropy.io import fits
import polars as pl
import numpy as np
from astropy.table import Table


def read_galaxies(galaxy_file: str) -> pl.DataFrame:
    return pl.read_csv(galaxy_file)
    # logmstar, Z, GroupID


def read_groups(group_file: str, n_members: int) -> pl.DataFrame:
    hdu_groups = fits.open(group_file)
    data = Table(hdu_groups[1].data)
    df_pandas = data.to_pandas()
    df_pandas = df_pandas[df_pandas["Nfof"] >= n_members]
    df_pandas = df_pandas[df_pandas["GroupID"] < 400000]
    return pl.from_pandas(df_pandas)


if __name__ == "__main__":
    GROUP_FILE = "data/G3CFoFGroupv10.fits"
    GALAXY_FILE = "data/GAMAGalsInGroups.csv"
    N_CUT = 3

    df_groups = read_groups(GROUP_FILE, N_CUT)
    df_galaxies = read_galaxies(GALAXY_FILE)
    df_group_galaxies = df_galaxies.join(
        df_groups, on="GroupID", maintain_order="right"
    )

    df_members = (
        df_group_galaxies.group_by("GroupID")
        .agg(
            pl.struct(["logmstar", "Z"])
            .sort_by("logmstar", descending=True)
            .get(N_CUT - 1)
            .alias("nth_member")
        )
        .unnest("nth_member")
    )
    df_members = df_members.sort("GroupID")
    df_groups = df_groups.sort("GroupID")

    group_masses = np.log10(np.array(df_groups["MassAfunc"]))
    galaxy_masses = np.array(df_members["logmstar"])
    galaxy_redshift = np.array(df_members["Z"])

    cut = np.where(galaxy_masses > 0)
    plt.scatter(group_masses[cut], galaxy_masses[cut], s=1, color="k", alpha=0.5)
    plt.xlabel(r"$\log M_{\rm halo}$", size=14)
    plt.ylabel(r"$\log M_{*}$", size=14)
    plt.savefig("HSM_relation.png")
    plt.show()

    plt.scatter(galaxy_redshift[cut], galaxy_masses[cut], s=1, color="k")
    plt.xlabel("Redshift", size=14)
    plt.ylabel(r"$\log M_{*}$", size=14)
    plt.savefig("SMZ_relation.png")
    plt.show()

    plt.scatter(galaxy_redshift[cut], group_masses[cut])
    plt.show(
