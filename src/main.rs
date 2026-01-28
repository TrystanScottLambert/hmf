use cosmoxide::Cosmology;
use polars::prelude::*;
use polars::{error::PolarsResult, frame::DataFrame, prelude::ParquetReader};
use std::fs::File;
use std::panic::panic_any;

fn calculate_max_redshift(
    redshift: f64,
    apparent_mag: f64,
    apparent_mag_limit: f64,
    cosmology: &Cosmology,
) -> f64 {
    let absolute_mag = apparent_mag - cosmology.distance_modulus(redshift);
    let max_distance = 10_f64.powf(((apparent_mag_limit - absolute_mag) / 5.) + 1.) * 1e-6;
    cosmology.inverse_lumdist(max_distance)
}

fn cosmic_variance(volume: f64, number_of_sightlines: f64) -> f64 {
    ((219.7 - 52.4 * volume.log10() + 3.21 * (volume.log10().powi(2)))
        / number_of_sightlines.powf(0.5))
        / 100.0
}

fn read_parquet(file_name: &str) -> PolarsResult<DataFrame> {
    let file = File::open(file_name).expect("File not found.");
    let reader = ParquetReader::new(file);
    reader.finish()
}

fn read_groups(
    file_name: &str,
    group_id_column: &str,
    redshift_column: &str,
    multi_column: &str,
    mass_column: &str,
) -> DataFrame {
    let data_frame = read_parquet(file_name).expect("Can't parse Group parquet file");
    let mut selected_df = data_frame
        .select([group_id_column, redshift_column, multi_column, mass_column])
        .expect("Column names not found");
    selected_df
        .set_column_names(vec!["group_id", "redshift", "multiplicity", "mass"])
        .unwrap();
    selected_df
}

fn read_galaxies(
    file_name: &str,
    redshift_column: &str,
    apparent_mag_colummn: &str,
    group_id_column: &str,
) -> DataFrame {
    let data_frame = read_parquet(file_name).expect("Can't parse the Galaxy data.");

    let mut selected_df = data_frame
        .select([redshift_column, apparent_mag_colummn, group_id_column])
        .unwrap();
    selected_df
        .set_column_names(vec!["redshift", "apparent_mag", "group_id"])
        .unwrap();
    selected_df
}

fn get_group_zmax(multiplicity: u32, galaxies: DataFrame) -> Result<DataFrame, PolarsError> {
    galaxies
        .clone()
        .lazy()
        .group_by([col("group_id")])
        .agg([col("zmax")
            .sort(SortOptions::default().with_order_descending(true))
            .get(
                when(col("zmax").len().gt(lit(multiplicity)))
                    .then(lit(multiplicity - 1))
                    .otherwise(col("zmax").len() - lit(1)),
            )
            .alias("group_zmax")])
        .collect()
}

fn main() {
    // Read in the group and galaxy data.
    let galaxy_file_name = "/Users/00115372/refactoring_simons_code/galaxies.parquet";
    let group_file_name = "/Users/00115372/refactoring_simons_code/groups.parquet";
    let multiplicity = 4;
    let zmin = 0.015;
    let fractional_area = 0.004361383875261386;
    let cosmo = Cosmology {
        h0: 70.,
        omega_m: 0.3,
        omega_k: 0.,
        omega_l: 0.7,
    };
    let min_volume = cosmo.comoving_volume(zmin);
    let apparent_mag_lim = 19.8;
    let mut groups = read_groups(group_file_name, "GroupID", "Zfof", "Nfof", "MassAfunc");
    let mut galaxies = read_galaxies(galaxy_file_name, "Z", "Rpetro", "GroupID");
    let redshifts: Vec<f64> = galaxies
        .column("redshift")
        .unwrap()
        .f64()
        .unwrap()
        .into_no_null_iter()
        .collect();
    let apparent_mags: Vec<f64> = galaxies
        .column("apparent_mag")
        .unwrap()
        .f64()
        .unwrap()
        .into_no_null_iter()
        .collect();
    let z_max_values: Vec<f64> = redshifts
        .iter()
        .zip(apparent_mags.iter())
        .map(|(&z, &mag)| calculate_max_redshift(z, mag, apparent_mag_lim, &cosmo))
        .collect();

    let galaxy_z_max = Series::new("zmax".into(), z_max_values);
    galaxies.with_column(galaxy_z_max).unwrap();
    let group_z_maxes = get_group_zmax(multiplicity, galaxies.clone()).unwrap();
    groups = groups
        .lazy()
        .join(
            group_z_maxes.lazy(),
            [col("group_id")],
            [col("group_id")],
            JoinArgs::new(JoinType::Inner), // Note the inner join here. Only
                                            // overlapping groups allowed.
        )
        .collect()
        .unwrap();
    let zmaxes = groups
        .column("group_zmax")
        .unwrap()
        .as_series()
        .unwrap()
        .f64()
        .unwrap();
    let max_volumes: Vec<f64> = zmaxes
        .into_iter()
        .map(|x| match x {
            Some(x) => fractional_area * (cosmo.comoving_volume(x) - min_volume),
            None => panic!("This value is not correct!!!"),
        })
        .collect();
    groups
        .with_column(Series::new("vmax".into(), max_volumes))
        .unwrap();

    //TODO: Calculate the volumes for each group.
    //TODO: Work out the cosmic variance.
    //TODO: "Monte Carlo" for the eddington bias.
    //TODO: Apply eddington bias.
    //TODO: Forward model fit.
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::iter::zip;
    use tempfile::NamedTempFile;

    #[test]
    fn test_max_redshift() {
        let limit = 19.8;
        let apparent_mags = [18.622816, 19.219088];
        let redshifts = [0.05054, 0.16029];
        let answers = [0.086145, 0.20101];
        let cosmo = Cosmology {
            h0: 67.37,
            omega_k: 0.,
            omega_m: 0.3147,
            omega_l: 1. - 0.3147,
        };
        let results: Vec<f64> = redshifts
            .iter()
            .zip(apparent_mags.iter())
            .map(|(&r, &m)| calculate_max_redshift(r, m, limit, &cosmo))
            .collect();
        for (r, a) in zip(results, answers) {
            println!("r:{r}, a:{a}");
            assert!((r - a).abs() < 1e-2)
        }
    }
    #[test]
    fn test_cosmic_variance() {
        let result = cosmic_variance(3., 3.);
        let answer = 1.128313;
        dbg!(result);
        dbg!(answer);
        assert!((result - answer).abs() < 1e-5)
    }

    #[test]
    fn test_read_groups_selects_and_renames_columns() {
        // Create input DataFrame with extra columns
        let mut df = df!(
            "z_obs" => &[0.1, 0.2],
            "group_mass" => &[1e12, 2e12],
            "junk" => &[42, 43],
            "n_members" => &[3, 5],
            "GroupID" => &[1001, 1002]
        )
        .unwrap();

        // Write to temporary parquet
        let file = NamedTempFile::new().unwrap();
        let path = file.path().to_str().unwrap();

        let mut f = File::create(path).unwrap();
        ParquetWriter::new(&mut f).finish(&mut df).unwrap();

        // Run function under test
        let out = read_groups(path, "GroupID", "z_obs", "n_members", "group_mass");

        // Assert shape
        assert_eq!(out.shape(), (2, 4));

        // Assert renamed columns
        assert_eq!(
            out.get_column_names(),
            &["group_id", "redshift", "multiplicity", "mass"]
        );

        // Assert values survived correctly
        assert_eq!(
            out.column("redshift").unwrap().f64().unwrap().get(0),
            Some(0.1)
        );
        assert_eq!(
            out.column("mass").unwrap().f64().unwrap().get(0),
            Some(1e12)
        );
        assert_eq!(
            out.column("group_id").unwrap().i32().unwrap().get(0),
            Some(1001)
        );
        assert_eq!(
            out.column("multiplicity").unwrap().i32().unwrap().get(1),
            Some(5)
        );
    }
    #[test]
    fn test_read_galaxies_selects_and_renames_columns() {
        let mut df = df!(
            "id" => &[1, 2],
            "z_obs" => &[0.1, 0.2],
            "junk" => &[2, 3],
            "mag" => &[12.2, 12.8],

        )
        .unwrap();
        let file = NamedTempFile::new().unwrap();
        let path = file.path().to_str().unwrap();
        let mut f = File::create(path).unwrap();
        ParquetWriter::new(&mut f).finish(&mut df).unwrap();

        let out = read_galaxies(path, "z_obs", "mag", "id");
        assert_eq!(out.shape(), (2, 3));
        assert_eq!(
            out.get_column_names(),
            &["redshift", "apparent_mag", "group_id"]
        );
        assert_eq!(
            out.column("redshift").unwrap().f64().unwrap().get(0),
            Some(0.1)
        );
        assert_eq!(
            out.column("apparent_mag").unwrap().f64().unwrap().get(0),
            Some(12.2)
        );
    }
    #[test]
    fn test_get_group_zmax_basic() -> Result<(), PolarsError> {
        let galaxies = df!(
            "group_id" => &[1, 1, 1, 1, 2, 2, 3],
            "zmax"     => &[0.50, 0.40, 0.30, 0.20, 0.60, 0.10, 0.25],
        )?;

        // multiplicity = 3
        let result = get_group_zmax(3, galaxies)?;

        // Sort for deterministic comparison
        let result = result.sort(["group_id"], SortMultipleOptions::default())?;

        let group_id = result.column("group_id")?.i32()?;
        let group_zmax = result.column("group_zmax")?.f64()?;

        // group 1: zmax = [0.50, 0.40, 0.30, 0.20] → 3rd largest = 0.30
        assert_eq!(group_id.get(0), Some(1));
        assert!((group_zmax.get(0).unwrap() - 0.30).abs() < 1e-12);

        // group 2: zmax = [0.60, 0.10] → smaller than multiplicity → last = 0.10
        assert_eq!(group_id.get(1), Some(2));
        assert!((group_zmax.get(1).unwrap() - 0.10).abs() < 1e-12);

        // group 3: zmax = [0.25] → single value
        assert_eq!(group_id.get(2), Some(3));
        assert!((group_zmax.get(2).unwrap() - 0.25).abs() < 1e-12);

        Ok(())
    }
}
