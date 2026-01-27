use cosmoxide::Cosmology;
use polars::prelude::*;
use polars::{error::PolarsResult, frame::DataFrame, prelude::ParquetReader};
use std::fs::File;

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
    redshift_column: &str,
    multi_column: &str,
    mass_column: &str,
) -> DataFrame {
    let data_frame = read_parquet(file_name).expect("Can't parse Group parquet file");
    let mut selected_df = data_frame
        .select([redshift_column, multi_column, mass_column])
        .expect("Column names not found");
    selected_df
        .set_column_names(vec!["redshift", "multiplicity", "mass"])
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

fn main() {
    // Read in the group and galaxy data.
    let galaxy_file_name = "/Users/00115372/refactoring_simons_code/galaxies.parquet";
    let group_file_name = "/Users/00115372/refactoring_simons_code/groups.parquet";
    let cosmo = Cosmology {
        h0: 70.,
        omega_m: 0.3,
        omega_k: 0.,
        omega_l: 0.7,
    };
    let apparent_mag_lim = 19.8;
    let groups = read_groups(group_file_name, "Zfof", "Nfof", "MassAfunc");
    let galaxies = read_galaxies(galaxy_file_name, "Z", "Rpetro", "GroupID");
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
    let z_max_values = redshifts
        .iter()
        .zip(apparent_mags.iter())
        .map(|(&z, &mag)| calculate_max_redshift(z, mag, apparent_mag_lim, &cosmo));

    //let galaxy_z_max = Series::new("zmax", calculate_max_redshift(, apparent_mag, apparent_mag_lim, cosmo))

    //TODO: Determine the maximum redshift for each group.
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
        )
        .unwrap();

        // Write to temporary parquet
        let file = NamedTempFile::new().unwrap();
        let path = file.path().to_str().unwrap();

        let mut f = File::create(path).unwrap();
        ParquetWriter::new(&mut f).finish(&mut df).unwrap();

        // Run function under test
        let out = read_groups(path, "z_obs", "n_members", "group_mass");

        // Assert shape
        assert_eq!(out.shape(), (2, 3));

        // Assert renamed columns
        assert_eq!(
            out.get_column_names(),
            &["redshift", "multiplicity", "mass"]
        );

        // Assert values survived correctly
        assert_eq!(
            out.column("redshift").unwrap().f64().unwrap().get(0),
            Some(0.1)
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
}
