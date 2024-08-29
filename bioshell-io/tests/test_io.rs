#[cfg(test)]
mod tests_io {
    use std::io::BufReader;
    use bioshell_io::{open_file, read_csv, read_tsv, split_into_strings};

    #[test]
    fn test_tsv_loading() {
        let txt_f64 = "1.0\t2.0\t3.0\t4.0
5.0\t6.0\t7.0\t8.0
9.0\t10.0\t11.0\t12.0
";
        let data_f64: Vec<Vec<f64>> = read_tsv(BufReader::new(txt_f64.as_bytes())).unwrap();
        assert_eq!(data_f64.len(), 3);
        assert_eq!(data_f64[0].len(), 4);
        let txt_u16 = "1\t2\t3\t4
5\t6\t7\t8
9\t10\t11\t12";
        let data_u16: Vec<Vec<u16>> = read_tsv(BufReader::new(txt_u16.as_bytes())).unwrap();
        assert_eq!(data_u16.len(), 3);
        assert_eq!(data_u16[2].len(), 4);
    }

    #[test]
    fn test_csv_loading() {
        let reader = open_file("tests/test_files/f64.csv").expect("Can't open f64.csv file!");
        let data_f64: Vec<Vec<f64>> = read_csv(reader).expect("Can't parse f64.csv file!");
        assert_eq!(data_f64.len(), 2);
        assert_eq!(data_f64[1].len(), 3);
    }


    #[test]
    fn test_split_string() {
        let examples = [
            (r#"5 non-polymer syn '2,2',2"-[1,2,3-BENZENE-TRIYLTRIS(OXY)]TRIS[N,N,N-TRIETHYLETHANAMINIUM]' 510.816"#, 5),
            ("foo bar", 2),
            ("'foo bar'", 1),
            ("'foo bar' 'foo\" bar'", 2),
        ];

        for (string, expected) in examples {
            let splitted: Vec<String> = split_into_strings(&string, false);
            // println!("{:?}", splitted);
            assert_eq!(splitted.len(), expected);
        }
    }

}