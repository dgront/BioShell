use std::io::{BufRead, Lines};

pub(crate) struct CifLineIterator<R: BufRead> {
    inner: Lines<R>,
}

impl<R: BufRead> CifLineIterator<R> {
    pub(crate) fn new(lines: Lines<R>) -> Self {
        Self { inner: lines, }
    }
}

impl<R: BufRead> Iterator for CifLineIterator<R> {
    type Item = String;

    fn next(&mut self) -> Option<Self::Item> {
        let mut current: Option<String> = None;

        while let Some(Ok(mut line)) = self.inner.next() {
            if line.starts_with(';') {                           // --- we have a multiline marker
                if let Some(ref mut content) = current {   // --- it's a closing semicolon
                    content.push_str(&line);
                    return current;
                } else {                                        // --- no, it's the semicolon opening a multiline
                    line.push_str("\n");
                    current = Some(line);
                }
            } else {
                if let Some(ref mut content) = current {   // --- it's a continuation of a multiline
                    content.push_str(&line);
                    content.push_str("\n");
                } else {
                    return Some(line);                          // --- it's a regular line
                }
            }
        }
        current
    }
}
#[cfg(test)]
mod tests {
    use std::io::{BufReader, Cursor};
    use super::*;

    fn run_cif_line_iterator(input_str: &str) {

        let cursor = Cursor::new(input_str);
        let reader = BufReader::new(cursor);

        // Get a Lines iterator
        let lines = reader.lines();
        let mut line_iter = CifLineIterator::new(lines);
        let mut result: String = String::new();
        while let Some(line) = line_iter.next() {
            result.push_str(&line);
            result.push_str("\n");
        }
        assert_eq!(result.trim(), input_str.trim());
    }

    #[test]
    fn test_cif_line_iterator() {
        let input1: &'static str = r#";
multiline
;
"#;
        let input2: &'static str = r#";multiline
;
"#;
        let input3: &'static str = r#"line1
line2
line3
"#;

        let input4: &'static str = r#";
multiline

multiline2
;
"#;
        run_cif_line_iterator(input1);
        run_cif_line_iterator(input2);
        run_cif_line_iterator(input3);
        run_cif_line_iterator(input4);
    }

}