
/// Splits a string by a whitespace into strings that can contain a whitespace within
///
/// Unlike the [`str::split_whitespace()`](str::split_whitespace) from the Rust standard library, this function
/// takes quotation marks (both single and double) into account; a quoted string, that may contain white space
/// characters, is returned as a single token. The function also allows for nested quotes.
/// The second parameter of the function determines whether the quotation mark are removed
/// from a substring token or not.
///
/// # Examples
///
/// ```rust
/// use bioshell_io::split_into_strings;
/// let tokens_easy = split_into_strings("The quick brown fox jumps over the lazy dog", false);
/// assert_eq!(tokens_easy.len(), 9);
/// let tokens_quoted = split_into_strings("The 'quick brown fox' jumps over the 'lazy dog'", false);
/// assert_eq!(tokens_quoted.len(), 6);
/// assert_eq!(tokens_quoted[1], "'quick brown fox'");
/// let tokens_quoted = split_into_strings("The 'quick brown fox' jumps over the 'lazy dog'", true);
/// assert_eq!(tokens_quoted[1], "quick brown fox");
/// let tokens_tricky = split_into_strings("O \"O5'\" \"O5'\"", false);
/// assert_eq!(tokens_tricky.len(), 3);
/// let tokens_tricky = split_into_strings("O \"O5'\" \"O1\"", false);
/// assert_eq!(tokens_tricky.len(), 3);
/// let cif_tokens = split_into_strings("A   'RNA linking'       y \"ADENOSINE-5'-MONOPHOSPHATE\" ? 'C10 H14 N5 O7 P' 347.221", false);
/// assert_eq!(cif_tokens.len(), 7);
/// assert_eq!(cif_tokens[3], "\"ADENOSINE-5'-MONOPHOSPHATE\"".to_string());
/// ```
pub fn split_into_strings(s: &str, if_remove_quotes: bool) -> Vec<String> {

    let mut tokens: Vec<String> = Vec::new();
    let mut iter = s.split_whitespace();
    tokens.push(iter.next().unwrap().to_string());
    for word in iter {
        let token = tokens.last_mut().unwrap();
        match quote_style(token) {
            QuoteStyle::None | QuoteStyle::End(_) | QuoteStyle::Both(_) => { tokens.push(word.to_string()); }
            QuoteStyle::Begin(_) => {
                token.push_str(" ");
                token.push_str(word);
            }
        }
    }

    if if_remove_quotes { remove_paired_quotes(&mut tokens) }

    return tokens;
}

/// helps to dispatch the right quoting variant
#[derive(PartialEq, Debug)]
enum QuoteStyle {
    None,
    Begin(char),
    End(char),
    Both(char),
}
fn quote_style(token: &str) -> QuoteStyle {
    // --- a special case when a token is a quote character itself
    if token.len() == 1 && is_quote_char(token.chars().next().unwrap()) {
        return QuoteStyle::Begin(token.chars().next().unwrap());
    }

    let first_char = token.chars().next();
    let last_char = token.chars().rev().next();

    match (first_char, last_char) {
        (Some(f), Some(l)) if f == l && is_quote_char(f) => QuoteStyle::Both(f),
        (Some(f), _) if is_quote_char(f) => QuoteStyle::Begin(f),
        (_, Some(l)) if is_quote_char(l) => QuoteStyle::End(l),
        _ => QuoteStyle::None,
    }
}

/// Checks if a character is a common quotation mark
fn is_quote_char(c: char) -> bool { c == '"' || c == '\'' || c == '“' || c == '”' || c == '‘' || c == '’' }

/// Removes paired quotes from each string in a vector
fn remove_paired_quotes(strings: &mut [String]) {

    for s in strings {
        if s.len() >= 2 {
            let first_char = s.chars().next().unwrap();
            let last_char = s.chars().last().unwrap();
            // Check if the first and last characters are the same type of quote
            if (first_char == '\'' || first_char == '"') && first_char == last_char {
                // Remove the first and last characters in place
                s.remove(0); // Remove first character
                s.pop();     // Remove last character
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::split_into_strings::{quote_style, QuoteStyle, remove_paired_quotes};

    #[test]
    fn test_quote_style() {
        let examples = [
            ("\"Hello\"", QuoteStyle::Both('"')),
            ("'World", QuoteStyle::Begin('\'')),
            ("Rust'", QuoteStyle::End('\'')),
            ("No quotes", QuoteStyle::None),
            ("\"Fancy\"", QuoteStyle::Both('"')),
            ("'Smart'", QuoteStyle::Both('\'')),
            ("'Different\"", QuoteStyle::Begin('\'')),
        ];

        for (input, expected) in examples {
            assert_eq!(quote_style(input), expected, "Failed on input: {}", input);
        }
    }


    #[test]
    fn test_remove_paired_quotes() {
        let mut strings = vec![
            String::from("'hello'"),
            String::from("\"world\""),
            String::from("no_quotes"),
            String::from("'mismatched\""),
            String::from("''"),
            String::from("\"\""),
            String::from("single'"),
            String::from("\"double\""),
        ];

        let expected_results = vec![
            String::from("hello"),
            String::from("world"),
            String::from("no_quotes"),
            String::from("'mismatched\""),
            String::from(""),
            String::from(""),
            String::from("single'"),
            String::from("double"),
        ];
        remove_paired_quotes(&mut strings);

        for (i, (result, expected)) in strings.iter().zip(expected_results.iter()).enumerate() {
            assert_eq!(result, expected, "Mismatch at index {}: got '{}', expected '{}'", i, result, expected);
        }
    }
}