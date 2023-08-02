use crate::sequence::Sequence;

pub trait SequenceFilter {
    fn filter(&self, sequence: &Sequence) -> bool;
}

pub struct DescriptionContains { pub substring: String}

impl SequenceFilter for DescriptionContains {
    fn filter(&self, sequence: &Sequence) -> bool {
        sequence.description().contains(self.substring.as_str())
    }
}

pub struct AlwaysTrue;

impl SequenceFilter for AlwaysTrue {
    fn filter(&self, _: &Sequence) -> bool { true }
}
