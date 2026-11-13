//! Phred+33 decoding and 3' low-quality tail trimming.
//!
//! Trimming walks windows from the end of the read and stops at the first window
//! that clears the threshold, so the kept prefix always ends in signal. Sorting
//! is not needed and memory use is a fixed-width window counter.

/// ASCII offset of Phred score 0 in Sanger/Illumina 1.8+ quality strings.
pub const PHRED_OFFSET: u8 = 33;

/// Mean Phred score of a quality string. An empty string scores `0.0`.
pub fn mean_phred(qual: &[u8]) -> f64 {
    if qual.is_empty() {
        return 0.0;
    }
    let sum: u32 = qual
        .iter()
        .map(|&q| u32::from(q.saturating_sub(PHRED_OFFSET)))
        .sum();
    f64::from(sum) / qual.len() as f64
}

/// Cuts trailing windows whose mean score is below `min_mean`, in parallel from
/// `seq` and `qual`. Returns the kept length.
///
/// The window shrinks at the front of the read, so a read shorter than `window`
/// is judged as a whole. A `min_mean` of `0.0` never trims.
pub fn trim_low_tail(seq: &mut Vec<u8>, qual: &mut Vec<u8>, window: usize, min_mean: f64) -> usize {
    debug_assert_eq!(seq.len(), qual.len());
    let window = window.max(1);
    let mut end = qual.len();
    while end > 0 {
        let start = end.saturating_sub(window);
        if mean_phred(&qual[start..end]) >= min_mean {
            break;
        }
        end = start;
    }
    seq.truncate(end);
    qual.truncate(end);
    end
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pair(seq: &str, qual: &str) -> (Vec<u8>, Vec<u8>) {
        (seq.as_bytes().to_vec(), qual.as_bytes().to_vec())
    }

    #[test]
    fn empty_quality_scores_zero() {
        assert_eq!(mean_phred(b""), 0.0);
    }

    #[test]
    fn decodes_offset_scores() {
        assert_eq!(mean_phred(b"IIII"), 40.0);
        assert_eq!(mean_phred(b"!I"), 20.0);
    }

    #[test]
    fn clamps_below_offset() {
        // 0x1f is below '!', so it must not wrap around.
        assert_eq!(mean_phred(&[0x1f, b'I']), 20.0);
    }

    #[test]
    fn leaves_good_read_alone() {
        let (mut seq, mut qual) = pair("ACGTACGT", "IIIIIIII");
        assert_eq!(trim_low_tail(&mut seq, &mut qual, 4, 20.0), 8);
        assert_eq!(seq.len(), 8);
    }

    #[test]
    fn strips_low_tail_in_whole_windows() {
        // 8 bp at Q40 then 4 bp at Q2; window 4 removes exactly the low part.
        let (mut seq, mut qual) = pair("ACGTACGTACGT", "IIIIIIII####");
        assert_eq!(trim_low_tail(&mut seq, &mut qual, 4, 20.0), 8);
        assert_eq!(seq, b"ACGTACGT".to_vec());
        assert_eq!(qual, b"IIIIIIII".to_vec());
    }

    #[test]
    fn trims_entirely_low_read() {
        let (mut seq, mut qual) = pair("ACGTA", "#####");
        assert_eq!(trim_low_tail(&mut seq, &mut qual, 4, 20.0), 0);
        assert!(seq.is_empty());
        assert!(qual.is_empty());
    }
}
