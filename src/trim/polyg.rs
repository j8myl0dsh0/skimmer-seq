//! Poly-G tail removal for two-colour chemistry.
//!
//! NovaSeq and NextSeq report an unamplified base on two-colour flow cells as
//! `G`, which turns into a low-quality poly-G waterfall at the 3' end of short
//! inserts. Those `G`s score well (the caller is confident about the wrong
//! base), so a quality-only filter keeps them; they need removing by sequence.
//!
//! One non-`G` base per window is tolerated, because a single dark cycle inside
//! the run also reads as `G` and a single sequencing error reads as `A/C/T`.

/// Number of trailing bases that look like a poly-G artefact.
///
/// Stops at the first window of `window` bases that holds more than `max_other`
/// non-`G` bases, so a genuine G-rich insert tail survives.
pub fn poly_g_tail_len(seq: &[u8], window: usize, max_other: usize) -> usize {
    let window = window.max(1);
    let mut end = seq.len();
    while end > 0 {
        let start = end.saturating_sub(window);
        let others = seq[start..end]
            .iter()
            .filter(|b| !matches!(b, b'G' | b'g'))
            .count();
        if others > max_other {
            break;
        }
        end = start;
    }
    seq.len() - end
}

/// Cuts a poly-G tail from `seq` and the matching `qual` slice. Returns the
/// number of bases removed.
pub fn trim_poly_g(seq: &mut Vec<u8>, qual: &mut Vec<u8>, window: usize, max_other: usize) -> usize {
    let cut = poly_g_tail_len(seq, window, max_other);
    if cut == 0 {
        return 0;
    }
    let keep = seq.len() - cut;
    seq.truncate(keep);
    qual.truncate(keep);
    cut
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn leaves_plain_sequence_alone() {
        assert_eq!(poly_g_tail_len(b"ACGTACGT", 10, 0), 0);
    }

    #[test]
    fn strips_trailing_run() {
        assert_eq!(poly_g_tail_len(b"ACGTGGGGGG", 6, 0), 6);
    }

    #[test]
    fn tolerance_extends_the_run() {
        // One A inside the run is a dark cycle on two-colour data.
        assert_eq!(poly_g_tail_len(b"ACGTGGGAGG", 4, 1), 8);
        // Without tolerance the walk stops at the first window holding the A.
        assert_eq!(poly_g_tail_len(b"ACGTGGGAGG", 4, 0), 0);
    }

    #[test]
    fn stops_before_a_non_g_tail() {
        // Trailing Ts are real sequence: nothing to remove.
        assert_eq!(poly_g_tail_len(b"CTTTGGGGTTTT", 4, 0), 0);
    }

    #[test]
    fn trims_sequence_and_quality_together() {
        let (mut seq, mut qual) = (b"ACGTGGGG".to_vec(), b"IIIIIIII".to_vec());
        assert_eq!(trim_poly_g(&mut seq, &mut qual, 4, 0), 4);
        assert_eq!(seq, b"ACGT".to_vec());
        assert_eq!(qual.len(), 4);
    }

    #[test]
    fn whole_read_of_g_collapses_to_empty() {
        let (mut seq, mut qual) = (b"GGGGG".to_vec(), b"IIIII".to_vec());
        assert_eq!(trim_poly_g(&mut seq, &mut qual, 4, 0), 5);
        assert!(seq.is_empty());
        assert!(qual.is_empty());
    }
}
