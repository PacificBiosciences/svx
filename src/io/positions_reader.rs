use crate::{
    core::containers::interval_tree::{Interval, IntervalTree},
    utils::util::Result,
};
use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq)]
pub struct PositionTuple {
    pub contig: String,
    pub start: i64,
    pub end: Option<i64>,
}

fn parse_position_tuple_content(content_str: &str) -> Result<PositionTuple> {
    let parts: Vec<&str> = content_str.splitn(2, ':').collect();
    if parts.len() != 2 {
        return Err(crate::svx_error!(
            "Invalid tuple content format: '{}'. Expected 'contig:start[-end]', if per contig filtering is wanted use the `--contig` flag",
            content_str
        ));
    }

    let contig_str = parts[0].trim();
    if contig_str.is_empty() {
        return Err(crate::svx_error!(
            "Invalid tuple content: Contig cannot be empty in '{}'",
            content_str
        ));
    }

    let position_part_str = parts[1].trim();
    if position_part_str.is_empty() {
        return Err(crate::svx_error!(
            "Invalid tuple content: Position part cannot be empty in '{}'",
            content_str
        ));
    }

    let start_pos: i64;
    let end_pos: Option<i64>;

    if position_part_str.contains('-') {
        let range_parts: Vec<&str> = position_part_str.splitn(2, '-').collect();
        if range_parts.len() != 2 {
            return Err(crate::svx_error!(
                "Invalid range format in '{}': Expected 'start-end' in position part '{}'",
                content_str,
                position_part_str
            ));
        }
        let start_str = range_parts[0].trim();
        let end_str = range_parts[1].trim();

        if start_str.is_empty() {
            return Err(crate::svx_error!(
                "Start part of range is empty in '{}' for position part '{}'",
                content_str,
                position_part_str
            ));
        }
        if end_str.is_empty() {
            return Err(crate::svx_error!(
                "End part of range is empty in '{}' for position part '{}'",
                content_str,
                position_part_str
            ));
        }

        let start_val = start_str.parse::<i64>().map_err(|error| {
            crate::svx_error!(
                "Failed to parse start position '{}' in tuple content '{}': {}",
                start_str,
                content_str,
                error
            )
        })?;
        if start_val < 1 {
            return Err(crate::svx_error!(
                "Start position {} must be greater than or equal to 1 (1-based) in tuple content: '{}'",
                start_val,
                content_str
            ));
        }
        start_pos = start_val - 1; // -1 for 0-based indexing

        let end_val = end_str.parse::<i64>().map_err(|error| {
            crate::svx_error!(
                "Failed to parse end position '{}' in tuple content '{}': {}",
                end_str,
                content_str,
                error
            )
        })?;
        if end_val < 1 {
            return Err(crate::svx_error!(
                "End position {} must be greater than or equal to 1 (1-based) in tuple content: '{}'",
                end_val,
                content_str
            ));
        }
        let end_val = end_val - 1; // -1 for 0-based indexing

        if start_pos > end_val {
            return Err(crate::svx_error!(
                "Start position {} (0-based) must be less than or equal to end position {} (0-based) in tuple content: '{}'",
                start_pos,
                end_val,
                content_str
            ));
        }
        end_pos = Some(end_val);
    } else {
        let start_val = position_part_str.parse::<i64>().map_err(|error| {
            crate::svx_error!(
                "Failed to parse position '{}' in tuple content '{}': {}",
                position_part_str,
                content_str,
                error
            )
        })?;
        if start_val < 1 {
            return Err(crate::svx_error!(
                "Position {} must be greater than or equal to 1 (1-based) in tuple content: '{}'",
                start_val,
                content_str
            ));
        }
        start_pos = start_val - 1; // -1 for 0-based indexing
        end_pos = None;
    }

    Ok(PositionTuple {
        contig: contig_str.to_string(),
        start: start_pos,
        end: end_pos,
    })
}

/// Parses a single position tuple in the format (contig:start) or (contig:start-end) for example: "(chr1:12345)"
pub fn parse_position_tuple(s: &str) -> Result<PositionTuple> {
    let s = s.trim();
    if !s.starts_with('(') || !s.ends_with(')') {
        return Err(crate::svx_error!(
            "Invalid format: Tuple must start with '(' and end with ')'. Got: {}",
            s
        ));
    }

    let content_str = &s[1..s.len() - 1];
    parse_position_tuple_content(content_str)
}

/// Parses a string of comma-separated position tuples in the format (contig:start) or (contig:start-end) for example: "(chr1:12345),(chr2:67890-67899)"
pub fn parse_position_tuples(s: &str) -> Result<Vec<PositionTuple>> {
    let mut positions = Vec::new();
    let s = s.trim();

    if s.is_empty() {
        return Ok(positions);
    }

    if !s.starts_with('(') || !s.ends_with(')') {
        return Err(crate::svx_error!(
            "Invalid format: String must start with '(' and end with ')'. Got: {}",
            s
        ));
    }

    let inner_s = &s[1..s.len() - 1];
    let tuple_contents: Vec<&str> = inner_s.split("),(").collect();

    for content_str in tuple_contents {
        positions.push(parse_position_tuple_content(content_str)?);
    }
    Ok(positions)
}

pub fn positions_to_interval_trees(
    positions: &Vec<PositionTuple>,
) -> Result<HashMap<String, IntervalTree<i64, ()>>> {
    let mut intervals_by_contig: HashMap<String, Vec<Interval<i64, ()>>> = HashMap::new();

    for pos_tuple in positions {
        let interval_start = pos_tuple.start;
        let interval_end = match pos_tuple.end {
            Some(e) => e + 1,            // End is inclusive, interval is [start, end+1)
            None => pos_tuple.start + 1, // Point interval [start, start+1)
        };
        let interval = Interval::new(interval_start, interval_end, ());
        intervals_by_contig
            .entry(pos_tuple.contig.clone())
            .or_default()
            .push(interval);
    }

    let interval_trees: HashMap<String, IntervalTree<i64, ()>> = intervals_by_contig
        .into_iter()
        .map(|(contig, intervals)| (contig, IntervalTree::new(intervals)))
        .collect();

    Ok(interval_trees)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::containers::interval_tree::Interval;

    #[test]
    fn test_parse_position_tuples_valid() {
        let s = "(chr1:12346),(chr2:67891-67900),(chrM:1000)";
        let expected = vec![
            PositionTuple {
                contig: "chr1".to_string(),
                start: 12345, // 12346-1
                end: None,
            },
            PositionTuple {
                contig: "chr2".to_string(),
                start: 67890,     // 67891-1
                end: Some(67899), // 67900-1
            },
            PositionTuple {
                contig: "chrM".to_string(),
                start: 999, // 1000-1
                end: None,
            },
        ];
        assert_eq!(parse_position_tuples(s).unwrap(), expected);

        let s_single = "(chrX:101)";
        let expected_single = vec![PositionTuple {
            contig: "chrX".to_string(),
            start: 100, // 101-1
            end: None,
        }];
        assert_eq!(parse_position_tuples(s_single).unwrap(), expected_single);

        let s_single_range = "(chrY:101-105)";
        let expected_single_range = vec![PositionTuple {
            contig: "chrY".to_string(),
            start: 100,     // 101-1
            end: Some(104), // 105-1
        }];
        assert_eq!(
            parse_position_tuples(s_single_range).unwrap(),
            expected_single_range
        );
    }

    #[test]
    fn test_parse_position_tuples_empty() {
        assert_eq!(parse_position_tuples("").unwrap(), vec![]);
        assert_eq!(parse_position_tuples("  ").unwrap(), vec![]);
    }

    #[test]
    fn test_parse_position_tuples_invalid() {
        assert!(parse_position_tuples("(chr1:12345").is_err()); // Missing closing paren
        assert!(parse_position_tuples("chr1:12345)").is_err()); // Missing opening paren
        assert!(parse_position_tuples("(chr1 12345)").is_err()); // Missing colon
        assert!(parse_position_tuples("(chr1:abc)").is_err()); // Invalid position
        assert!(parse_position_tuples("(chr1:123-abc)").is_err()); // Invalid end range
        assert!(parse_position_tuples("(chr1:123-45-67)").is_err()); // Too many hyphens
        assert!(parse_position_tuples("(chr1:123:45)").is_err()); // Too many colons (caught by splitn)
        assert!(parse_position_tuples("(:123)").is_err()); // Empty contig
        assert!(parse_position_tuples("(chr1:)").is_err()); // Empty position part
        assert!(parse_position_tuples("()").is_err()); // Empty tuple parts
        assert!(parse_position_tuples("(:)").is_err()); // Empty contig and position part
        assert!(parse_position_tuples("(chr1:),(chr2:2)").is_err()); // Empty position part in first tuple
        assert!(parse_position_tuples("(chr1:1),(:2)").is_err()); // Empty contig in second tuple
        assert!(parse_position_tuples("(chr1:50-10)").is_err()); // Start > end
        assert!(parse_position_tuples("(chr1:123-)").is_err()); // Empty end part of range
        assert!(parse_position_tuples("(chr1:-123)").is_err()); // Empty start part of range
        assert!(parse_position_tuples("(chr1:123),(chr1:900-1200),(chr1:9000-4200)").is_err());
        // Start > end in last entry
    }

    #[test]
    fn test_parse_position_tuples_reject_non_positive_coordinates() {
        assert!(parse_position_tuples("(chr1:0)").is_err());
        assert!(parse_position_tuples("(chr1:0-1)").is_err());
        assert!(parse_position_tuples("(chr1:0-0)").is_err());
    }

    #[test]
    fn test_positions_to_interval_trees() {
        let positions = vec![
            PositionTuple {
                // chr1:101
                contig: "chr1".to_string(),
                start: 100,
                end: None,
            },
            PositionTuple {
                // chr2:201-205 (0-based: 200-204)
                contig: "chr2".to_string(),
                start: 200,
                end: Some(204),
            },
            PositionTuple {
                // chr1:151
                contig: "chr1".to_string(),
                start: 150,
                end: None,
            },
        ];

        let trees = positions_to_interval_trees(&positions).unwrap();

        assert_eq!(trees.len(), 2);
        assert!(trees.contains_key("chr1"));
        assert!(trees.contains_key("chr2"));

        let chr1_tree = trees.get("chr1").unwrap();
        // Intervals: [100, 101) and [150, 151)
        let expected_chr1_intervals =
            vec![Interval::new(100, 101, ()), Interval::new(150, 151, ())];
        let mut found_chr1_intervals = chr1_tree.find_overlapping(0, 200); // Query range that covers both
        found_chr1_intervals.sort_by_key(|iv| iv.start);
        assert_eq!(found_chr1_intervals, expected_chr1_intervals);
        assert_eq!(chr1_tree.find_overlapping(100, 101).len(), 1); // Exact point
        assert_eq!(chr1_tree.find_overlapping(150, 151).len(), 1); // Exact point
        assert_eq!(chr1_tree.find_overlapping(120, 130).len(), 0); // No overlap

        let chr2_tree = trees.get("chr2").unwrap();
        // Interval: [200, 204+1) = [200, 205)
        let expected_chr2_intervals = vec![Interval::new(200, 205, ())];
        let found_chr2_intervals = chr2_tree.find_overlapping(0, 300); // Query range that covers it
        assert_eq!(found_chr2_intervals, expected_chr2_intervals);
        assert_eq!(chr2_tree.find_overlapping(200, 201).len(), 1); // Overlaps start
        assert_eq!(chr2_tree.find_overlapping(204, 205).len(), 1); // Overlaps end
        assert_eq!(chr2_tree.find_overlapping(199, 200).len(), 1); // Touches start (exclusive end of query)
        assert_eq!(chr2_tree.find_overlapping(205, 206).len(), 1); // Touches end (exclusive start of query)
        assert_eq!(chr2_tree.find_overlapping(200, 205).len(), 1); // Exact range
    }

    #[test]
    fn test_positions_to_interval_trees_empty() {
        let positions: Vec<PositionTuple> = vec![];
        let trees = positions_to_interval_trees(&positions).unwrap();
        assert!(trees.is_empty());
    }
}
