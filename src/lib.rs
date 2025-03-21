pub mod dir;
pub mod zoc;

pub const MAX_DEPTH: u8 = 29;

pub const fn is_valid_depth(depth: u8) -> bool {
    depth <= MAX_DEPTH
}
