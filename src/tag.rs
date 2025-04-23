//! # Tagged Hashes
//!
//! I've found it really easy to accidentally mix up hashes when working with varying depths. To help fix this, I've added some basic compile-time validation for these hashes.
//! In the simplest form, these tags force you to track "where" a hash or range came from.
use std::fmt::Debug;


pub trait Tag {
    fn 
}