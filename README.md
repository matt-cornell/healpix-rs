# HEALPix manipulation

A HEALPix is a way of dividing a sphere such that every point has an equal area. This library gives utilities for working with the NESTED scheme, which is optimal for neighbor checks and subdivision. A good deal of code has been taken from the [`cdshealpix`](https://crates.io/crates/cdshealpix) crate, but the API has been reworked to be cleaner and easier to use.

The majority of the implementation logic was forked from commit `5e58d17` of <https://github.com/cds-astro/cds-healpix-rust.git>.
